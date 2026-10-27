/*
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
**  using class RSVGH.ox
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<RSVGH.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs, nest, nprd, vburn;
	decl flt, flb, flbc, flmm, smodel, sdir;
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu;
	decl dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, vnu0;
	decl beg, end, time, timetmp, t, stemp, mytmp, vposm;

	// set data
	sdata = "DJIA";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vrv = mdat[][2]*10000; //Realized variance calculated from 5-minute returns
	nobs = sizer(vdr~vrv);
	my = vdr ~ log(vrv);

	// Setup for rolling estimation and prediction
	nobs = sizer(my);
	nprd = 253;  // NYSE 2020暦年の取引日数（N225は242）
	nest = nobs - nprd;
	vburn = 5000 | constant(1500, nprd-1, 1); //burn-in: 5000 for the first estimation, 1500 from then on

	// initial values
	dphi    = 0.95;
	dsigeta = 0.2;
	drho    = -0.3;
	dmu     = 0;
	dbeta   = -0.7;
	dnu     = 10;
	dxi     = -0.8;
	dpsi    = 1;
	dsigu   = 0.3;
	
	// hyper parameters for priors
	dmm0 = 0;		//mu ~ N(mm0, sm0^2)
	dsm0 = 10;	
	dmx0  = 0;	//xi ~ N(mx0, sx0^2)
	dsx0  = 1;
	dn0  = 5;		//sigma_eta^{-2} ~ Gamma(n0, S0)
	dS0  = 0.05;
	dsun0	= 5;	//sigma_u^{-2} ~ Gamma(sun0, suS0)
	dsuS0 = 0.2;
	dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
	dbp0 = 1.5;
	dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
	dbr0 = 2;
	db0   = 0;		//beta ~ N(b0, sigb02)
	dsb02 = 1;
	dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
	dnS0  = 0.5;
	dmp0 = 1;	//psi ~ N(mp0, svp0^2)
	dsp0 = 1;
	
	// set estimated model
	flt = 1; //heavy-tail
	flb = 1; //skewness
	flbc = 0; //slope coefficient for log-RV equation
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_gh_", flb?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	sdir = sprint(sdata, "_prediction_2020");

	// rolling estimation and prediction
	beg = 0;
	end = nprd;
	time = timer();
	for(t=beg; t<end; t++)
	{
		// set dates
		timetmp = timer();
		println("\n\n=========================================");
		println("Iteration: ", t);
		print("Data period: ");
		print("%8.0f", vdate[t]);
		print(" - ");
		print("%8.0f", vdate[t+nest-1]);
		print("\n");
		stemp = sprint(smodel,"_",t);
		println(stemp);

		// set data	and show summary statistics
		mytmp =	my[t:t+nest-1][];
		print("\nSummary statistics:");
		println("%c", {"Return", "log-RV"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(mytmp));
		
		// MCMC estimation
		decl rsvgh = new RSVGH(mytmp);	//input data (0:simulation)
		rsvgh.SetFilename(stemp);
		rsvgh.SetDirectory(sdir);
		rsvgh.SetDist(flt + flb - 1); // -1:normal, 0:t, 1:skew-t
		rsvgh.SetCount(-1);
		rsvgh.SetBC(0);	// 1: both, 0: constant only
		rsvgh.SetMultiMove(flmm); //0:single-move; 1:multi-move
		rsvgh.SetForecast(1); // 1:Volatility and quantile forecasts; 0:No forecast
		rsvgh.SetIni(dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu, 0);
		rsvgh.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, (t==beg));
		rsvgh.MCMC(15000,vburn[t-beg]);
		delete rsvgh;

		// 1-day ahead return and RK
		println("\nOne-day-ahead return, log-RV, and RV", my[t+nest][]~exp(my[t+nest][1]));

		// update initial values
		vposm = loadmat(sprint(sdir, "/posstat_", stemp, ".csv"))[][0];
		dphi    = vposm[0];
		dsigeta = vposm[1];
		drho    = vposm[2];
		dmu     = vposm[3];
		dbeta   = flb ? vposm[4] : 0;
		dnu     = flt ? vposm[4+flb] : M_INF_POS;
		dxi     = vposm[4+flt+flb];
		dpsi    = flbc ? vposm[4+flt+flb+1] : 1;
		dsigu   = vposm[5+flt+flb+flbc];

		println("\nComputation time for this period: ", timespan(timetmp));
	}

	println("\n\n=========================================");
	println("\nComputation time for entire period: ", timespan(time));
}
