/*
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
**  using class SVGH.ox
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<SVGH.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vy;
	decl nobs, nest, nprd, vburn;
	decl flt, flb, smodel, sdir;
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu;
	decl dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0;
	decl beg, end, time, timetmp, t, stemp, vytmp, vposm;

	// set data
	sdata = "N225";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv5_rkth2_bv_medrv_2012_2020.csv"));
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vy = vdr;

	// Setup for rolling estimation and prediction
	nobs = sizer(vy);
	nprd = 242;
	nest = nobs - nprd;
	vburn = 5000 | constant(1500, nprd-1, 1); //burn-in: 5000 for the first estimation, 1500 from then on

	// initial values
	dphi    = 0.95;
	dsigeta = 0.2;
	drho    = -0.3;
	dmu     = 0.4;
	dbeta   = 0;
	dnu     = 10;
	
	// hyper parameters for priors
	dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
	dbp0 = 1.5;
	dn0  = 5;		//sigma_eta^{-2} ~ Gamma(n0, S0)
	dS0  = 0.05;
	dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
	dbr0 = 2;
	dmm0 = 0;		//mu ~ N(mm0, sm0^2)
	dsm0 = 10;        // 2026-09-16: mu ~ N(0,100) as in the paper and all other drivers (was sqrt(10), PROGRESS §54)	
	db0   = 0;		//beta ~ N(b0, sigb02)
	dsb02 = 1;
	dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
	dnS0  = 0.5;
	
	// set estimated model
	flt = 1; //heavy-tail
	flb = 1; //skewness
	smodel = sprint("sv_gh_", flb?"s":"", flt?"t":"n");
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
		vytmp =	vy[t:t+nest-1];
		print("\nSummary statistics:");
		println("%c", {"Return"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(vytmp));
		
		// MCMC estimation
		decl svgh = new SVGH(vytmp);	//input data (0:simulation)
		svgh.SetFilename(stemp);
		svgh.SetDirectory(sdir);
		svgh.SetDist(flt + flb - 1); // -1:normal, 0:t, 1:skew-t
		svgh.SetCount(-1);
		svgh.SetForecast(1); // 1:Volatility and quantile forecasts; 0:No forecast
		svgh.SetIni(dphi, dsigeta, drho, dmu, dbeta, dnu, 0);
		svgh.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, (t==beg));
		svgh.MCMC(15000,vburn[t-beg]);
		delete svgh;

		// 1-day ahead return and RK
		println("\nOne-day-ahead return", vy[t+nest][]);

		// update initial values
		vposm = loadmat(sprint(sdir, "/posstat_", stemp, ".csv"))[][0];
		dphi    = vposm[0];
		dsigeta = vposm[1];
		drho    = vposm[2];
		dmu     = vposm[3];
		dbeta   = flb ? vposm[4] : 0;
		dnu     = flt ? vposm[4+flb] : M_INF_POS;

		println("\nComputation time for this period: ", timespan(timetmp));
	}

	println("\n\n=========================================");
	println("\nComputation time for entire period: ", timespan(time));
}
