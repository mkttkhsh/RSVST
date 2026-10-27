/*
**  MCMC estimation for realized stochastic volatility model
**  with skew t-distribution by Azzalini (1985)
**  using class RSVAZ.ox
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<RSVAZ.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs, nest, nprd, vburn;
	decl dmu, dxi, dsig, dsigu, dphi, drho, ddelta, dnu;
	decl dmm0, dsm0, dmx0, dsx0, dn0, dS0, dsun0, dsuS0, dap0, dbp0, dar0, dbr0, dad0, dbd0, dnn0, dnS0;
	decl fls, flt, flmm, smodel, sdir;
	decl beg, time, timetmp, t, stemp, mytmp, vposm, flresume, nburn, fh, asarg, flsmoke, end;

	// set data
	sdata = "N225";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv_200906_201909.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vrv = mdat[][2]*10000; //Realized variance calculated from 5-minute returns
	nobs = sizer(vdr~vrv);
	my = vdr ~ log(vrv);

	// Setup for rolling estimation and prediction
	nobs = sizer(my);
	nprd = 590;
	nest = nobs - nprd;
	vburn = 5000 | constant(1500, nprd-1, 1); //burn-in: 5000 for the first estimation, 1500 from then on
	
	// initial values
	dphi = 0.95;
	dsig = 0.2;
	drho = -0.3;
	dmu = 0;
	ddelta = -0.7;
	dnu = 10;
	dxi = -0.8;
	dsigu = 0.3;
	
	// hyper parameters for priors
	dmm0 = 0;		//mu ~ N(mm0, sm0^2)
	dsm0 = 10;	
	dmx0 = 0;	//xi ~ N(mx0, sx0^2)
	dsx0 = 1;
	dn0 = 5;		//sigma^2 ~ InvGamma(n0/2, S0/2)
	dS0 = 0.05;
	dsun0 = 5;	//sigma_u^2 ~ InvGamma(sun0/2, suS0/2)
	dsuS0 = 0.2;
	dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
	dbp0 = 1.5;
	dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
	dbr0 = 2;
	dad0 = 1;		//(delta+1)/2 ~ Beta(ad0, bd0)
	dbd0 = 2;
	dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
	dnS0  = 0.5;

	// set estimated model
	fls = 1; //skewness
	flt = 0; //heavy-tail
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_az_", fls?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	sdir = "out/n225_main_rsv_az_sn";

	// rolling estimation and prediction
	ranseed(20261065);
	beg = 295;
	end = 354;
	flresume = 0;
	asarg = arglist(); flsmoke = sizeof(asarg) > 1 && asarg[1] == "smoke";
	if(flsmoke){ end = beg + 2; sdir = sprint(sdir, "_smoke"); }
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
		// resume: a window whose h1_y already exists is skipped; its posterior means become the initial values of the next window
		fh = fopen(sprint(sdir, "/h1_y_", stemp, ".csv"), "r");
		if(isfile(fh)){
			fclose(fh);
			println("already done, skipped");
			vposm = loadmat(sprint(sdir, "/posstat_", stemp, ".csv"))[][0];
			dphi   = vposm[0];
			dsig   = vposm[1];
			drho   = vposm[2];
			dmu    = vposm[3];
			ddelta = fls ? vposm[4] : 0;
			dnu    = flt ? vposm[4+fls] : M_INF_POS;
			dxi    = vposm[4+fls+flt];
			dsigu  = vposm[5+fls+flt];
			flresume = 1;
			continue;
		}
		nburn = (t == beg || flresume) ? 5000 : 1500;  // 5000 at the start of a chunk or right after a resume, 1500 otherwise (as in the original scripts)
		if(flsmoke) nburn = (nburn == 5000) ? 100 : 50;
		flresume = 0;

		// set data	and show summary statistics
		mytmp =	my[t:t+nest-1][];
		print("\nSummary statistics:");
		println("%c", {"Return", "log-RK"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(mytmp));

		// MCMC estimation
		decl rsvaz = new RSVAZ(mytmp);
		rsvaz.SetFilename(stemp);
		rsvaz.SetDirectory(sdir);
		rsvaz.SetDist(fls, flt); // (1,1):Skew t; (1,0):Skew normal; (0,1):Student's t; (0,0):Normal
		rsvaz.SetMultiMove(flmm); //0:single-move; 1:multi-move
		rsvaz.SetCount(-1); // positive: print counter, negative: noprint
		rsvaz.SetSampSigRho(1); //0:separate sampling; 1:joint sampling
		rsvaz.SetSampNu(0); //0:M-H sampling; 1:A-R M-H sampling
		rsvaz.SetForecast(1);	// 1:Volatility and quantile forecasts; 0:No forecast
		rsvaz.SetIni(dphi, dsig, drho, dmu, ddelta, dnu, dxi, dsigu, 1);
		rsvaz.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, dad0, dbd0, dnn0, dnS0, dmx0, dsx0, dsun0, dsuS0, (t==beg));
		rsvaz.MCMC(flsmoke ? 300 : 15000, nburn);
		delete rsvaz;

		// 1-day ahead return and RK
		println("\nOne-day-ahead return, log-RK, and RK", my[t+nest][]~exp(my[t+nest][1]));

		// update initial values
		vposm = loadmat(sprint(sdir, "/posstat_", stemp, ".csv"))[][0];
		dphi   = vposm[0];
		dsig   = vposm[1];
		drho   = vposm[2];
		dmu    = vposm[3];
		ddelta = fls ? vposm[4] : 0;
		dnu    = flt ? vposm[4+fls] : M_INF_POS;
		dxi    = vposm[4+fls+flt];
		dsigu  = vposm[5+fls+flt];

		println("\nComputation time for this period: ", timespan(timetmp));
	}

	println("\n\n=========================================");
	println("\nComputation time for entire period: ", timespan(time));
}
