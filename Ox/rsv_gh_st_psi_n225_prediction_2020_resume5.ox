/*
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
**  using class RSVGH.ox
**
**  2026 revision: resume script for the N225 COVID-period psi-free GH-ST
**  rolling run, which diverged at window t=5 (2020-01-14) -- mu jumped
**  0.25->6.7 and nu pinned at exactly 4.0000 with near-zero posterior
**  variance, a clear MCMC failure, then propagated forward through every
**  subsequent window via the usual carry-forward initial-value mechanism
**  (windows t=0..4 were fine). No ranseed() was ever called in the
**  original script, so Ox's default (deterministic) seed was used --
**  a plain re-run of the whole thing would reproduce the identical
**  failure. This script instead:
**    - starts from t=5 instead of t=0 (windows 0-4 don't need re-running)
**    - seeds explicitly with a different value, to get a different draw
**      sequence at the point of failure
**    - uses window 4's posterior means (read from the existing
**      posstat_rsv_gh_st_psi_4.csv, produced by the original run) as the
**      initial values, exactly matching what the original loop would have
**      passed forward into t=5
**  Output overwrites the existing (corrupted) var_/es_/volf_/posstat_/etc.
**  files for t=5..241 in place -- t=0..4 files are untouched.
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<RSVGH.ox>

main()
{
	ranseed(20260721); // 2026 revision: explicit distinct seed for the resume, see header comment

	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs, nest, nprd, vburn;
	decl flt, flb, flbc, flmm, smodel, sdir;
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu;
	decl dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, vnu0;
	decl beg, end, time, timetmp, t, stemp, mytmp, vposm;

	// set data
	sdata = "N225";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv5_rkth2_bv_medrv_2012_2020.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vrv = mdat[][2]*10000; //Realized variance calculated from 5-minute returns
	nobs = sizer(vdr~vrv);
	my = vdr ~ log(vrv);

	// Setup for rolling estimation and prediction
	nobs = sizer(my);
	nprd = 242;
	nest = nobs - nprd;
	vburn = 5000 | constant(1500, nprd-1, 1); //burn-in: 5000 for the first estimation, 1500 from then on

	// initial values: window 4's posterior means (last clean window before the
	// t=5 divergence), read from posstat_rsv_gh_st_psi_4.csv Var1 column.
	// Index mapping matches the original script's "update initial values" block.
	dphi    = 0.9263776986448594;
	dsigeta = 0.34620141618203953;
	drho    = -0.3222246710956261;
	dmu     = 0.2528427780636074;
	dbeta   = -0.1701983115900296;
	dnu     = 16.901759457951556;
	dxi     = -1.0535157055140305;
	dpsi    = 0.8114206168005168;
	dsigu   = 0.436160778860555;

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
	flbc = 1; //slope coefficient for log-RV equation
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_gh_", flb?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	if(flbc) smodel = sprint(smodel, "_psi");
	sdir = sprint(sdata, "_prediction_2020");

	// rolling estimation and prediction -- resume from t=5
	beg = 5;
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
		rsvgh.SetBC(flbc);	// 1: both, 0: constant only
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
