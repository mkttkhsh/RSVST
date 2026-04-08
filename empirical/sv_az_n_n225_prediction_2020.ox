/*
**  MCMC estimation for realized stochastic volatility model
**  with skew t-distribution by Azzalini (1985)
**  using class RSVAZ.ox
**
**  Reference:
**    Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).
**    Realized Stochastic Volatility Model with Skew-t Distributions for
**    Improved Volatility and Quantile Forecasting. arXiv:2501.13179
**    https://arxiv.org/abs/2501.13179
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<SVAZ.ox>

main()
{
	decl sdata_dir;
	decl sdata, mdat, vdate, vdr, vy;
	decl nobs, nest, nprd, vburn;
	decl dmu, dsig, dphi, drho, ddelta, dnu;
	decl dmm0, dsm0, dn0, dS0, dap0, dbp0, dar0, dbr0, dad0, dbd0, dnn0, dnS0;
	decl fls, flt, flmm, smodel, sdir;
	decl beg, end, time, timetmp, t, stemp, vytmp, vposm;

	// set data
	sdata = "N225";
	sdata_dir = "../data/";  // <- Set this to your data directory
	mdat = loadmat(sprint(sdata_dir, sdata, "_rcc_rv5_rkth2_bv_medrv_2012_2020.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vy = vdr;

	// Setup for rolling estimation and prediction
	nobs = sizer(vy);
	nprd = 242;
	nest = nobs - nprd;
	vburn = 5000 | constant(1500, nprd-1, 1); //burn-in: 5000 for the first estimation, 1500 from then on
	
	// initial values
	dphi = 0.95;
	dsig = 0.2;
	drho = -0.3;
	dmu = 0.4;
	ddelta = -0.7;
	dnu = 10;
	
	// hyper parameters for priors
	dmm0 = 0;		//mu ~ N(mm0, sm0^2)
	dsm0 = 10;	
	dn0 = 5;		//sigma^2 ~ InvGamma(n0/2, S0/2)
	dS0 = 0.05;
	dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
	dbp0 = 1.5;
	dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
	dbr0 = 2;
	dad0 = 1;		//(delta+1)/2 ~ Beta(ad0, bd0)
	dbd0 = 2;
	dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
	dnS0  = 0.5;

	// set estimated model
	fls = 0; //skewness
	flt = 0; //heavy-tail
	flmm = 1; //block or single-move sampler
	smodel = sprint("sv_az_", fls?"s":"", flt?"t":"n");
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
		vytmp =	vy[t:t+nest-1];
		print("\nSummary statistics:");
		println("%c", {"Return"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(vytmp));

		// MCMC estimation
		decl svaz = new SVAZ(vytmp);
		svaz.SetFilename(stemp);
		svaz.SetDirectory(sdir);
		svaz.SetDist(fls, flt); // (1,1):Skew t; (1,0):Skew normal; (0,1):Student's t; (0,0):Normal
		svaz.SetMultiMove(flmm); //0:single-move; 1:multi-move
		svaz.SetCount(-1); // positive: print counter, negative: noprint
		svaz.SetSampSigRho(1); //0:separate sampling; 1:joint sampling
		svaz.SetSampNu(0); //0:M-H sampling; 1:A-R M-H sampling
		svaz.SetForecast(1);	// 1:Volatility and quantile forecasts; 0:No forecast
		svaz.SetIni(dphi, dsig, drho, dmu, ddelta, dnu,1);
		svaz.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, dad0, dbd0, dnn0, dnS0, (t==beg));
		svaz.MCMC(15000,vburn[t-beg]);
		delete svaz;

		// 1-day ahead return
		println("\nOne-day-ahead return", vy[t+nest][]);

		// update initial values
		vposm = loadmat(sprint(sdir, "/posstat_", stemp, ".csv"))[][0];
		dphi   = vposm[0];
		dsig   = vposm[1];
		drho   = vposm[2];
		dmu    = vposm[3];
		ddelta = fls ? vposm[4] : 0;
		dnu    = flt ? vposm[4+fls] : M_INF_POS;

		println("\nComputation time for this period: ", timespan(timetmp));
	}

	println("\n\n=========================================");
	println("\nComputation time for entire period: ", timespan(time));
}
