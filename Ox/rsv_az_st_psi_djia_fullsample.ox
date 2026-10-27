/*
**  MCMC estimation for realized stochastic volatility model
**  with skew t-distribution by Azzalini (1985)
**  using class RSVAZ.ox
**
**  IJF R2 (2026): single full-sample in-sample psi-free estimation
**  (main-sample period, not rolling) -- reports the posterior distribution
**  / 95% CI of psi for the sensitivity analysis (Comment 2). Distinct from
**  the 590-window rolling *_prediction.ox script, which is not re-run this
**  revision; this is a single MCMC run over the whole main sample instead.
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<RSVAZ.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs;
	decl dmu, dxi, dpsi, dsig, dsigu, dphi, drho, ddelta, dnu;
	decl dmm0, dsm0, dmx0, dsx0, dn0, dS0, dsun0, dsuS0, dap0, dbp0, dar0, dbr0, dad0, dbd0, dnn0, dnS0, dmp0, dsp0;
	decl fls, flt, flbc, flmm, smodel, sdir;
	decl time, mytmp;

	// set data
	sdata = "DJIA";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv_200906_201909.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vrv = mdat[][2]*10000; //Realized variance calculated from 5-minute returns
	nobs = sizer(vdr~vrv);
	my = vdr ~ log(vrv);

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
	flt = 1; //heavy-tail
	flbc = 1; //slope coefficient psi of log-RV equation (free: sensitivity/robustness)
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_az_", fls?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	if(flbc) smodel = sprint(smodel, "_psi");
	smodel = sprint(smodel, "_fullsample");
	dpsi = 1;  //initial value for psi
	dmp0 = 1; dsp0 = 1; //psi ~ N(mp0, sp0^2)
	sdir = sprint(sdata, "_fullsample_psi");

	// single full-sample in-sample estimation (no rolling window)
	mytmp = my;
	print("\nData period: ");
	print("%8.0f", vdate[0]);
	print(" - ");
	print("%8.0f", vdate[nobs-1]);
	print("\n");
	print("\nSummary statistics:");
	println("%c", {"Return", "log-RK"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(mytmp));

	time = timer();
	decl rsvaz = new RSVAZ(mytmp);
	rsvaz.SetFilename(smodel);
	rsvaz.SetDirectory(sdir);
	rsvaz.SetDist(fls, flt); // (1,1):Skew t; (1,0):Skew normal; (0,1):Student's t; (0,0):Normal
	rsvaz.SetBC(flbc); // 1: free psi; 0: psi fixed at 1
	rsvaz.SetMultiMove(flmm); //0:single-move; 1:multi-move
	rsvaz.SetCount(-1); // positive: print counter, negative: noprint
	rsvaz.SetSampSigRho(1); //0:separate sampling; 1:joint sampling
	rsvaz.SetSampNu(0); //0:M-H sampling; 1:A-R M-H sampling
	rsvaz.SetForecast(1);	// 1:Volatility and quantile forecasts; 0:No forecast (simulated, no held-out data needed)
	rsvaz.SetIni(dphi, dsig, drho, dmu, ddelta, dnu, dxi, dsigu, 1, dpsi);
	rsvaz.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, dad0, dbd0, dnn0, dnS0, dmx0, dsx0, dsun0, dsuS0, 1, dmp0, dsp0);
	rsvaz.MCMC(40000, 10000);
	delete rsvaz;

	println("\n\n=========================================");
	println("\nComputation time: ", timespan(time));
}
