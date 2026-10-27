/*
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
**  using class RSVGH.ox
**
**  2026 revision: single full-sample in-sample psi-free estimation
**  (main-sample period, not rolling) -- reports the posterior distribution
**  / 95% CI of psi for the sensitivity analysis (Comment 2). Distinct from
**  the 590-window rolling *_prediction.ox script, which is not re-run this
**  revision; this is a single MCMC run over the whole main sample instead.
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#include<RSVGH.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs;
	decl flt, flb, flbc, flmm, smodel, sdir;
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu;
	decl dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, vnu0;
	decl time, mytmp;

	// set data
	sdata = "N225";
	mdat = loadmat(sprint("../Data/", sdata, "_rcc_rv_200906_201909.csv")); //Date, daily return, RV
	vdate = mdat[][0]; //dates yyyymmdd
	vdr = mdat[][1]*100; //close-to-close (log) return (%)
	vrv = mdat[][2]*10000; //Realized variance calculated from 5-minute returns
	nobs = sizer(vdr~vrv);
	my = vdr ~ log(vrv);

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
	dmx0 = 0;	//xi ~ N(mx0, sx0^2)
	dsx0 = 1;
	dn0 = 5;		//sigma_eta^2 ~ InvGamma(n0/2, S0/2)
	dS0 = 0.05;
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
	flbc = 1; //slope coefficient for log-RV equation (free psi: sensitivity/robustness analysis)
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_gh_", flb?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	if(flbc) smodel = sprint(smodel, "_psi");
	smodel = sprint(smodel, "_fullsample");
	sdir = sprint(sdata, "_fullsample_psi");

	// single full-sample in-sample estimation (no rolling window)
	mytmp = my;
	print("\nData period: ");
	print("%8.0f", vdate[0]);
	print(" - ");
	print("%8.0f", vdate[nobs-1]);
	print("\n");
	print("\nSummary statistics:");
	println("%c", {"Return", "log-RV"}, "%r", {"NOBS", "Mean", "SD", "Skew", "Kurt"}, moments(mytmp));

	time = timer();
	decl rsvgh = new RSVGH(mytmp);	//input data (0:simulation)
	rsvgh.SetFilename(smodel);
	rsvgh.SetDirectory(sdir);
	rsvgh.SetDist(flt + flb - 1); // -1:normal, 0:t, 1:skew-t
	rsvgh.SetCount(-1);
	rsvgh.SetBC(flbc);	// 1: both (free psi), 0: constant only
	rsvgh.SetMultiMove(flmm); //0:single-move; 1:multi-move
	rsvgh.SetForecast(1); // 1:Volatility and quantile forecasts; 0:No forecast (simulated, no held-out data needed)
	rsvgh.SetIni(dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu, 0);
	rsvgh.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, 1);
	rsvgh.MCMC(40000, 10000);
	delete rsvgh;

	println("\n\n=========================================");
	println("\nComputation time: ", timespan(time));
}
