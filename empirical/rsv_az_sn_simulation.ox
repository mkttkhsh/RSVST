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
#include<RSVAZ.ox>

main()
{
	decl sdata, mdat, vdate, vdr, vrv, vrk, my;
	decl nobs, nest, nprd, vburn;
	decl dmu, dxi, dsig, dsigu, dphi, drho, ddelta, dnu;
	decl dmm0, dsm0, dmx0, dsx0, dn0, dS0, dsun0, dsuS0, dap0, dbp0, dar0, dbr0, dad0, dbd0, dnn0, dnS0;
	decl fls, flt, flmm, smodel, sdir;
	decl beg, end, time, timetmp, t, stemp, mytmp, vposm;

	// Setup for simulation study
	nobs = 2000;
	
	// initial values
	dphi = 0.9;
	dsig = 0.3;
	drho = -0.4;
	dmu = -0.5;
	ddelta = -0.9;
	dnu = 10;
	dxi = -0.5;
	dsigu = 0.2;
	
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
	sdir = "simulation";

	timetmp = timer();

    // MCMC estimation
    decl rsvaz = new RSVAZ(0); // simulation study
    rsvaz.SetRanseed(13);
    rsvaz.SetFilename(smodel);
    rsvaz.SetDirectory(sdir);
    rsvaz.SetDist(fls, flt); // (1,1):Skew t; (1,0):Skew normal; (0,1):Student's t; (0,0):Normal
    rsvaz.SimObs(nobs, dphi, dsig, drho, dmu, ddelta, dnu, dxi, dsigu);
	rsvaz.DrawData();
    rsvaz.SetMultiMove(flmm); //0:single-move; 1:multi-move
    rsvaz.SetCount(5000); // positive: print counter, negative: noprint
    rsvaz.SetSampSigRho(1); //0:separate sampling; 1:joint sampling
    rsvaz.SetSampNu(0); //0:M-H sampling; 1:A-R M-H sampling
    rsvaz.SetForecast(0);	// 1:Volatility and quantile forecasts; 0:No forecast
    rsvaz.SetIni(dphi, dsig, drho, dmu, ddelta, dnu, dxi, dsigu, 1);
    rsvaz.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, dad0, dbd0, dnn0, dnS0, dmx0, dsx0, dsun0, dsuS0, 1);
    rsvaz.MCMC(50000,10000);
    delete rsvaz;

	println("\nComputation time: ", timespan(timetmp));
}
