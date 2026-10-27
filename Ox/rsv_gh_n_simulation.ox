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

	// Setup for simulation study
	nobs = 2000;

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
	flt = 0; //heavy-tail
	flb = 0; //skewness
	flbc = 0; //slope coefficient for log-RV equation
	flmm = 1; //block or single-move sampler
	smodel = sprint("rsv_gh_", flb?"s":"", flt?"t":"n");
	if(!flmm) smodel = sprint(smodel, "_sm");
	sdir = "simulation";

	timetmp = timer();

    // MCMC estimation
    decl rsvgh = new RSVGH(0);	//input data (0:simulation)
    rsvgh.SetFilename(smodel);
    rsvgh.SetDirectory(sdir);
    rsvgh.SetDist(flt + flb - 1); // -1:normal, 0:t, 1:skew-t
    rsvgh.SimObs(nobs, dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu);
	rsvgh.DrawData();
    rsvgh.SetCount(5000);
    rsvgh.SetBC(0);	// 1: both, 0: constant only
    rsvgh.SetMultiMove(flmm); //0:single-move; 1:multi-move
    rsvgh.SetForecast(0); // 1:Volatility and quantile forecasts; 0:No forecast
    rsvgh.SetIni(dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu, 1);
    rsvgh.SetPri(dap0, dbp0, dn0, dS0, dar0, dbr0, dmm0, dsm0, db0, dsb02, dnn0, dnS0, dmx0, dsx0, dmp0, dsp0, dsun0, dsuS0, 1);
    rsvgh.MCMC(50000,10000);
    delete rsvgh;

	println("\nComputation time: ", timespan(timetmp));
}
