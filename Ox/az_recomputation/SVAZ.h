/*----------------------------------------------------------
                        SVAZ.h
----------------------------------------------------------*/
											 
/*
**  HEADER FILE
**
**  MCMC estimation for stochastic volatility model
**  with skew t-distribution by Azzalini (1985)
*/

class SVAZ
{
	decl m_iRun, m_vy, m_vyorg, m_ns, m_fls, m_flt, m_flsigrho, m_flnuarmh, m_flmultimove, m_sfilename, m_sdirectory;
	decl m_nKm, m_vh_t, m_vlambda_t, m_vz0_t, m_vini, m_ncount, m_k, m_vhini, m_vpri;
	decl m_dap0, m_dbp0, m_dn0, m_dS0, m_dar0, m_dbr0;		// hyper parameters for phi, sigma, and rho
	decl m_dmm0, m_dsm0, m_dad0, m_dbd0, m_dnn0, m_dnS0;	// hyper parameters for mu, delta, and nu
	decl m_cmp, m_cms, m_cmr, m_cmd, m_cmn, m_can, m_cml, m_cah, m_cmh;
	decl m_dphi, m_dsig, m_drho, m_dmu, m_ddelta, m_dnu;
	decl m_vpar, m_vz0, m_vh, m_vlambda, m_vy1, m_vhbd;
	decl m_vposm, m_flML, m_kk, m_mmls, m_flf, m_flmode;
 
	SVAZ(const mdata);	//constructer

	SetDist(const fls, const flt);
	SetFilename(const sfilename);
	SetDirectory(const sdirectory);
	SetRanseed(const irs);
	SetSampSigRho(const flsigrho);
	SetSampNu(const flnuarmh);
	SetMultiMove(const flmultimove);
	SetBlocks(const nKm);
	SetCount(const ncount);
	SetForecast(const iForecast);
	
	SimObs(const ns, const dphi, const dsig, const drho, const dmu, const ddelta, const dnu); //simulate sample
	SetIni(const dphi, const dsig, const drho, const dmu, const ddelta, const dnu, const vh); //set initial
	SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0,
		   const dad0, const dbd0, const dnn0, const dnS0, const flprint); // set hyper parameters
 	DrawData(); //draw data (& state)

	MCMC(const nsim, const nburn);
	SampMu(const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	SampPhi(const dmu, const dphio, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	SampSigma(const dmu, const dphi, const dsigo, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	fLogPosRho(const dp, const dfunc, const dsc, const dhe);
	SampRho(const dmu, const dphi, const dsig, const drhoo, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	fLogPosSigRho(const dp, const dfunc, const dsc, const dhe);
	SampSigRho(const dmu, const dphi, const dsigo, const drhoo, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	fLogPosDelta(const dp, const dfunc, const dsc, const dhe);
	SampDelta(const dmu, const dphi, const dsig, const drho, const ddeltao, const dnu, const vz0, const vh, const vlambda, const vy);
	fPosNut(const vp, const dfunc, const dsc, const dhe);
	fPosNu(const vp, const dfunc, const dsc, const dhe);
	SampNu(const dmu, const dphi, const dsig, const drho, const ddeltao, const dnu, const vz0, const vh, const vlambda, const vy);
	SampZ0(const dmu, const dphi, const dsig, const drho, const ddeltao, const dnu, const vh, const vlambda, const vy);
	SampLambda(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambdao, const vy);
	SingleMove(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);
	MultiMove(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy);

	fLogPri(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu);

	ML(const nI, const nR, const nsim);
	PFilter(const nI, const nR);

	PPA(const nsim);
	
	fTsvar(const vx, const iBm);
	fGeweke(const vx, const iBm);
	
	fBetaMoment(const da, const db, const cm);
	fGammaMoment(const da, const db, const cm);
	fInvGammaMoment(const da, const db, const cm);
};
