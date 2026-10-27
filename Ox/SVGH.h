/*----------------------------------------------------------
                        SVGH.h
----------------------------------------------------------*/
											 
/*
**  HEADER FILE
**
**  MCMC estimation for stochastic volatility model
**  with generalized hyperbolic distribution
*/

class SVGH
{
	decl m_iRun, m_vy, m_vyorg, m_ns, m_flt, m_flb, m_sdirectory, m_sfilename;
	decl m_nKm, m_vh_t, m_vz_t, m_vini, m_ncount, m_k, m_vhini, m_vpri;
	decl m_dap0, m_dbp0, m_dn0, m_dS0, m_dar0, m_dbr0;
	decl m_dmm0, m_dsm0, m_db0, m_dsb02, m_dnn0, m_dnS0;
	decl m_cmp, m_cms, m_cmn, m_cmz, m_cmb, m_cmh, m_cah;
	decl m_cnu2;	//counter of breaks in loop for sampling nu
	decl m_dphi, m_dsig, m_drho, m_dmu, m_dbeta, m_dnu;
	decl m_vh, m_vhbd, m_vz, m_vyb0;
	decl m_vposm, m_flML, m_kk, m_mmls, m_flf;
	decl m_ni;
 
	SVGH(const vdata);	//constructer

	SetDirectory(const sdirectory);
	SetDist(const iDist);
	SetFilename(const sfilename);
	SetRanseed(const irs);
	SetBlocks(const nKm);
	SetCount(const ncount);
	SetForecast(const iForecast);
	
	SimObs(const ns, const dphi, const dsig, const drho,
		   const dmu, const dbeta_a, const dnu);
			 			//simulate sample
	SetIni(const dphi, const dsig, const drho,
		   const dmu, const dbeta, const dnu, const vh);
		   				//set initial
	SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0,
		   const db0, const dsb02, const dnn0, const dnS0, const flprint);	// set hyper parameters
 	DrawData(); //draw data (& state)

	MCMC(const nsim, const nburn);
	SampPhi(const vyb, const dphio, const dsig,
			const drho, const vhb);
	fPosSigRho(const vp, const dfunc, const vsc,
			   const mhe);
	SampSigRho(const dphi, const dsigo, const drhoo,
			   const vh, const vyb0);
	SampMu(const dphi, const dsig, const drho, const dbeta,
		   const vyb, const vh);
	fPosBeta(const dp, const dfunc, const dsc, const dhe);
	SampBeta(const dphi, const dsig, const drho, const dmu,	const dbetao, const dnu,
			 const vh, const vz, const vzb);
	fPosNu(const vp, const dfunc, const dsc, const dhe);
	SampNu(const dphi, const dsig, const drho, const dmu,
		   const dbeta, const dnuo, const vh, const vz);
	SampZ(const dphi, const dsig, const drho, const dbeta,
		  const dlam, const ddel, const vh, const vzo);
	MultiMove(const dphi, const dsig, const drho,
			  const dgam, const dnu, const dbeta, const vh, const vz,
			  const vzb);

	ML(const nI, const nR, const nsim);
	PFilter(const nI, const nR);

	fTsvar(const vx, const iBm);
	fGeweke(const vx, const iBm);
	
	fBetaMoment(const da, const db, const cm);
	fGammaMoment(const da, const db, const cm);
	fInvGammaMoment(const da, const db, const cm);
};
