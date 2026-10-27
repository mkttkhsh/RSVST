/*----------------------------------------------------------
                        RSVGH.h
----------------------------------------------------------*/
											 
/*
**  HEADER FILE
**
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
*/

class RSVGH
{
	decl m_iRun, m_my, m_myorg, m_ns, m_flt, m_flb, m_flbc, m_flmultimove, m_sfilename, m_sdirectory;
	decl m_nKm, m_vh_t, m_vz_t, m_vini, m_ncount, m_k, m_vhini, m_vpri;
	decl m_dap0, m_dbp0, m_dn0, m_dS0, m_dar0, m_dbr0;		// hyper parameters for phi, sigma_eta, and rho
	decl m_dmm0, m_dsm0, m_db0, m_dsb02, m_dnn0, m_dnS0;	// hyper parameters for mu, beta, and nu
	decl m_dmx0, m_dsx0, m_dmp0, m_dsp0, m_dsun0, m_dsuS0;	// hyper parameters for xi, varphi, and sigma_u
	decl m_cmp, m_cms, m_cmn, m_cmb, m_cmz, m_cmh, m_cah;
	decl m_dphi, m_dsigeta, m_drho, m_dmu, m_dbeta, m_dnu;
	decl m_vh, m_vhbd, m_vz, m_vy1b0;
	decl m_vposm, m_flML, m_kk, m_mmls, m_flf;
	decl m_cnu1, m_cnu2;
 
	RSVGH(const mdata);	//constructer

	SetDist(const iDist);
	SetBC(const iBC);
	SetMultiMove(const flmultimove);
	SetFilename(const sfilename);
	SetDirectory(const sdirectory);
	SetRanseed(const irs);
	SetBlocks(const nKm);
	SetCount(const ncount);
	SetForecast(const iForecast);
	
	SimObs(const ns, const dphi, const dsigeta, const drho,
		   const dmu, const dbeta_a, const dnu, const dxi, const dpsi, const dsigu);
			 			//simulate sample
	SetIni(const dphi, const dsigeta, const drho,
		   const dmu, const dbeta, const dnu, const dxi, const dpsi, const dsigu, const vh);
		   				//set initial
	SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0,
		   const db0, const dsb02, const dnn0, const dnS0, const dmx0, const dsx0, const dmp0, const dsp0,
		   const dsun0, const dsuS0, const flprint);
						// set hyper parameters
 	DrawData(); //draw data (& state)

	MCMC(const nsim, const nburn);
	SampPhi(const vyb, const dphio, const dsigeta,
			const drho, const vhb);
	fPosSigRho(const vp, const dfunc, const vsc,
			   const mhe);
	SampSigRho(const dphi, const dsigetao, const drhoo,
			   const vh, const vyb0);
	SampMu(const dphi, const dsigeta, const drho, const dbeta,
		   const vyb, const vh);
	fPosBeta(const dp, const dfunc, const dsc, const dhe);
	SampBeta(const dphi, const dsigeta, const drho, const dmu, const dbetao, const dnu,
			 const vh, const vz, const vzb);
	fPosNut(const vp, const dfunc, const dsc, const dhe);
	fPosNu(const vp, const dfunc, const dsc, const dhe);
	SampNu(const dphi, const dsigeta, const drho, const dmu,
		   const dbeta, const dnuo, const vh, const vz);
	SampZ(const dphi, const dsigeta, const drho, const dbeta,
		  const dlam, const ddel, const vh, const vzo);
	SampXi(const dpsi, const dsigu, const vh);
	SampPsi(const dxi, const dsigu, const vh);
	SampSigu(const dxi, const dpsi, const vh);
	MultiMove(const dphi, const dsigeta, const drho,
			  const dgam, const dbeta, const dnu, const dxb, const dpsi, const dsigu, const vh, const vz,
			  const vzb);
	SingleMove(const dmu, const dphi, const dsigeta, const drho, const dbeta, const dnu, const dxi, const dsigu,
				const vho, const vz, const my);
	ML(const nI, const nR, const nsim);
	PFilter(const nI, const nR);

	fTsvar(const vx, const iBm);
	fGeweke(const vx, const iBm);
	
	fBetaMoment(const da, const db, const cm);
	fGammaMoment(const da, const db, const cm);
	fInvGammaMoment(const da, const db, const cm);
};
