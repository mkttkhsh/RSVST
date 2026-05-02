/*----------------------------------------------------------
                        RSV_skewFS.h
----------------------------------------------------------*/
											 
/*
**  HEADER FILE
**
*/

class RSV_skewFS
{
	decl m_mdata;
	decl m_cT;
	decl m_prediction_size, m_prediction_number;
	decl m_isini;
	decl m_k;
	decl m_isT, m_isSkew, m_isRSV;
	decl m_my, m_mx;
	decl m_vh, m_vlambda;
	decl m_dmu, m_dphi, m_drho, m_dgamma, m_dsigeta, m_dxi, m_dsigu, m_dnu;
	
	decl m_dm0_mu, m_dsig0_mu, m_dm0_xi, m_dsig0_xi;
	decl m_dn0_eta, m_dS0_eta, m_dn0_u, m_dS0_u;
	decl m_da0_phi, m_db0_phi, m_da0_rho, m_db0_rho;
	decl m_dn0_gamma, m_dS0_gamma, m_dn0_nu, m_dS0_nu;

	decl m_ML, m_BM;
	decl m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast;
	
	RSV_skewFS(const mdata, const cT, const prediction_size, const prediction_number, const isini);
	fMyWishart(const n, const mS);
	fcholeski(const mS);
	SetML(const BM, const dmu, const dphi, const drho, const dgamma, const dsigeta,
					 const dxi, const dsigu, const dnu, const vh);
	SetParam(const dmu, const dphi, const drho, const dgamma, const dsigeta,
					 const dxi, const dsigu, const dnu, const vh, const vlambda);
	SetIsT(const isT);
	SetIsSkew(const isSkew);
	SetIsRSV(const isRSV);
	GetParam(const addmu, const addphi, const addrho, const addgamma, const addsigeta,
					 const addxi, const addsigu, const addnu, const advh, const advlambda);
	SetPrior(const dm0_mu, const dsig0_mu, const dm0_xi, const dsig0_xi,
					 const dn0_eta, const dS0_eta, const dn0_u, const dS0_u,
					 const da0_phi, const db0_phi, const da0_rho, const db0_rho,
					 const dn0_gamma, const dS0_gamma, const dn0_nu, const dS0_nu);
	PrintAll();
	MCMC(const nsim, const cburn);
	CalcMu_Lambda(const dnu);
	CalcM1(const dnu);
	CalcM2(const dnu);
	CalcC_nu(const dnu);
	CalcMu_Ast(const dgamma, const dnu);
	CalcSig_Ast(const dgamma, const dnu);
	CalcMu_Ast_SN(const dgamma);
	CalcSig_Ast_SN(const dgamma);
	densFSST(const dz, const dgamma, const dnu);
	densFSSN(const vz, const dgamma);
	densNormal(const vX, const vMu, const vS2);
	logdensGa(const x, const a, const b);
	logdensIGa(const x, const a, const b);
	SampleMu();
	SamplePhi(const iProp = 0);
	SampleRho(const iProp = 0);
	Posterior_Rho_Max(const dg, const adFunc, const avScore, const amHessian);
	SampleGamma(const iProp = 0);
	Posterior_Gamma_Max(const dlog_gamma, const adFunc, const avScore, const amHessian);
	SampleSigeta();
	Prob_SampleSigeta(const dsigeta);
	SampleSigeta2(const iProp = 0);
	Posterior_SigEta_Max(const dlog_sigeta, const adFunc, const avScore, const amHessian);
	SampleXi();
	SampleSigu();
	SampleNu(const iProp = 0);
	Posterior_Nu_Max(const dlog_nu, const adFunc, const avScore, const amHessian);
	SampleNu_RW();
	SampleH();
	ProbH(const vh, const t);
	SampleLambda();
	ProbLambda(const dlambda, const t);
	ranFSN(const dgamma);
	ranFS(const dgamma, const dnu, const dlambda);
	CalculateNumRho();
	CalculateDenomRho();
	CalculateNumGamma();
	CalculateDenomGamma();
	CalculateNumSigeta();
	CalculateDenomSigeta();
	CalculateNumNu();
	CalculateDenomNu();
	CalculateML(const nsim, const nburn);
	SearchMin(const vpi, const du);
	WeightToRatio(const vw);
	CalculateLikelihood(const nI);
	CalculateLikelihood_PF(const nI);
	CalculatePrior();
	SimObs();
};