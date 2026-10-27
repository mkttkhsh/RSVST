#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<RSV_skewFS.h>
#include<ReportMCMC.ox>

RSV_skewFS::RSV_skewFS(const mdata, const cT, const prediction_size, const prediction_number, const isini){
	decl time_str = "summary_";
	time_str += date();
	if(prediction_size > 0){
		time_str += sprint("_",prediction_number);
	}
	time_str += ".txt"; 
	fopen(time_str, "l");
	println( "%80s", time() ~ "/" ~ date() );

	m_cT = cT;
	
	m_mdata = mdata;
	if(mdata == "sim"){
		m_my = zeros(1, cT);
		m_mx = zeros(1, cT);
	}else{
		println("real data");
		if(mdata == "DJIA-RV"){
			println("DJIA data (RV)");
			m_my = loadmat("../Data/DJIA_rcc_rv_rk_2000_2013.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/DJIA_rcc_rv_rk_2000_2013.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "DJIA-RK"){
			println("DJIA data (RK)");
			m_my = loadmat("../Data/DJIA_rcc_rv_rk_2000_2013.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/DJIA_rcc_rv_rk_2000_2013.csv")[0:(m_cT-1)][3]' * (100^2));
		}
		if(mdata == "Nikkei"){
			println("Nikkei data");
			m_my = loadmat("../Data/Nikkei_rcc_rk_19960311_20140930.csv")[305:305+(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/Nikkei_rcc_rk_19960311_20140930.csv")[305:305+(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "Nikkei-RV"){
			println("Nikkei data (RV)");
			m_my = loadmat("../Data/Nikkei_rcc_rv_19960311_20140930.csv")[305:305+(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/Nikkei_rcc_rv_19960311_20140930.csv")[305:305+(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "DJIA-RV-2009"){
			println("DJIA data 2009-2019 (RV)");
			m_my = loadmat("../Data/DJIA_rcc_rv_2009_2019.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/DJIA_rcc_rv_2009_2019.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "Nikkei-RV-2009"){
			println("Nikkei225 data 2009-2019 (RV)");
			m_my = loadmat("../Data/Nikkei_rcc_rv_2009_2019.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/Nikkei_rcc_rv_2009_2019.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "Nikkei-RK-2009"){
			println("Nikkei225 data 2009-2019 (RK)");
			m_my = loadmat("../Data/Nikkei_rcc_rk_2009_2019.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/Nikkei_rcc_rk_2009_2019.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		// 2026 revision: COVID-period (2012-2020) samples, added for the psi-free sensitivity analysis.
		if(mdata == "Nikkei-RV-2020"){
			println("Nikkei225 data 2012-2020 (RV)");
			m_my = loadmat("../Data/N225_rcc_rv5_rkth2_bv_medrv_2012_2020.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/N225_rcc_rv5_rkth2_bv_medrv_2012_2020.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "DJIA-RV-2020"){
			println("DJIA data 2012-2020 (RV, close-grid)");
			m_my = loadmat("../Data/DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv")[0:(m_cT-1)][1]' * 100;
			m_mx = log(loadmat("../Data/DJIA_rcc_rv5_rkth2_bv_medrv_2012_2020_tick.csv")[0:(m_cT-1)][2]' * (100^2));
		}
		if(mdata == "DJIA-RV-APPLE"){
			println("APPLE_2010-2019 (RV)");
			m_my = loadmat("../Data/AAPL_DR_2010-2019.csv")[0:(m_cT-1)][]';
			m_mx = log(loadmat("../Data/AAPL_RV_2010-2019.csv")[0:(m_cT-1)][]');
		}
		if(mdata == "DJ-RV"){
			println("DJ_2010-2019 (RV)");
			m_my = loadmat("../Data/DJ_DR2010-2019.csv")[0:(m_cT-1)][]';
			m_mx = log(loadmat("../Data/DJ_RCOV2010-2019.csv")[0:(m_cT-1)][]');
		}
		if(mdata == "TOPIX-RV"){
			println("TOPIX_2010-2020 (RV)");
			m_my = loadmat("../Data/TOPIX_DR2010-2020.csv")[0:(m_cT-1)][]';
			m_mx = log(loadmat("../Data/TOPIX_RCOV2010-2020.csv")[0:(m_cT-1)][]');
		}
		if(mdata == "SP500-RV"){
			println("SP500_2010-2019 (RV)");
			m_my = loadmat("../Data/SP500_DR2010-2019.csv")[0:(m_cT-1)][]';
			m_mx = log(loadmat("../Data/SP500_RCOV2010-2019.csv")[0:(m_cT-1)][]');
		}
	}

	if(prediction_size > 0){
		println(sprint("prediction_", prediction_number));
		m_my = m_my[][prediction_number:prediction_number + (m_cT - prediction_size - 1)];
		m_mx = m_mx[][prediction_number:prediction_number + (m_cT - prediction_size - 1)];
		m_cT = m_cT - prediction_size;
	}
	
	m_prediction_size = prediction_size;
	m_prediction_number = prediction_number;
	m_isini = isini;

	m_isT = 1;
	m_isSkew = 1;
	m_isRSV = 1;
	m_isPsi = 0;		// 0: psi fixed at 1 (original); 1: psi estimated (free)
	m_dpsi = 1;			// slope coefficient of log-RV equation: x_t = xi + psi*h_t + u_t
	m_dm0_psi = 1;		// psi ~ N(m0_psi, sig0_psi^2)
	m_dsig0_psi = 1;

	m_ML = 0;
	decl LINE = "--------------------------------------------";
	println(LINE, LINE);
}


RSV_skewFS::fMyWishart(const n, const mS){
	// X~W(n,S), X:(pxp) symmetrx matrix
	// f(X|n,S)
	// =  *|X|^{(n-p-1)/2}*exp-0.5*tr(S^{-1}X)
	// n: degrees of freedom
	// S: (p x p) symmetric parameter matrix
	// E(X_{ij})=n*s_{ij}, s_{ij}: (i,j) element of S.
	decl ci,cj,cp,mA,mL,mT,mX;
	cp=rows(mS);mT=zeros(cp,cp);
	if(n<cp){print("Error! d.f. too small.");}
	// generate A~W(n,I)
	for(ci=0;ci<cp;++ci){
	mT[ci][ci]=sqrt(ranchi(1,1,n-ci));
	for(cj=0;cj<ci;++cj){mT[ci][cj]=rann(1,1);}
	}
	mA=mT*mT'; // A=TT' ~W(n,I)
	mL=choleski(mS); // S=LL'(Choleski decomposition)
	mX=mL*mA*mL'; // X=LAL'~W(n,S)
return(mX);
}

RSV_skewFS::fcholeski(const mS)
{
	decl mA = choleski(mS);
	if(mA == 0){
		//println("Singular matrix!");
		decl vlambda, mP;
		eigensym(mS, &vlambda, &mP);
		mP ./= sqrt( sumsqrc(mP) );
		mA = mP * diag( sqrt(vlambda) );
	}
	return mA;
}

RSV_skewFS::SetML(const BM, const dmu, const dphi, const drho, const dgamma, const dsigeta,
					 const dxi, const dsigu, const dnu, const vh){
	m_ML = 1;
	m_BM = BM;
	m_dmu_ast = dmu;
	m_dphi_ast = dphi;
	m_drho_ast = drho;
	m_dgamma_ast = dgamma;
	m_dsigeta_ast = dsigeta;
	m_dxi_ast = dxi;
	m_dsigu_ast = dsigu;
	m_dnu_ast = dnu;
	m_vh_ast = vh;
}

//Setter of Parameters
RSV_skewFS::SetParam(const dmu, const dphi, const drho, const dgamma, const dsigeta,
					 const dxi, const dsigu, const dnu, const vh, const vlambda, const dpsi){
	m_dpsi = m_isPsi ? dpsi : 1;
	m_dmu = dmu;
	m_dphi = dphi;
	m_drho = drho;
	if(m_isSkew == 1){
		m_dgamma = dgamma;
	}else{
		m_dgamma = 1;
	}
	m_dsigeta = dsigeta;
	m_dxi = dxi;
	m_dsigu = dsigu;
	m_dnu = dnu;
	if(m_mdata == "sim"){
		m_vh = vh;
	}else{
		m_vh = m_mx;
	}
	if(m_prediction_size > 0 && m_isini == 0){
		m_vh = vh;
	}
	if(m_ML == 1){
		m_vh = vh;
	}
	m_vlambda = vlambda;
}

RSV_skewFS::SetIsT(const isT){
	m_isT = isT;
}

RSV_skewFS::SetIsSkew(const isSkew){
	m_isSkew = isSkew;
}

RSV_skewFS::SetIsRSV(const isRSV){
	m_isRSV = isRSV;
}

RSV_skewFS::SetIsPsi(const isPsi){
	m_isPsi = isPsi;	// 1: estimate slope coefficient psi (free); 0: fix psi at 1
	if(!m_isPsi) m_dpsi = 1;
}

//Getter of Parameters
RSV_skewFS::GetParam(const addmu, const addphi, const addrho, const addgamma, const addsigeta,
					 const addxi, const addsigu, const addnu, const advh, const advlambda, const adpsi){
	addmu[0] = m_dmu;
	addphi[0] = m_dphi;
	addrho[0] = m_drho;
	addgamma[0] = m_dgamma;
	addsigeta[0] = m_dsigeta;
	addxi[0] = m_dxi;
	addsigu[0] = m_dsigu;
	addnu[0] = m_dnu;
	advh[0] = m_vh;
	if(adpsi != 0) adpsi[0] = m_dpsi;	// optional: return psi when a pointer is supplied
}

RSV_skewFS::SetPrior(const dm0_mu, const dsig0_mu, const dm0_xi, const dsig0_xi,
					 const dn0_eta, const dS0_eta, const dn0_u, const dS0_u,
					 const da0_phi, const db0_phi, const da0_rho, const db0_rho,
					 const dn0_gamma, const dS0_gamma, const dn0_nu, const dS0_nu){
	println("PRIOR SETTING");
	println("-------------");
	m_dm0_mu = dm0_mu;
	println("dm0_mu = ", dm0_mu);
	m_dsig0_mu = dsig0_mu;
	println("dsig0_mu = ", dsig0_mu);
	m_dm0_xi = dm0_xi;
	println("dm0_xi = ", dm0_xi);
	m_dsig0_xi = dsig0_xi;
	println("dsig0_xi = ", dsig0_xi);

	m_dn0_eta = dn0_eta;
	println("dn0_eta = ", dn0_eta);
	m_dS0_eta = dS0_eta;
	println("dS0_eta = ", dS0_eta);
	m_dn0_u = dn0_u;
	println("dn0_u = ", dn0_u);
	m_dS0_u = dS0_u;
	println("dS0_u = ", dS0_u);
	
	m_da0_phi = da0_phi;
	println("da0_phi = ", da0_phi);
	m_db0_phi = db0_phi;
	println("db0_phi = ", db0_phi);
	m_da0_rho = da0_rho;
	println("da0_rho = ", da0_rho);
	m_db0_rho = db0_rho;
	println("db0_rho = ", db0_rho);
	
	m_dn0_gamma = dn0_gamma;
	println("dn0_gamma = ", dn0_gamma);
	m_dS0_gamma = dS0_gamma;
	println("dS0_gamma = ", dS0_gamma);
	m_dn0_nu = dn0_nu;
	println("dn0_nu = ", dn0_nu);
	m_dS0_nu = dS0_nu;
	println("dS0_nu = ", dS0_nu);
	
	decl LINE = "--------------------------------------------";
	println(LINE, LINE);
}

RSV_skewFS::PrintAll(){
	println("mu = ", m_dmu);
	println("phi = ", m_dphi);
	println("rho = ", m_drho);
	println("gamma = ", m_dgamma);
	println("sigma_eta = ", m_dsigeta);
	println("xi = ", m_dxi);
	println("sigu = ", m_dsigu);
	println("nu = ", m_dnu);
	println("mean of h = ", meanr(m_vh));
	println("mean of lambda = ", meanr(m_vlambda));

}

RSV_skewFS::MCMC(const nsim, const cburn){

	println("INITIAL VALUES SETTING");
	println("----------------------");
	println("mu_ini = ", m_dmu);
	println("phi_ini = ", m_dphi);
	println("rho_ini = ", m_drho);
	println("gamma_ini = ", m_dgamma);
	println("sigeta_ini = ", m_dsigeta);
	println("xi_ini = ", m_dxi);
	println("sigu_ini = ", m_dsigu);
	println("nu_ini = ", m_dnu);

	decl LINE = "--------------------------------------------";
	println(LINE, LINE);
	
	decl dmu_s = zeros(nsim, 1);
	decl dphi_s = zeros(nsim, 1);
	decl drho_s = zeros(nsim, 1);
	decl dgamma_s = zeros(nsim, 1);
	decl dsigeta_s = zeros(nsim, 1);
	decl dxi_s = zeros(nsim, 1);
	decl dsigu_s = zeros(nsim, 1);
	decl dpsi_s = zeros(nsim, 1);
	decl dnu_s = zeros(nsim, 1);
	decl vh_s = zeros(nsim, 1);
	decl vlambda_s = zeros(nsim, 1);
	decl mh_full_s = zeros(nsim, m_cT);

	decl vlambda_pred, veps_pred, veta_pred, vu_pred, vh_pred, mx_pred, my_pred;
	decl vh_pred_mean = 0;
	decl mx_pred_mean = 0;
	decl my_pred_mean = 0;

	decl vh_pred_s = zeros(nsim, 1);
	decl mx_pred_s = zeros(nsim, 1);
	decl my_pred_s = zeros(nsim, 1);
	
	decl progress_bar = 0;

	if(m_isSkew != 1){
		m_dgamma = 1;
	}
	
	for(m_k=-cburn;m_k<nsim+1;m_k++){
		SampleMu();
		SamplePhi();
		SampleRho();
		if(m_isSkew == 1){
			SampleGamma();
		}
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SamplePsi();
			SampleSigu();
		}
		SampleNu();
		SampleH();

		if(m_k == -cburn){
			println("Progress");
			println("--|--|--|--|--|--|--|--|--|--|");
		}
		if(100 * ((m_k + cburn) / (nsim + cburn)) >= (progress_bar + 1) * 5){
			print("=");
			progress_bar += 1;
			if(imod(progress_bar, 2) == 0 && progress_bar != 0){
				print("|");
			}
		}
		if(m_k > 0){
			dmu_s[m_k - 1][] = m_dmu;
			dphi_s[m_k - 1][] = m_dphi;
			drho_s[m_k - 1][] = m_drho;
			dgamma_s[m_k - 1][] = m_dgamma;
			dsigeta_s[m_k - 1][] = m_dsigeta;
			dxi_s[m_k - 1][] = m_dxi;
			dsigu_s[m_k - 1][] = m_dsigu;
			dpsi_s[m_k - 1][] = m_dpsi;
			dnu_s[m_k - 1][] = m_dnu;
			vh_s[m_k - 1][] = meanr(m_vh)';
			mh_full_s[m_k - 1][] = m_vh;
			if(m_prediction_size > 0){
				veta_pred = m_drho * m_dsigeta * (m_my[][m_cT-1]) * exp(-m_vh[][(m_cT-1)]/2) + sqrt((1 - m_drho^2) * m_dsigeta^2) * rann(1, 1);
				vu_pred = m_dsigu * rann(1, 1);
				vh_pred = m_dmu + m_dphi * (m_vh[][(m_cT-1)] - m_dmu) + veta_pred;
				mx_pred = m_dxi + m_dpsi * vh_pred + vu_pred;
				vlambda_pred = invert(rangamma(1, 1, m_dnu/2, m_dnu/2));
				veps_pred = ranFS(m_dgamma, m_dnu, vlambda_pred);
				my_pred = veps_pred * exp(vh_pred/2);

				vh_pred_mean += (vh_pred - vh_pred_mean)/m_k;
				mx_pred_mean += (mx_pred - mx_pred_mean)/m_k;
				my_pred_mean += (my_pred - my_pred_mean)/m_k;

				vh_pred_s[m_k - 1][] = vh_pred;
				mx_pred_s[m_k - 1][] = mx_pred;
				my_pred_s[m_k - 1][] = my_pred;
			}
		}
	}
	if(m_prediction_size == 0){
		savemat(sprint("h_quantile.csv"), quantiler(mh_full_s', <0.025,0.975>));
		savemat(sprint("h_mean.csv"), meanr(mh_full_s'));
		println("\n");
		decl out;
		out = new ReportMCMC(dmu_s);
		out.SetOutfileName("dmu");
		println("--------------------------\n\n-------------Mu--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(dphi_s);
		out.SetOutfileName("dphi");
		println("--------------------------\n\n-------------Phi--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(drho_s);
		out.SetOutfileName("drho");
		println("--------------------------\n\n-------------Rho--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(dgamma_s);
		out.SetOutfileName("dgamma");
		println("--------------------------\n\n-------------Gamma--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(dsigeta_s);
		out.SetOutfileName("dsigeta");
		println("--------------------------\n\n-------------SigEta--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(dxi_s);
		out.SetOutfileName("dxi");
		println("--------------------------\n\n-------------Xi--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(dsigu_s);
		out.SetOutfileName("dsigu");
		println("--------------------------\n\n-------------SigU--------------");
		out.Report();
		delete out;
		if(m_isPsi == 1){
			out = new ReportMCMC(dpsi_s);
			out.SetOutfileName("dpsi");
			println("--------------------------\n\n-------------Psi--------------");
			out.Report();
			delete out;
		}
		out = new ReportMCMC(dnu_s);
		out.SetOutfileName("dnu");
		println("--------------------------\n\n-------------Nu--------------");
		out.Report();
		delete out;
		out = new ReportMCMC(vh_s);
		out.SetOutfileName("vh");
		println("--------------------------\n\n-------------H--------------");
		out.Report();
		delete out;
	}
	if(m_prediction_size > 0){
		savemat(sprint("h_pred_s_", m_prediction_number,".csv"), vh_pred_s);
		savemat(sprint("x_pred_s_", m_prediction_number,".csv"), mx_pred_s);
		savemat(sprint("y_pred_s_", m_prediction_number,".csv"), my_pred_s);
		savemat(sprint("par_est_s_", m_prediction_number,".csv"), dmu_s ~ dphi_s ~ drho_s ~ dgamma_s ~ dsigeta_s ~ dsigeta_s ~ dxi_s ~ dsigu_s ~ dnu_s ~ dpsi_s);
		savemat(sprint("h_est_s_", m_prediction_number,".csv"), meanc(mh_full_s)' ~ quantilec(mh_full_s, <0.025,0.5,0.975>)');
		m_dmu = meanc(dmu_s);
		m_dphi = meanc(dphi_s);
		m_drho = meanc(drho_s);
		m_dgamma = meanc(dgamma_s);
		m_dsigeta = meanc(dsigeta_s);
		m_dxi = meanc(dxi_s);
		m_dsigu = meanc(dsigu_s);
		m_dnu = meanc(dnu_s);
		m_dpsi = m_isPsi ? meanc(dpsi_s) : 1;
		m_vh = meanc(mh_full_s)[][1:(m_cT-1)] ~ vh_pred_mean;
	}
	if(m_ML == 1){
		m_dmu = meanc(dmu_s);
		m_dphi = meanc(dphi_s);
		m_drho = meanc(drho_s);
		m_dgamma = meanc(dgamma_s);
		m_dsigeta = meanc(dsigeta_s);
		m_dxi = meanc(dxi_s);
		m_dsigu = meanc(dsigu_s);
		m_dnu = meanc(dnu_s);
		m_dpsi = m_isPsi ? meanc(dpsi_s) : 1;
		m_vh = meanc(mh_full_s);
	}
}

RSV_skewFS::CalcMu_Lambda(const dnu){
	return dnu / (dnu - 2);
}

RSV_skewFS::CalcM1(const dnu){
	return (2 * CalcC_nu(dnu) * dnu) / (dnu - 1);
}

RSV_skewFS::CalcM2(const dnu){
	return dnu / (dnu - 2);
}

RSV_skewFS::CalcC_nu(const dnu){
	return gammafact((dnu + 1)/2) / (gammafact(dnu/2) * gammafact(1/2) * sqrt(dnu));
}

RSV_skewFS::CalcMu_Ast(const dgamma, const dnu){
	return CalcM1(dnu) * (dgamma - dgamma^(-1));
}

RSV_skewFS::CalcSig_Ast(const dgamma, const dnu){
	return sqrt( CalcM2(dnu) * ( (dgamma^3 + dgamma^(-3))/(dgamma + dgamma^(-1)) ) - CalcMu_Ast(dgamma, dnu)^2 );
}

RSV_skewFS::CalcMu_Ast_SN(const dgamma){
	return sqrt(2/M_PI) * (dgamma - dgamma^(-1));
}

RSV_skewFS::CalcSig_Ast_SN(const dgamma){
	return sqrt( 1 * ( (dgamma^3 + dgamma^(-3))/(dgamma + dgamma^(-1)) ) - CalcMu_Ast_SN(dgamma)^2 );
}

RSV_skewFS::densFSST(const vz, const dgamma, const dnu){
	decl dsig_ast =	CalcSig_Ast(dgamma, dnu);
	decl dmu_ast = CalcMu_Ast(dgamma, dnu);
	return (2 * dgamma * dsig_ast * CalcC_nu(dnu))/(1 + dgamma^2)
		  * (1 + ( (dsig_ast .* vz + dmu_ast).^2 / dnu )
		   	     .* ( dgamma^(-2) .* ( vz .>= - dmu_ast / dsig_ast )
					 + dgamma^2 .* (vz .< - dmu_ast / dsig_ast )
				    ) 
			).^(-(dnu+1)/2)
			; 

}

RSV_skewFS::densFSSN(const vz, const dgamma){
	decl dsig_ast =	CalcSig_Ast_SN(dgamma);
	decl dmu_ast = CalcMu_Ast_SN(dgamma);
	return (2 * dgamma * dsig_ast)/(1 + dgamma^2)
		  *	(1/sqrt(2*M_PI))
		  * exp( -( (dsig_ast .* vz + dmu_ast).^2 / 2 )
		   	     .* ( dgamma^(-2) .* ( vz .>= - dmu_ast / dsig_ast )
					 + dgamma^2 .* (vz .< - dmu_ast / dsig_ast )
				    ) 
			)
			;
}

RSV_skewFS::densNormal(const vX, const vMu, const vS2){
    decl s = sqrt( vS2 );
    decl p = densn( ( vX - vMu ) ./ s ) ./ s;
    
    return p;
}

RSV_skewFS::logdensGa(const x, const a, const b){
    return a .* log( b ) - loggamma( a ) + ( a - 1 ) .* log(x) - b .* x;
}

RSV_skewFS::logdensIGa(const x, const a, const b){
   return a*log(b) - loggamma(a) - (a + 1.0)*log(x) - b./x;
}

RSV_skewFS::SampleMu(){
	decl dsig_mu = (m_dsig0_mu^(-2) + m_dsigeta^(-2) * (1 - m_dphi^2 + (m_cT - 1) * (1 - m_dphi)^2 * (1 - m_drho^2)^(-1) ) )^(-1/2);
	decl dm_mu = dsig_mu^2 * (m_dsig0_mu^(-2) * m_dm0_mu
								+ m_dsigeta^(-2) *
									(  (1 - m_dphi^2) * m_vh[][0]
									+ ((1 - m_dphi) / (1 - m_drho^2))
										* sumr(m_vh[][1:(m_cT-1)] - m_dphi * m_vh[][0:(m_cT-2)]
											- m_drho * m_dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)]
						      	          	)
									)
					  		  );
	m_dmu = dm_mu + dsig_mu * rann(1, 1);
}

RSV_skewFS::SamplePhi(const iProp){
	decl dsig_phi = sqrt( (m_dsigeta^2 * (1 - m_drho^2)) / (m_drho^2 * (m_vh[][0] - m_dmu)^2 + sumr((m_vh[][1:(m_cT-1)] - m_dmu).^2) ) );
	decl dm_phi = (dsig_phi^2 / (m_dsigeta^2 * (1 - m_drho^2)) )
					* sumr(
						(m_vh[][1:(m_cT-1)] - m_dmu - m_drho .* m_dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)]) 
						.* (m_vh[][0:(m_cT-2)] - m_dmu)
						);
	decl dphio = m_dphi;
	decl dphin = dm_phi + dsig_phi * rann(1, 1);
	decl dfrac = ((1 - dphin)^(m_da0_phi - (1/2)) * (1 + dphin)^(m_db0_phi - (1/2)))
				 / (1 - dphio)^(m_da0_phi - (1/2)) * (1 + dphio)^(m_db0_phi - (1/2));
	if((ranu(1, 1) < dfrac || iProp == 1) && fabs(dphin) < 1 ){
		m_dphi = dphin;
	}
}

RSV_skewFS::SampleRho(const iProp){
	decl dghat, drhohat;
	decl Func, drv1, drv2;
	dghat = (log(1+m_drho) - log(1-m_drho));
	MaxBFGS(Posterior_Rho_Max, &dghat, &Func, 0, 1);
	Num2Derivative(Posterior_Rho_Max, dghat, &drv2);
	drv2 *= 20000;// To prevent the number getting too small 
	decl dsig_rho = sqrt(-invert(drv2));

	decl dm_rho = dghat;
	decl dgo = (log(1+m_drho) - log(1-m_drho));
	decl dgn = dm_rho + dsig_rho * rann(1, 1);
	Posterior_Rho_Max(dgn, &Func, 0, 0);
	decl dln = 20000 * Func;
	Posterior_Rho_Max(dgo, &Func, 0, 0);
	decl dlo = 20000 * Func;	
	decl dfrac = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_rho^2) * (dgo - dghat)^2)
				* exp((1/2) * invert(dsig_rho^2) * (dgn - dghat)^2);
	if(ranu(1, 1) < dfrac || iProp == 1){
		m_drho = (exp(dgn) - 1)/(exp(dgn) + 1);
	}	
}

RSV_skewFS::Posterior_Rho_Max(const dg, const adFunc, const avScore, const amHessian){
	decl drho = (exp(dg) - 1)/(exp(dg) + 1);
	adFunc[0] =	log( (1 + drho)^(m_da0_rho - 1) * (1 - drho)^(m_db0_rho - 1) )
				 + log(1 - drho^2) * (-(m_cT-1)/2)
				 -(1/(2 * (1 - drho^2) * m_dsigeta^2))
						* sumr( ( (m_vh[][1:(m_cT-1)] - m_dmu - m_dphi * (m_vh[][0:(m_cT-2)] - m_dmu) )
								  - drho .* m_dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)]
								).^2
							)
				;
	adFunc[0] += log(2 * exp(dg) * (exp(dg) + 1)^(-2));//Jacobian
	adFunc[0] *= (1/20000);// restore the scale
	return 1;
}

RSV_skewFS::SampleGamma(const iProp){
	if(m_isSkew != 1){
		return 1;
	}
	decl dloghat, dgammahat;
	decl Func, drv1, drv2;
	dloghat = log(m_dgamma);
	MaxBFGS(Posterior_Gamma_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Gamma_Max, dloghat, &drv2);

	decl dsig_gamma = sqrt(invert(-drv2));
	decl dm_gamma = dloghat;
	decl dlogo = log(m_dgamma);
	decl dlogn = dm_gamma + dsig_gamma * rann(1, 1);
	Posterior_Gamma_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Gamma_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl dfrac = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_gamma^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_gamma^2) * (dlogn - dloghat)^2);
	if(ranu(1, 1) < dfrac || iProp == 1){
		m_dgamma = exp(dlogn);
	}

}

RSV_skewFS::Posterior_Gamma_Max(const dlog_gamma, const adFunc, const avScore, const amHessian){
	decl dgamma = exp(dlog_gamma);
	decl dmu_ast, dsig_ast;
	if(m_isT == 1){
		dmu_ast = CalcMu_Ast(dgamma, m_dnu);
		dsig_ast =	CalcSig_Ast(dgamma, m_dnu);
	}else{
		dmu_ast = CalcMu_Ast_SN(dgamma);
		dsig_ast =	CalcSig_Ast_SN(dgamma);
	}
	
	if(m_isT == 1){
		adFunc[0] =	log(dgamma^(m_dn0_gamma - 1)) - m_dS0_gamma * dgamma
					 + sumr(log(densFSST(m_my .* exp(-m_vh./2), dgamma, m_dnu)))
				 ;
	}else{
		adFunc[0] =	log(dgamma^(m_dn0_gamma - 1)) - m_dS0_gamma * dgamma
					 + sumr(log(densFSSN(m_my .* exp(-m_vh./2), dgamma)))
				 ;
	}
	adFunc[0] += log(exp(dlog_gamma));
	return 1;
}

RSV_skewFS::SampleSigeta(){
	decl dn1_eta = m_dn0_eta + m_cT;
	decl dS1_eta = m_dS0_eta + (1 - m_dphi^2) * (m_vh[][0] - m_dmu)^2
							 + (1 / (1 - m_drho^2)) * sumr( ( m_vh[][1:(m_cT-1)] - m_dmu - m_dphi * (m_vh[][0:(m_cT-2)] - m_dmu) ).^2 );
	decl dsigetao = m_dsigeta;
	decl dsigetan = sqrt(invert(rangamma(1, 1, dn1_eta/2, dS1_eta/2)));
	decl dfrac = exp(Prob_SampleSigeta(dsigetan) - Prob_SampleSigeta(dsigetao));
	if(ranu(1, 1) < dfrac){
		m_dsigeta = dsigetan;
	}
}

//Log value
RSV_skewFS::Prob_SampleSigeta(const dsigeta){
	return dsigeta^(-1)
				* sumr(((m_drho * m_my[][0:(m_cT-2)] .* exp(-m_vh[][0:(m_cT-2)]./2)) / (1 - m_drho^2))
							.* (m_vh[][1:(m_cT-1)] - m_dmu - m_dphi .* (m_vh[][0:(m_cT-2)] - m_dmu))
					)
			;
}

RSV_skewFS::SampleSigeta2(const iProp){
	decl dloghat, dsighat;
	decl Func, drv1, drv2;
	dloghat = log(m_dsigeta);
	MaxBFGS(Posterior_SigEta_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_SigEta_Max, dloghat, &drv2);

	decl dsig_sigeta = sqrt(invert(-drv2));
	decl dm_sigeta = dloghat;
	decl dlogo = log(m_dsigeta);
	decl dlogn = dm_sigeta + dsig_sigeta * rann(1, 1);
	Posterior_SigEta_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_SigEta_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl dfrac = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_sigeta^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_sigeta^2) * (dlogn - dloghat)^2);
	if(ranu(1, 1) < dfrac || iProp == 1){
		m_dsigeta = exp(dlogn);
	}
}

RSV_skewFS::Posterior_SigEta_Max(const dlog_sigeta, const adFunc, const avScore, const amHessian){
	decl dsigeta = exp(dlog_sigeta);
	adFunc[0] =	(-(m_dn0_eta + m_cT)/2 - 1) * log(dsigeta^2) + (-m_dS0_eta / (2 * dsigeta^2))
			    + (-1 / (2 * (1 - m_drho^2) *dsigeta^2) )
				* sumr( ( m_vh[][1:(m_cT-1)] - m_dmu - m_dphi .* (m_vh[][0:(m_cT-2)] - m_dmu) - m_drho .* dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)] ).^2 )
				 ;
	adFunc[0] += log(exp(dlog_sigeta));
	return 1;
}

RSV_skewFS::SampleXi(){
	if(m_isRSV != 1){
		return 1;
	}
	decl dsig_xi = sqrt(invert(m_dsig0_xi^(-2) + m_cT * m_dsigu^(-2)));
	decl dmu_xi = dsig_xi^2 * (m_dsig0_xi^(-2) * m_dm0_xi + m_dsigu^(-2) * sumr(m_mx - m_dpsi * m_vh));
	m_dxi = dmu_xi + dsig_xi * rann(1, 1);
}

// SamplePsi samples the slope coefficient psi of the log-RV equation
//   x_t = xi + psi*h_t + u_t,  u_t ~ N(0, sigu^2),  psi ~ N(m0_psi, sig0_psi^2)
RSV_skewFS::SamplePsi(){
	if(m_isRSV != 1 || m_isPsi != 1){
		return 1;
	}
	decl dsig_psi = sqrt(invert(m_dsig0_psi^(-2) + m_dsigu^(-2) * sumr(m_vh.^2)));
	decl dmu_psi = dsig_psi^2 * (m_dsig0_psi^(-2) * m_dm0_psi + m_dsigu^(-2) * sumr((m_mx - m_dxi) .* m_vh));
	m_dpsi = dmu_psi + dsig_psi * rann(1, 1);
}

RSV_skewFS::SampleSigu(){
	if(m_isRSV != 1){
		return 1;
	}
	decl dn1_u = m_dn0_u + m_cT;
	decl dS1_u = m_dS0_u + sumr((m_mx - m_dxi - m_dpsi * m_vh).^2);
	m_dsigu = sqrt(invert(rangamma(1, 1, dn1_u/2, dS1_u/2)));
}

RSV_skewFS::SampleLambda(){
	decl da, db, dg;
	decl t;
	decl dlambdao, dlambdan;
	decl dsig_ast =	CalcSig_Ast(m_dgamma, m_dnu);
	decl dmu_ast = CalcMu_Ast(m_dgamma, m_dnu);
	decl dfrac;
	for(t = 0; t < m_cT; t++){
		da = m_dnu + 1;
		db = m_dnu + (CalcMu_Lambda(m_dnu) * m_drho^2 * exp(-m_vh[][t]/2) * m_my[][t] / (1 - m_drho^2)) * (t < m_cT - 1)
				   + CalcMu_Lambda(m_dnu) * (dsig_ast * m_my[][t] * exp(-m_vh[][t]/2) + dmu_ast)^2
				   					 * (m_dgamma^(-2) * (m_my[][t] * exp(-m_vh[][t]/2) >= -dmu_ast/dsig_ast)
									 	 + m_dgamma^2 * (m_my[][t] * exp(-m_vh[][t]/2) < -dmu_ast/dsig_ast)
									 	);
		dlambdao = m_vlambda[][t];
		dlambdan = invert(rangamma(1, 1, da/2, db/2));
		if(t < m_cT - 1){
			dfrac = ProbLambda(dlambdan, t) / ProbLambda(dlambdao, t); 
		}else{
			dfrac = 1;
		}
		if(ranu(1, 1) < dfrac){
			 m_vlambda[][t] = dlambdan;
		}
	}

}

RSV_skewFS::ProbLambda(const dlambda, const t){
	return exp( sqrt(CalcMu_Lambda(m_dnu)/dlambda)
						* (m_drho * m_my[][t]) * exp(-m_vh[][t]/2) / ((1 - m_drho^2) * m_dsigeta)
						* (m_vh[][t+1] - m_dmu - m_dphi * (m_vh[][t] - m_dmu))
			);
}

RSV_skewFS::SampleNu(const iProp){
	if(m_isT != 1){
		return 1;
	}
	decl dloghat, dnuhat;
	decl Func, drv1, drv2;
	dloghat = log(m_dnu - 2);
	MaxBFGS(Posterior_Nu_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Nu_Max, dloghat, &drv2);
	decl dsig_nu = invert(sqrt(-drv2));
	decl dm_nu = dloghat;
	decl dlogo = log(m_dnu - 2);
	decl dlogn = dm_nu + dsig_nu * rann(1, 1);
	Posterior_Nu_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Nu_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;
	decl dfrac = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_nu^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_nu^2) * (dlogn - dloghat)^2);
	if(ranu(1, 1) < dfrac || iProp == 1){
		m_dnu = exp(dlogn) + 2;
	}
}

RSV_skewFS::Posterior_Nu_Max(const dlog_nu, const adFunc, const avScore, const amHessian){
	decl dnu = exp(dlog_nu) + 2;
	decl dsig_ast =	CalcSig_Ast(m_dgamma, dnu);
	decl dmu_ast = CalcMu_Ast(m_dgamma, dnu);
	adFunc[0] =	log(dnu^(m_dn0_nu - 1)) - m_dS0_nu * dnu
				 + sumr(log(densFSST(m_my .* exp(-m_vh./2), m_dgamma, dnu)))
				 ;
	adFunc[0] += log(exp(dlog_nu));
	return 1;
}

RSV_skewFS::SampleH(){
	decl dm_h, dsig_h;
	decl dho, dhn;
	decl vhn, vho;
	decl dfrac;
	decl num_accept = 0;
	decl t;
	for(t=0;t<m_cT;t++){
		dsig_h = sqrt(invert(
				 m_dsigeta^(-2)
				 * ( (t == 0) ? (1 - m_dphi^2 + (m_dphi^2/(1 - m_drho^2))) : ( (t < m_cT - 1) ? ((1 + m_dphi^2) / (1 - m_drho^2)) : (1 / (1 - m_drho^2)) ) )
				 + ((m_isRSV == 1) ? m_dpsi^2 * m_dsigu^(-2) : 0)
				 ));
		dm_h = dsig_h^2 * (
						  (-1/2)
						  + ( (t == m_cT - 1) ? 0 : (( m_dphi * (m_vh[][t+1] - (1 - m_dphi) * m_dmu) ) / ( m_dsigeta^2 * (1 - m_drho^2) )) )
						  + ( (t == 0) ? ((1 - m_dphi^2) * m_dmu / m_dsigeta^2) : ( ( (1 - m_dphi) * m_dmu + m_dphi * m_vh[][t-1] ) / ( m_dsigeta^2 * (1 - m_drho^2) )) )
						  + ( (t == 0) ? 0 : ( m_drho * m_my[][t-1] * exp(-m_vh[][t-1]/2) / ( (1 - m_drho^2) * m_dsigeta ) ) )
						  + ( (m_isRSV == 1) ? (m_dpsi * (m_mx[][t] - m_dxi)/m_dsigu^2) : 0 )
						  );
		dho = m_vh[][t];
		dhn = dm_h + dsig_h * rann(1, 1);
		vho = m_vh;
		vhn = vho;
		vhn[][t] = dhn;
		dfrac = exp(ProbH(vhn, t) - ProbH(vho, t));
		if(ranu(1, 1) < dfrac){
			m_vh[][t] = dhn;
			num_accept += 1;
		}
	}												  
}

RSV_skewFS::ProbH(const vh, const t){
	decl dmu_ast, dsig_ast;
	if(m_isT == 1){
		dmu_ast = CalcMu_Ast(m_dgamma, m_dnu);
		dsig_ast =	CalcSig_Ast(m_dgamma, m_dnu);
	}else{
		dmu_ast = CalcMu_Ast_SN(m_dgamma);
		dsig_ast =	CalcSig_Ast_SN(m_dgamma);
	}
	return	(
			(m_isT == 1) ? 
			(-(m_dnu + 1)/2) * log(1 + ((dsig_ast * m_my[][t] * exp(-vh[][t]/2) + dmu_ast)^2 / m_dnu)
			   * ( ( m_dgamma^(-2) * (m_my[][t] * exp(-vh[][t]/2) >= - dmu_ast/dsig_ast) )
				  + ( m_dgamma^2 * (m_my[][t] * exp(-vh[][t]/2) < - dmu_ast/dsig_ast) )
				 )
		    ) :
		    ((-(dsig_ast * m_my[][t] * exp(-vh[][t]/2) + dmu_ast)^2 / 2)
			   * ( ( m_dgamma^(-2) * (m_my[][t] * exp(-vh[][t]/2) >= - dmu_ast/dsig_ast) )
				  + ( m_dgamma^2 * (m_my[][t] * exp(-vh[][t]/2) < - dmu_ast/dsig_ast) )
				 )
			)
			)
		 + ((t == m_cT - 1) ?  0 :
		 	 - m_drho^2 * m_my[][t]^2 * (exp(-vh[][t]) / (2 * (1 - m_drho^2)))
			   + ((m_drho * m_my[][t] * exp(-vh[][t]/2)) / ((1 - m_drho^2) * m_dsigeta))
			     * (vh[][t+1] - m_dmu - m_dphi * (vh[][t] - m_dmu))
			 
		   );
}

RSV_skewFS::ranFSN(const dgamma){
	decl dsig_ast =	CalcSig_Ast_SN(dgamma);
	decl dmu_ast = CalcMu_Ast_SN(dgamma);
	decl deps0;
	//generate deps0 ~ FSSN
	if(ranu(1, 1) <  (1/dgamma)/(dgamma + (1/dgamma)) ){
		deps0 = (1 / dgamma) * (-fabs(rann(1, 1)));
	}else{
		deps0 = (1 * dgamma) * (fabs(rann(1, 1)));
	}
	return (deps0 - dmu_ast) / dsig_ast;	
}

RSV_skewFS::ranFS(const dgamma, const dnu, const dlambda){
	decl dsig_ast =	CalcSig_Ast(dgamma, dnu);
	decl dmu_ast = CalcMu_Ast(dgamma, dnu);
	decl deps0;
	//generate deps0 ~ FSST, dlambda ~ IG(dnu/2, dnu/2)
	if(ranu(1, 1) <  (1/dgamma)/(dgamma + (1/dgamma)) ){
		deps0 = (sqrt(dlambda) / dgamma) * (-fabs(rann(1, 1)));
	}else{
		deps0 = (sqrt(dlambda) * dgamma) * (fabs(rann(1, 1)));
	}
	return (deps0 - dmu_ast) / dsig_ast;
}

RSV_skewFS::CalculateNumRho(){
	decl dghat, drhohat;
	decl Func, drv1, drv2;
	dghat = (log(1+m_drho) - log(1-m_drho));
	MaxBFGS(Posterior_Rho_Max, &dghat, &Func, 0, 1);
	Num2Derivative(Posterior_Rho_Max, dghat, &drv2);
	decl dsig_rho = sqrt(-invert(drv2));
	decl dm_rho = dghat;

	decl dgn, dgo;
	decl dln, dlo;
	decl dnum;

	dgn = log(1+m_drho_ast) - log(1-m_drho_ast);
	Posterior_Rho_Max(dgn, &Func, 0, 0);
	dln = Func;
	dgo = log(1+m_drho) - log(1-m_drho);
	Posterior_Rho_Max(dgo, &Func, 0, 0);
	dlo = Func;

	dnum = exp(dln - dlo)
			* exp(-(1/2) * invert(dsig_rho^2) * (dgo - dghat)^2)
			* exp((1/2) * invert(dsig_rho^2) * (dgn - dghat)^2)
			* densNormal(dgn, dm_rho, dsig_rho^2)
			* 2/((1 + m_drho_ast) * (1 - m_drho_ast));
	return dnum;
}

RSV_skewFS::CalculateDenomRho(){
	decl dghat, drhohat;
	decl Func, drv1, drv2;
	dghat = (log(1+m_drho_ast) - log(1-m_drho_ast));
	MaxBFGS(Posterior_Rho_Max, &dghat, &Func, 0, 1);
	Num2Derivative(Posterior_Rho_Max, dghat, &drv2);
	drv2 *= 20000;
	decl dsig_rho = sqrt(-invert(drv2));
	decl dm_rho = dghat;

	decl dgn, dgo;
	decl dln, dlo;
	decl dnum, ddenom;
	
	dgn = log(1+m_drho) - log(1-m_drho);
	Posterior_Rho_Max(dgn, &Func, 0, 0);
	dln = 20000 * Func;
	dgo = log(1+m_drho_ast) - log(1-m_drho_ast);
	Posterior_Rho_Max(dgo, &Func, 0, 0);
	dlo = 20000 * Func;
	ddenom = exp(dln - dlo)
			* exp(-(1/2) * invert(dsig_rho^2) * (dgo - dghat)^2)
			* exp((1/2) * invert(dsig_rho^2) * (dgn - dghat)^2);
	return ddenom;			
}

RSV_skewFS::CalculateNumGamma(){
	if(m_isSkew != 1){
		return 1;
	}
	decl dloghat, dgammahat;
	decl Func, drv1, drv2;
	dloghat = log(m_dgamma);
	MaxBFGS(Posterior_Gamma_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Gamma_Max, dloghat, &drv2);

	decl dsig_gamma = sqrt(invert(-drv2));
	decl dm_gamma = dloghat;
	decl dlogo = log(m_dgamma);
	decl dlogn = log(m_dgamma_ast);
	
	Posterior_Gamma_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Gamma_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl dnum = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_gamma^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_gamma^2) * (dlogn - dloghat)^2);
	dnum *= densNormal(dlogn, dm_gamma, dsig_gamma^2);
	dnum *= 1/m_dgamma_ast;
	return dnum;
}

RSV_skewFS::CalculateDenomGamma(){
	if(m_isSkew != 1){
		return 1;
	}
	decl dloghat, dgammahat;
	decl Func, drv1, drv2;
	dloghat = log(m_dgamma_ast);
	MaxBFGS(Posterior_Gamma_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Gamma_Max, dloghat, &drv2);

	decl dsig_gamma = sqrt(invert(-drv2));
	decl dm_gamma = dloghat;
	decl dlogo = log(m_dgamma_ast);
	decl dlogn = log(m_dgamma);
	
	Posterior_Gamma_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Gamma_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl ddenom = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_gamma^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_gamma^2) * (dlogn - dloghat)^2);
	return ddenom;
}

RSV_skewFS::CalculateNumSigeta(){
	decl dloghat, dsighat;
	decl Func, drv1, drv2;
	dloghat = log(m_dsigeta);
	MaxBFGS(Posterior_SigEta_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_SigEta_Max, dloghat, &drv2);

	decl dsig_sigeta = sqrt(invert(-drv2));
	decl dm_sigeta = dloghat;
	decl dlogo = log(m_dsigeta);
	decl dlogn = log(m_dsigeta_ast);
	Posterior_SigEta_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_SigEta_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl dnum = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_sigeta^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_sigeta^2) * (dlogn - dloghat)^2);
	dnum *= densNormal(dlogn, dm_sigeta, dsig_sigeta^2);
	dnum *= 1/m_dsigeta_ast;
	return dnum;
}

RSV_skewFS::CalculateDenomSigeta(){
	decl dloghat, dsighat;
	decl Func, drv1, drv2;
	dloghat = log(m_dsigeta_ast);
	MaxBFGS(Posterior_SigEta_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_SigEta_Max, dloghat, &drv2);

	decl dsig_sigeta = sqrt(invert(-drv2));
	decl dm_sigeta = dloghat;
	decl dlogo = log(m_dsigeta_ast);
	decl dlogn = log(m_dsigeta);
	Posterior_SigEta_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_SigEta_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;	
	decl ddenom = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_sigeta^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_sigeta^2) * (dlogn - dloghat)^2);
	return ddenom;
}

RSV_skewFS::CalculateNumNu(){
	if(m_isT != 1){
		return 1;
	}
	decl dloghat, dnuhat;
	decl Func, drv1, drv2;
	dloghat = log(m_dnu - 2);
	MaxBFGS(Posterior_Nu_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Nu_Max, dloghat, &drv2);
	decl dsig_nu = invert(sqrt(-drv2));
	decl dm_nu = dloghat;
	decl dlogo = log(m_dnu - 2);
	decl dlogn = log(m_dnu_ast - 2);
	Posterior_Nu_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Nu_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;
	decl dnum = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_nu^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_nu^2) * (dlogn - dloghat)^2);
	dnum *= densNormal(dlogn, dm_nu, dsig_nu^2);
	dnum *= 1/(m_dnu_ast - 2);

	return dnum;
}

RSV_skewFS::CalculateDenomNu(){
	if(m_isT != 1){
		return 1;
	}
	decl dloghat, dnuhat;
	decl Func, drv1, drv2;
	dloghat = log(m_dnu_ast - 2);
	MaxBFGS(Posterior_Nu_Max, &dloghat, &Func, 0, 1);
	Num2Derivative(Posterior_Nu_Max, dloghat, &drv2);
	decl dsig_nu = invert(sqrt(-drv2));
	decl dm_nu = dloghat;
	decl dlogo = log(m_dnu_ast - 2);
	decl dlogn = log(m_dnu - 2);
	Posterior_Nu_Max(dlogn, &Func, 0, 0);
	decl dln = Func;
	Posterior_Nu_Max(dlogo, &Func, 0, 0);
	decl dlo = Func;
	decl ddenom = exp(dln - dlo)
				* exp(-(1/2) * invert(dsig_nu^2) * (dlogo - dloghat)^2)
				* exp((1/2) * invert(dsig_nu^2) * (dlogn - dloghat)^2);
	return ddenom;
}

RSV_skewFS::CalculateML(const nsim, const nburn){
	
	decl dm_temp, dsig_temp;
	decl dn1_temp, dS1_temp;
	decl dfrac;
	decl l, m;
	
	decl mnum = zeros(nsim, 8);
	decl mdenom = zeros(nsim, 8);

	//MH
	//Phi, Rho, Gamma, SigEta, Nu
	//Gibbs
	//Mu, Xi, SigU
	//Latent
	//H

	//Mu
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SampleMu();
		SamplePhi();
		SampleRho();
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		dsig_temp = (m_dsig0_mu^(-2) + m_dsigeta^(-2) * (1 - m_dphi^2 + (m_cT - 1) * (1 - m_dphi)^2 * (1 - m_drho^2)^(-1) ) )^(-1/2);
		dm_temp = dsig_temp^2 * (m_dsig0_mu^(-2) * m_dm0_mu
								+ m_dsigeta^(-2) *
									(  (1 - m_dphi^2) * m_vh[][0]
									+ ((1 - m_dphi) / (1 - m_drho^2))
										* sumr(m_vh[][1:(m_cT-1)] - m_dphi * m_vh[][0:(m_cT-2)]
											- m_drho * m_dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)]
						      	          	)
									)
					  		  );
		if(l>=0){
			mnum[l][0] = densNormal(m_dmu_ast, dm_temp, dsig_temp^2);
			mdenom[l][0] = 1;
		}
	}
	println("Mu complete");

	//Phi
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SamplePhi();
		SampleRho();
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		dsig_temp = sqrt( (m_dsigeta^2 * (1 - m_drho^2)) / (m_drho^2 * (m_vh[][0] - m_dmu)^2 + sumr((m_vh[][1:(m_cT-1)] - m_dmu).^2) ) );
		dm_temp = (dsig_temp^2 / (m_dsigeta^2 * (1 - m_drho^2)) )
						* sumr(
							(m_vh[][1:(m_cT-1)] - m_dmu - m_drho .* m_dsigeta .* exp(-m_vh[][0:(m_cT-2)]/2) .* m_my[][0:(m_cT-2)]) 
							.* (m_vh[][0:(m_cT-2)] - m_dmu)
							);
		dfrac = ((1 - m_dphi_ast)^(m_da0_phi - (1/2)) * (1 + m_dphi_ast)^(m_db0_phi - (1/2)))
				 / (1 - m_dphi)^(m_da0_phi - (1/2)) * (1 + m_dphi)^(m_db0_phi - (1/2));
		if(l>=0){
			mnum[l][1] = dfrac * densNormal(m_dphi_ast, dm_temp, dsig_temp^2);
		}
	}
	for(m=-nburn;m<nsim;m++){
		m_dphi = m_dphi_ast;
		SampleRho();
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		SamplePhi(1);
		dfrac = ((1 - m_dphi)^(m_da0_phi - (1/2)) * (1 + m_dphi)^(m_db0_phi - (1/2)))
				 / (1 - m_dphi_ast)^(m_da0_phi - (1/2)) * (1 + m_dphi_ast)^(m_db0_phi - (1/2));
		if(m>=0){
			mdenom[m][1] = dfrac;
		}
	}
	println("Phi complete");

	//Rho
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SampleRho();
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		if(l>=0){
			mnum[l][2] = CalculateNumRho();
		}
	}
	for(m=-nburn;m<nsim;m++){
		m_drho = m_drho_ast;
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		SampleRho(1);
		if(m>=0){
			mdenom[m][2] = CalculateDenomRho();
		}
	}
	println("Rho complete");

	//Gamma
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SampleGamma();
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		if(l>=0){
			mnum[l][3] = CalculateNumGamma();
		}
	}
	for(m=-nburn;m<nsim;m++){
		m_dgamma = m_dgamma_ast;
		SampleSigeta2();
		SampleXi();
		if(m_isRSV == 1){
			SampleSigu();
			SampleNu();
		}
		SampleH();
		SampleGamma(1);
		if(m>=0){
			mdenom[m][3] = CalculateDenomGamma();
		}
	}
	println("Gamma complete");

	//Sigeta
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SampleSigeta2();
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		if(l>=0){
			mnum[l][4] = CalculateNumSigeta();
		}
	}
	for(m=-nburn;m<nsim;m++){
		m_dsigeta = m_dsigeta_ast;
		if(m_isRSV == 1){
			SampleXi();
			SampleSigu();
		}
		SampleNu();
		SampleH();
		SampleSigeta2(1);
		if(m>=0){
			mdenom[m][4] = CalculateDenomSigeta();
		}
	}
	println("Sigeta complete");

	//Xi
	if(m_isRSV == 1){
		SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
		for(l=-nburn;l<nsim;l++){
			SampleXi();
			SampleSigu();
			SampleNu();
			SampleH();
			dsig_temp = sqrt(invert(m_dsig0_xi^(-2) + m_cT * m_dsigu^(-2)));
			dm_temp = dsig_temp^2 * (m_dsig0_xi^(-2) * m_dm0_xi + m_dsigu^(-2) * sumr(m_mx - m_dpsi * m_vh));
			if(l>=0){
				mnum[l][5] = densNormal(m_dxi_ast, dm_temp, dsig_temp^2);
				mdenom[l][5] = 1;
			}
		}
		println("Xi complete");
	}else{
		println("Xi skipped");
	}
	//Sigu
	if(m_isRSV == 1){
		SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
		for(l=-nburn;l<nsim;l++){
			SampleSigu();
			SampleNu();
			SampleH();
			dn1_temp = m_dn0_u + m_cT;
			dS1_temp = m_dS0_u + sumr((m_mx - m_dxi - m_dpsi * m_vh).^2);
			if(l>=0){
				mnum[l][6] = exp(logdensIGa(m_dsigu_ast^2, dn1_temp/2, dS1_temp/2));
				mdenom[l][6] = 1;
			}
		}
		println("Sigu complete");
	}else{
		println("Sigu skipped");
	}

	//Nu
	SetParam(m_dmu_ast, m_dphi_ast, m_drho_ast, m_dgamma_ast, m_dsigeta_ast, m_dxi_ast, m_dsigu_ast, m_dnu_ast, m_vh_ast, m_vlambda);
	for(l=-nburn;l<nsim;l++){
		SampleNu();
		SampleH();
		if(l>=0){
			mnum[l][7] = CalculateNumNu();
		}
	}
	for(m=-nburn;m<nsim;m++){
		m_dnu = m_dnu_ast;
		SampleH();
		SampleNu(1);
		if(m>=0){
			mdenom[m][7] = CalculateDenomNu();
		}
	}
	println("Nu complete");

	println("Numerator = ", mnum);
	println("Denominator = ",mdenom);
	println("NumMean = ", meanc(mnum)');
	println("DenomMean = ", meanc(mdenom)');
	
	decl LogPostDens = sumc(log(meanc(mnum)' ./ meanc(mdenom)'));
	decl PostvecDens = zeros(nsim, sizec(mnum) * 2);
	decl j;
	for(j=0;j<sizec(mnum);j++){
		PostvecDens[][2*j] = mnum[][j];
		PostvecDens[][2*j+1] = mdenom[][j];
	}
	decl PostbarvecDens = meanc(PostvecDens)';
	decl hatc = PostbarvecDens.^(-1);
	decl i;
	for(i=0;i<sizer(hatc);i++){
		hatc[i] *= (-1)^(i+2);
	}
	decl s;
	decl VarPost, Omega_s;
	VarPost = zeros(sizec(mnum) * 2, sizec(mnum) * 2);
	for(m=0;m<nsim;m++){
		VarPost += (1/nsim) * ( PostvecDens[m][] - PostbarvecDens' )' * ( PostvecDens[m][] - PostbarvecDens' );
	}
	for(s=1;s<=m_BM;s++){
		Omega_s = zeros(sizec(mnum) * 2, sizec(mnum) * 2);
		for(m=s;m<nsim;m++){
			Omega_s += ( PostvecDens[m][] - PostbarvecDens' )' * ( PostvecDens[m-s][] - PostbarvecDens' ) ;
		}
		Omega_s *= (1/nsim);
		VarPost	+= (1 - (s/(m_BM + 1)) * (Omega_s + Omega_s'));
	}
	VarPost *= (1/nsim);

	println("posterior density = ", LogPostDens);
	println("numerical standard error of estimation of posterior density = ", sqrt(hatc' * VarPost * hatc));
	//savemat("post_hbar.csv", PostbarvecDens);
	//savemat("post_hs.csv", PostvecDens);
	return LogPostDens;
}

RSV_skewFS::SearchMin(const vpi, const du){
	decl i;
	decl vpi_cum = zeros(vpi);
	for(i=0;i<sizer(vpi);i++){
		vpi_cum[i] = sumc(vpi[0:i]);
	}
	return sizer(vpi) - sizer(selectr((vpi_cum .< du), 0));
}

RSV_skewFS::WeightToRatio(const vw){
	decl vpi = vw ./ sumc(vw);
	return vpi;
}

RSV_skewFS::CalculateLikelihood(const nI){
	decl i;//Particle id
	decl t;
	decl j, nJ;
	decl vw, vpi;
	decl vw2, vpi2;
	decl mw, mpi;
	decl vw_rec;
	decl vh_current, vh_fwd, vh_temp;
	decl mh_current, mh_fwd, mh_temp;
	decl vg;
	decl mg;
	decl dmu_fwd;
	decl du;

	println(sprint("The Number of Particles: ", nI));
	nJ = 1;//10//R = I * J
	println(sprint("J = ", nJ));

	decl k;
	decl nrep = 1;//10;
	decl likelihood_box = zeros(nrep, 1);
	
	vw = vw2 = zeros(nI, 1);
	vpi = vpi2 = zeros(nI, 1);
	mw = zeros(nI, nJ);
	mpi = zeros(nI, nJ);
	vw_rec = zeros(1, m_cT);

	vh_current = vh_fwd = zeros(nI, 1);
	mh_current = mh_fwd = zeros(nI, nJ);
	vg = zeros(nI, 1);
	mg = zeros(nI, nJ);

	decl vg_picked = zeros(nI, 1);
	decl vpi_picked = zeros(nI, 1);
	decl pick_index;
	
	decl progress_bar = 0;
	println("Progress");
	println("--|--|--|--|--|--|--|--|--|--|");
	decl ddens_temp;

	for(k=0;k<nrep;k++){
		for(t=0;t<m_cT-1;t++){
		if(100 * (t / m_cT) >= (progress_bar + 1) * 5){
				print("=");
				progress_bar += 1;
				if(imod(progress_bar, 2) == 0 && progress_bar != 0){
					print("|");
				}
			}
			for(i=0;i<nI;i++){	
				if(t == 0){
					mh_current[i][] = ones(1, nJ) * (m_dmu + sqrt(m_dsigeta^2 / (1 - m_dphi^2)) * rann(1, 1));
					for(j=0;j<nJ;j++){
						if(m_isT == 1 && m_isSkew == 1){
							ddens_temp = densFSST(m_my[][0] * exp(-mh_current[i][j]/2), m_dgamma, m_dnu);
						}
						if(m_isT == 0 && m_isSkew == 1){
							ddens_temp = densFSSN(m_my[][0] * exp(-mh_current[i][j]/2), m_dgamma);
						}
						if(m_isT == 1 && m_isSkew == 0){
							ddens_temp = densFSST(m_my[][0] * exp(-mh_current[i][j]/2), 1, m_dnu);
						}
						if(m_isT == 0 && m_isSkew == 0){
							ddens_temp = densNormal(m_my[][0] * exp(-mh_current[i][j]/2), 0, 1);
						}
						mw[i][j] = ddens_temp  * exp(-mh_current[i][j]/2)
											   * densNormal(m_mx[][0] - m_dxi - m_dpsi*mh_current[i][j], 0, m_dsigu^2);
					}
					vw_rec[][0] = sumc(sumr(mw)) / (nI * nJ);
				}else{
					mh_current = mh_fwd;
				}
			}
			mpi = shape(WeightToRatio(vec(mw)), nI, nJ);
			for(i=0;i<nI;i++){
				for(j=0;j<nJ;j++){
					dmu_fwd = m_dmu + m_dphi * (mh_current[i][j] - m_dmu) + m_drho * m_dsigeta * exp(-mh_current[i][j]/2) * m_my[][t];
					if(m_isT == 1 && m_isSkew == 1){
						ddens_temp = densFSST(m_my[][t+1] * exp(-dmu_fwd/2), m_dgamma, m_dnu);
					}
					if(m_isT == 0 && m_isSkew == 1){
						ddens_temp = densFSSN(m_my[][t+1] * exp(-dmu_fwd/2), m_dgamma);
					}
					if(m_isT == 1 && m_isSkew == 0){
						ddens_temp = densFSST(m_my[][t+1] * exp(-dmu_fwd/2), 1, m_dnu);
					}
					if(m_isT == 0 && m_isSkew == 0){
						ddens_temp = densNormal(m_my[][t+1] * exp(-dmu_fwd/2), 0, 1);
					}
					mg[i][j] = ddens_temp * exp(-dmu_fwd/2)
								* densNormal(m_mx[][t+1] - m_dxi - m_dpsi*dmu_fwd, 0, m_dsigu^2) * mpi[i][j];
				}
			}
			mg = shape(WeightToRatio(vec(mg)), nI, nJ);
			//Sampling from g
			mh_temp = mh_current;
			for(i=0;i<nI;i++){
				pick_index = ranmultinomial(1, vec(mg))';
				vh_current[i] = vec(mh_temp)' * pick_index;
				vg_picked[i] = vec(mg)' * pick_index;
				vpi_picked[i] = vec(mpi)' * pick_index;
				mh_fwd[i][] = m_dmu + m_dphi * (vh_current[i] - m_dmu) + m_drho * m_dsigeta * exp(-vh_current[i]/2) * m_my[][t]
								+ sqrt((1 - m_drho^2)) * m_dsigeta * rann(1, nJ);
			}

			//Calculate Weight
			for(i=0;i<nI;i++){
				if(m_isT == 1 && m_isSkew == 1){
					ddens_temp = densFSST(m_my[][t+1] * exp(-mh_fwd[i][]/2), m_dgamma, m_dnu);
				}
				if(m_isT == 0 && m_isSkew == 1){
					ddens_temp = densFSSN(m_my[][t+1] * exp(-mh_fwd[i][]/2), m_dgamma);
				}
				if(m_isT == 1 && m_isSkew == 0){
					ddens_temp = densFSST(m_my[][t+1] * exp(-mh_fwd[i][]/2), 1, m_dnu);
				}
				if(m_isT == 0 && m_isSkew == 0){
					ddens_temp = densNormal(m_my[][t+1] * exp(-mh_fwd[i][]/2), 0, 1);
				}
				mw[i][] = (ddens_temp .* exp(-mh_fwd[i][]/2)
						.* ((m_isRSV == 1) ? densNormal(m_mx[][t+1] - m_dxi - m_dpsi*mh_fwd[i][], zeros(1, nJ), m_dsigu^2 * ones(1, nJ)) : 1)
						.* vpi_picked[i]) ./ vg_picked[i]
							;
			}
			vw_rec[][t+1] = sumr(sumc(mw))/(nI * nJ);
		}
		println("likelihood = ", sumr(log(vw_rec)));
		likelihood_box[k][] = sumr(log(vw_rec));
	}
	savemat("likelihood.csv", likelihood_box);
	return sumr(log(vw_rec));
}

RSV_skewFS::CalculateLikelihood_PF(const nI){
	decl i;//Particle id
	decl t;
	decl j, nJ;
	decl vw, vpi;
	decl vw2, vpi2;
	decl mw, mpi;
	decl vw_rec;
	decl vh_current, vh_fwd, vh_temp;
	decl mh_current, mh_fwd, mh_temp;
	decl vg;
	decl mg;
	decl dmu_fwd;
	decl du;

	nJ = 10;//R = I * J
	
	vw = vw2 = zeros(nI, 1);
	vpi = vpi2 = zeros(nI, 1);
	mw = zeros(nI, nJ);
	mpi = zeros(nI, nJ);
	vw_rec = zeros(1, m_cT);

	vh_current = vh_fwd = zeros(nI, 1);
	mh_current = mh_fwd = zeros(nI, nJ);
	vg = zeros(nI, 1);
	mg = zeros(nI, nJ);

	decl progress_bar = 0;
	println("Progress");
	println("--|--|--|--|--|--|--|--|--|--|");

	for(t=0;t<m_cT-1;t++){
	//println("t=", t);
	if(100 * (t / m_cT) >= (progress_bar + 1) * 5){
			print("=");
			progress_bar += 1;
			if(imod(progress_bar, 2) == 0 && progress_bar != 0){
				print("|");
			}
		}
		for(i=0;i<nI;i++){	
			if(t == 0){
				mh_current[i][] = ones(1, nJ) * (m_dmu + sqrt(m_dsigeta^2 / (1 - m_dphi^2)) * rann(1, 1));
				for(j=0;j<nJ;j++){
					mw[i][j] = densFSST(m_my[][0] * exp(-mh_current[i][j]/2), m_dgamma, m_dnu) * exp(-mh_current[i][j]/2)
										   * densNormal(m_mx[][0] - m_dxi - m_dpsi*mh_current[i][j], 0, m_dsigu^2);
				}
				vw_rec[][0] = sumc(sumr(mw)) / (nI * nJ);
			}else{
				mh_current = mh_fwd;
			}
		}
		mpi = shape(WeightToRatio(vec(mw)), nI, nJ);
		//Sampling from pi
		mh_temp = mh_current;
		for(i=0;i<nI;i++){
			vh_current[i] = vec(mh_temp)' * ranmultinomial(1, vec(mpi))';
			mh_fwd[i][] = m_dmu + m_dphi * (vh_current[i] - m_dmu) + m_drho * m_dsigeta * exp(-vh_current[i]/2) * m_my[][t]
							+ sqrt((1 - m_drho^2)) * m_dsigeta * rann(1, nJ);
		}
		
		for(i=0;i<nI;i++){
			mw[i][] = (densFSST(m_my[][t+1] * exp(-mh_fwd[i][]/2), m_dgamma, m_dnu) .* exp(-mh_fwd[i][]/2)
					.* densNormal(m_mx[][t+1] - m_dxi - m_dpsi*mh_fwd[i][], zeros(1, nJ), m_dsigu^2 * ones(1, nJ)))
						;
		}
		vw_rec[][t+1] = sumr(sumc(mw))/(nI * nJ);
	}
	println("likelihood = ", sumr(log(vw_rec)));
}

RSV_skewFS::CalculatePrior(){
	decl dprior = 0;

	dprior += log(densNormal(m_dmu, m_dm0_mu, m_dsig0_mu^2));
	if(m_isRSV == 1){
		dprior += log(densNormal(m_dxi, m_dm0_xi, m_dsig0_xi^2));
	}
	dprior += log(densbeta((m_dphi+1)/2, m_da0_phi, m_db0_phi) * (1/2));
	dprior += log(densbeta((m_drho+1)/2, m_da0_rho, m_db0_rho) * (1/2));
	dprior += logdensIGa(m_dsigeta^2, m_dn0_eta/2, m_dS0_eta/2);
	if(m_isRSV == 1){
		dprior += logdensIGa(m_dsigu^2, m_dn0_u/2, m_dS0_u/2);
	}
	if(m_isRSV == 1 && m_isPsi == 1){
		dprior += log(densNormal(m_dpsi, m_dm0_psi, m_dsig0_psi^2));
	}
	if(m_isSkew == 1){
		dprior += logdensGa(m_dgamma, m_dn0_gamma, m_dS0_gamma);
	}
	if(m_isT == 1){
		dprior += logdensGa(m_dnu, m_dn0_nu, m_dS0_nu);
	}
	return dprior;
}

RSV_skewFS::SimObs(){
	decl t;
	decl veps = zeros(1, m_cT);
	decl veta = zeros(1, m_cT);
	decl vu = zeros(1, m_cT);
	
	//initial state
	m_vh[][0] = m_dmu + sqrt(m_dsigeta^2 / (1 - m_dphi^2)) * rann(1, 1);
	
	for(t=0;t<m_cT;t++){
		m_vlambda[][t] = invert(rangamma(1, 1, m_dnu/2, m_dnu/2));
		if(m_isT == 1){
			veps[][t] = ranFS(m_dgamma, m_dnu, m_vlambda[][t]);
		}else{
			veps[][t] = ranFSN(m_dgamma);
		}
		veta[][t] = m_drho * m_dsigeta * veps[][t] + sqrt((1 - m_drho^2) * m_dsigeta^2) * rann(1, 1);
		vu[][t] = m_dsigu * rann(1, 1);
		if(t < m_cT - 1){
			m_vh[][t+1] = m_dmu + m_dphi * (m_vh[][t] - m_dmu) + veta[][t];
		}
		m_mx[][t] = m_dxi + m_dpsi * m_vh[][t] + vu[][t];
		m_my[][t] = veps[][t] * exp(m_vh[][t]/2);
	}																	  
	savemat("eps.csv", veps');
	savemat("eta.csv", veta');
	savemat("h.csv", m_vh');
	savemat("x.csv", m_mx');
	savemat("y.csv", m_my');
}
