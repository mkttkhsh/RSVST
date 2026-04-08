/*----------------------------------------------------------
                        SVAZ.ox
----------------------------------------------------------**
**  Reference:
**    Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).
**    Realized Stochastic Volatility Model with Skew-t Distributions for
**    Improved Volatility and Quantile Forecasting. arXiv:2501.13179
**    https://arxiv.org/abs/2501.13179
*/
											 
/*
**  CLASS FILE
**
**  MCMC estimation for stochastic volatility model
**  with skew t-distribution discussed in Azzalini (1985)
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<SVAZ.h>

SVAZ::SVAZ(const vdata)	//vdata=0 means simulation study
{												   
	m_iRun = ismatrix(vdata) ? 1 : 0; //0:simulation, 1:empirical
	m_vy = m_vyorg = vdata;		  

	if(m_iRun){
		m_ns = sizer(vdata);
		m_nKm = int(m_ns/25)-1;	//# of blocks in multi-move sampler
	}

	m_fls = 1;			//model specification: skew
	m_flt = 1;			//model specification: t
	m_vini = 0;			//initial value
	m_ncount = 1000;	//counter of iteration
	m_flML = 0;			//flag for ML estimation
	m_flf  = 0;			//flag for volatility and quantile forecasts
}

/*
** SetFilename set main name for output files
*/
SVAZ::SetFilename(const sfilename)
{
	m_sfilename = sfilename;
}

/*
** SetDirectory set directory for output files
*/
SVAZ::SetDirectory(const sdirectory)
{
	m_sdirectory = sdirectory;
}

/*
** SetDist set distribution of return
*/
SVAZ::SetDist(const fls, const flt)
{
	m_fls = fls; //0:Symmetric; 1:Skew
	m_flt = flt; //0:Normal; 1:Student's t
}

/*
** SetRanseed set ranseed for random generater
*/
SVAZ::SetRanseed(const irs)
{
	ranseed(irs);
}

/*
** SetSampSigRho set sampling method for sigma and rho
*/
SVAZ::SetSampSigRho(const flsigrho)
{
	m_flsigrho = flsigrho; //0:separate sampling; 1:joint sampling
}

/*
** SetSampNu set sampling method for nu
*/
SVAZ::SetSampNu(const flnuarmh)
{
	m_flnuarmh = flnuarmh; //0:M-H sampling; 1:A-R M-H sampling
}

/*
** SetMultiMove set sampling method for latent volatilities h
*/
SVAZ::SetMultiMove(const flmultimove)
{
	m_flmultimove = flmultimove; //0:single-move; 1:multi-move
}

/*
** SetBlocks set # of blocks in multi-move sampler
*/
SVAZ::SetBlocks(const nKm)
{
	m_nKm = nKm-1;
}

/*
** SetCount set counter of iteration
*/
SVAZ::SetCount(const ncount)
{
	m_ncount = ncount;
}

/*
** SetForecast set VaR and ES backtest
*/
SVAZ::SetForecast(const iForecast)
{
	m_flf = iForecast; //0:off; 1:on
}

/*
**  set true parameter for simulation
*/
SVAZ::SimObs(const ns, const dphi, const dsig, const drho, const dmu, const ddelta_a, const dnu_a)
{
	decl mSigsq, me, vh, vz0, ddelta, dc, veps, dcov, dcorr, dskew, dkurt;
	decl dnu, vlambda, dmu_lambda, vy, dvar_y, dskew_y, dkurt_y, t;
	
	m_ns = ns;

	// Generate log-volatilities
	mSigsq = choleski(1 ~ drho*dsig | drho*dsig ~ dsig^2);
	me = (mSigsq * rann(2,ns))';
	vh = zeros(ns, 1);
	vh[0] = dmu + rann(1, 1) * dsig / sqrt(1-dphi^2);
	for(t=0; t<ns-1; t++){
		vh[t+1] = dmu + dphi * (vh[t] - dmu) + me[t][1];
	}

	// Generate errors in daily returns
	vz0 = zeros(ns,1);
	if(m_fls){
		ddelta = ddelta_a;
		for(t=0; t<ns; t++){
			while(vz0[t]<=0){
				vz0[t] = rann(1,1);
			}
		}
	}
	else ddelta = 0;
	dc = sqrt(2 / M_PI);
	veps = (ddelta.*(vz0-dc) + sqrt(1-ddelta^2).*me[][0]) ./ sqrt(1 - dc^2*ddelta^2);
	dcov = drho * dsig * sqrt((1-ddelta^2)/(1-dc^2*ddelta^2)); //theoretical covariance between epsion and eta
	dcorr = drho * sqrt((1-ddelta^2)/(1-dc^2*ddelta^2)); //theoretical correlation between epsion and eta
	dskew = 0.5 * (4 - M_PI) * dc^3 * ddelta^3 * (1-dc^2*ddelta^2)^(-1.5); //theoretical skewness of epsilon
	dkurt = 3 + 2 * (M_PI - 3) * dc^4 * ddelta^4 / (1 - dc^2 * ddelta^2)^2; //theoretical kurtosis of epsilon
	dvar_y = exp(dmu + 0.5 * dsig^2 / (1 - dphi^2)); //theoretical variance of returns
	
	if(m_flt){ // (skew) t
		dnu = dnu_a;
		vlambda = 1 ./ rangamma(ns, 1, dnu/2, dnu/2);
		dmu_lambda = dnu / (dnu - 2);
		dskew_y = dskew * (dnu - 2)^1.5 * gammafact(0.5*dnu - 0.5) / (sqrt(2) * (dnu - 3) * gammafact(0.5*dnu))
					* exp(3/8 * dsig^2 / (1 - dphi^2)); //theoretical skewness of returns
		dkurt_y = dkurt * (dnu - 2) / (dnu - 4) * exp(dsig^2 / (1 - dphi^2)); //theoretical kurtosis of retuns
	}
	else{ // (skew) normal
		dnu = M_INF_POS;
		vlambda = ones(ns, 1);
		dmu_lambda = 1;
		dskew_y = dskew * exp(3 * dsig^2 / (8 * (1 - dphi^2))); //theoretical skewness of returns
		dkurt_y = dkurt * exp(dsig^2 / (1 - dphi^2)); //theoretical kurtosis of retuns
	}

	// Generate daily returns
	vy = veps .* sqrt(vlambda/dmu_lambda) .* exp(vh/2);
	
	m_vy = m_vyorg = vy;
	m_vh_t = vh;
	m_vlambda_t = vlambda;
	m_vz0_t = vz0;
	
	print("\n>>> Data is simulated. [SV-AZ-");
	if(m_fls)
		print("S");
	if(m_flt){
		print("T");
	}else{
		print("N");
	}
	print(" model]\n");
	print("Parameters:");
	println("%r", {"    phi", "    sigma", "    rho", "    mu", "    delta", "    nu"},
			"%cf", {"%9.2f"}, dphi | dsig | drho | dmu | ddelta | dnu);
	println("%r", {"E(eps)", "Var(eps)", "Skew(eps)", "Kurt(eps)", "E(y)", "Var(y)", "Skew(y)", "Kurt(y)",
					"E(eta)", "Var(eta)", "Corr(z,eta)", "Cov(eps,eta)", "Corr(eps,eta)"},
			"%c", {"True", "Sample"}, "%cf", {"%9.4f"},
			0 ~ meanc(veps) | 1 ~ varc(veps) | dskew ~ moments(veps,3)[3] | dkurt ~ moments(veps,4)[4] |
			0 ~ meanc(vy) | dvar_y ~ varc(vy) | dskew_y ~ moments(vy,3)[3] | dkurt_y ~ moments(vy,4)[4] |
			0 ~ meanc(me[][1]) | dsig^2 ~ varc(me[][1]) |
			drho ~ correlation(me)[0][1] | dcov ~ variance(veps~me[][1])[0][1] | dcorr ~ correlation(veps~me[][1])[0][1]);

	savemat(sprint(m_sdirectory, "/simdata_", m_sfilename, ".csv"), vy~vh~vlambda~vz0);
	savemat(sprint(m_sdirectory, "/simpar_", m_sfilename, ".csv"), dphi|dsig|drho|dmu|ddelta|dnu);

	m_nKm = int(ns/25);
}

/*
** SetIni set initial value for MCMC sampling
*/
SVAZ::SetIni(const dphi, const dsig, const drho, const dmu, const ddelta, const dnu, const vh)
{
	m_vini = dphi ~ dsig ~ drho ~ dmu ~ ddelta ~ dnu;
	m_vhini = vh;
}

/*
** SetPri set hyper parameters for priors
*/
SVAZ::SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0, const dad0, const dbd0,
	const dnn0, const dnS0, const flprint)
{
	m_vpri = dap0 ~ dbp0 ~ dn0 ~ dS0 ~ dar0 ~ dbr0 ~ dmm0 ~ dsm0 ~ dad0 ~ dbd0 ~ dnn0 ~ dnS0;
	if(flprint)
	{
		print("\n");
		println("Priors:");
		print("    (phi+1)/2   ~ Beta(", dap0, ", ", dbp0, ")\n");
		print("    sigma^2     ~ InvGamma(", dn0/2, ", ", dS0/2, ")\n");
		print("    (rho+1)/2   ~ Beta(", dar0, ", ", dbr0, ")\n");
		print("    mu          ~ N(", dmm0, ", ", dsm0^2, ")\n");
		if(m_fls)
		{
		print("    (delta+1)/2 ~ Beta(", dad0, ", ", dbd0, ")\n");
		}
		if(m_flt)
		{
		print("    nu          ~ Gamma(", dnn0, ", ", dnS0, ")\n");
		}
		print("\n");
		print("Moments of Priors:");
		println("%c", {"Mean", "Stdev"},
			  "%r", {"    phi", "    sigma^2", "    rho", "    nu", "    delta"}, "%cf", {"%9.2f"},
				fBetaMoment(dap0,dbp0,1)*2-1 ~ sqrt(fBetaMoment(dap0,dbp0,2))*2 |
				fInvGammaMoment(dn0/2,2/dS0,1) ~ sqrt(fInvGammaMoment(dn0/2,2/dS0,2)) |
				fBetaMoment(dar0,dbr0,1)*2-1 ~ sqrt(fBetaMoment(dar0,dbr0,2))*2 |
				fGammaMoment(dnn0,dnS0,1) ~ sqrt(fGammaMoment(dnn0,dnS0,2)) |
				fBetaMoment(dad0,dbd0,1)*2-1 ~ sqrt(fBetaMoment(dad0,dbd0,2))*2);
	}
}

/*
** fBetaMoment returns the first or second moment of Beta distribution
*/
SVAZ::fBetaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(da+db));
	else if(cm==2) return(da*db/((da+db)^2*(da+db+1)));
}

/*
** fGammaMoment returns the first or second moment of Gamma distribution
*/
SVAZ::fGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/db);
	else if(cm==2) return(da/db^2);
}

/*
** fInvGammaMoment returns the first or second moment of inverse-Gamma distribution
*/
SVAZ::fInvGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(db-1));
	else if(cm==2) return(db^2/(sqr(da-1)*(da-2)));
}

/*
** DrawData draw graph of data (and state if m_iRun=0)
*/
SVAZ::DrawData()
{
	Draw(0, m_vyorg[][0]');
	DrawTitle(0, "Daily return");
	if(!m_iRun){
		Draw(1, m_vh_t');
		DrawTitle(1, "Simulated log volatility ($h_t$)");
	}
	DrawDensity(1+!m_iRun, m_vyorg[][0]', 1, 1, 0, 1, 1);
	DrawTitle(1+!m_iRun, "Density of daily return");	
	DrawTitle(2+!m_iRun, "QQ-plot to normal density");	
	DrawDensity(3+!m_iRun, m_vyorg[][1]', 1, 1, 0, 1, 1);
	DrawTitle(3+!m_iRun, "Density of log realized volatility");	
	DrawTitle(4+!m_iRun, "QQ-plot to normal density");	
  	SaveDrawWindow(sprint(m_sdirectory, "/datadraw_",m_sfilename,".eps"));
  	CloseDrawWindow();
}


/*
** MCMC computes MCMC estimation
*/

SVAZ::MCMC(const nsim, const nburn)
{
	decl dphi, dsig, drho, dmu, ddelta, dnu, ncyc;
	decl mstate, vh, vz0, dc, dk_delta, vlambda, dmu_lambda, npmt, msamp, vsamph, vsampvol, msamph, vsampz0, vsamplambda;
	decl vsamp, iBm, iacf, aspar, aspar2, mpost, mposh, ascf, vaccptprob;
	decl vlogpri, vm, mV, vlogg, p, vp, dW, sign, dmgfp, dlogml;
	decl vsamph1, vsampy1, T, dlambda1, dz01, vprob, vvar, mfly, ves, vvolf, mpred, aspred;
	decl time;
	decl i, j, t;

	if(m_ncount>0) time = timer();
	
	/*--- Set Initial Values ---*/
	// Parameters
	if(!m_flML){
		if(ismatrix(m_vini)){
			dphi   = m_vini[0];
			dsig   = m_vini[1];
			drho   = m_vini[2];
			dmu    = m_vini[3];
			ddelta = m_vini[4];
			dnu    = m_vini[5];
		}
		else{
			dphi = 0.97;
			dsig = 0.2;
			drho = -0.3;
			dmu = 0.4;
			ddelta = -0.5;
			dnu = 10;
		}
	}
	else{
		dphi = m_vposm[0];
		dsig = m_vposm[1];
		drho = m_vposm[2];
		dmu  = m_vposm[3];
		ddelta = m_fls ? m_vposm[4] : 0;
		dnu    = m_flt ? m_vposm[4+m_fls] : M_INF_POS;
	}

	// Auxiliary varialbe: z_{0t}
	vz0 = zeros(m_ns, 1)-1;
	for(t=0; t<m_ns; t++){
		while(vz0[t]<=0){
			vz0[t] = rann(1,1);
		}
	}
	dc = sqrt(2 / M_PI);
	if(!m_fls){
		ddelta = 0;
		vz0 = ones(m_ns,1) * dc;
	}
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	

	// Auxiliary variable: lambda_t
	if(m_flt){ // (skew) t
		vlambda = 1 ./ rangamma(m_ns, 1, dnu/2, dnu/2);
		dmu_lambda = dnu / (dnu - 2);
	}
	else{ // normal
		dnu = M_INF_POS;
		vlambda = ones(m_ns, 1);
		dmu_lambda = 1;
	}

	// Log-volatilities							
	if(!m_flML){
		if(ismatrix(m_vhini)) vh = m_vhini;
		else vh = dmu + rann(m_ns, 1) * dsig;
	}
	else{
		mstate = loadmat(sprint(m_sdirectory, "/posh_",m_sfilename,".csv"));
		vh = mstate[][0];
		vz0 = mstate[][2];
		vlambda = mstate[][3];
	}
	if(!m_iRun){ // if simulation
		vh = m_vh_t;
		vlambda = m_vlambda_t;
		vz0 = m_vz0_t;
	}

	// One-day-ahead prediction
	if(m_flf) aspred = {"Log-volatility", "Return"};
		


	/*--- Set priors ---*/
	if(ismatrix(m_vpri)){
		m_dap0 = m_vpri[0];		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = m_vpri[1];
		m_dn0  = m_vpri[2];		//sigma^2 ~ InvGamma(n0/2, S0/2)
		m_dS0  = m_vpri[3];
		m_dar0 = m_vpri[4];			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = m_vpri[5];
		m_dmm0 = m_vpri[6];		//mu ~ N(mm0, sm0^2)
		m_dsm0 = m_vpri[7];	
		m_dad0 = m_vpri[8];		//(delta+1)/2 ~ Beta(ad0, bd0)
		m_dbd0 = m_vpri[9];
		m_dnn0  = m_vpri[10];		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = m_vpri[11];
	}
	else{
		m_dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = 1.5;
		m_dn0  = 5;		//sigma^2 ~ InvGamma(n0/2, S0/2)
		m_dS0  = 0.2;
		m_dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = 2;
		m_dmm0 = dmu;		//mu ~ N(mm0, sm0^2)
		m_dsm0 = 1;	
		m_dad0 = 1;		//(delta+1)/2 ~ Beta(ad0, bd0)
		m_dbd0 = 1;
		m_dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = 0.5;
	}

	/*--- Set sample option ---*/
	npmt = 4+m_fls+m_flt; //number of parameters

	if(!m_flML)	/*--------- parameter estimation --------*/
	{
		msamp = zeros(nsim, npmt);	//sample box
		vsamph = vsampvol = zeros(m_ns, 1);
		msamph = zeros(nsim, m_ns);
		vsampz0 = zeros(m_ns, 1);
		vsamplambda = zeros(m_ns, 1);
		vsamph1 = zeros(nsim, 1);
		vsampy1 = zeros(nsim, 1);
		m_cmp = m_cms = m_cmr = m_cmd = m_can = m_cmn = m_cml = m_cah = m_cmh = 0; //acceptance counter
	
		/*--- MCMC sampling ---*/
		if(m_ncount > 0) println("\n\nIteration:");
	
		/*----------- S A M P L I N G   S T A R T -----------*/
		for(m_k=-nburn ; m_k<nsim ; m_k++)
		{
			/*--- sampling mu ---*/
			dmu = SampMu(dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		
			/*--- sampling phi ---*/
			dphi = SampPhi(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		
			if(m_flsigrho){ // joint sampling
				/*--- sampling sigma and rho ---*/
				[dsig,drho] = SampSigRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			}
			else{ // separate sampling
				/*--- sampling rho ---*/
				drho = SampRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

				/*--- sampling sigma ---*/
				dsig = SampSigma(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			}
			
			/*--- sampling skewness parameter: delta ---*/
			if(m_fls) //skew distribution
				ddelta = SampDelta(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

			/*--- sampling nu ---*/
			if(m_flt) // if Student's t
				dnu = SampNu(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

			/*--- sampling z0 ---*/
			if(m_fls) // if skewed distribution
				vz0 = SampZ0(dmu, dphi, dsig, drho, ddelta, dnu, vh, vlambda, m_vy);

			/*--- sampling h ---*/
			if(m_flmultimove)
				vh = MultiMove(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			else
				vh = SingleMove(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

			/*--- sampling lambda ---*/
			if(m_flt) // if Student's t
				vlambda = SampLambda(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

	
			if(m_k >= 0)
			{
				/*--- storing sample ---*/
				vsamp = dphi ~ dsig ~ drho ~ dmu;
				if(m_fls)
					vsamp ~= ddelta;
				if(m_flt)
					vsamp ~= dnu;
				
				msamp[m_k][] = vsamp;
				vsamph += vh;
				vsampvol += exp(vh/2);
				msamph[m_k][] = vh';
				vsampz0 += vz0;
				vsamplambda += vlambda;

				/*--- one-day-ahead prediction ---*/
				if(m_flf)
				{
					// Student's t or skew-t
					if(m_flt){
						dlambda1 = 1 / rangamma(1, 1, dnu/2, dnu/2);
						dmu_lambda = dnu / (dnu - 2);
					}
					else{
						dlambda1 = 1;
						dmu_lambda = 1;
					}

					// skew-normal or skew-t
					if(m_fls){
						dz01 = rann(1,1);
						while(dz01 <= 0){
							dz01 = rann(1,1);
						}
						dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
					}
					else dz01 = dc;

					// one-day-ahead log-volatility, return, and log-RV
					T = m_ns-1;
					vsamph1[m_k] = dmu + dphi * (vh[T] - dmu)
								   + sqrt(dmu_lambda/vlambda[T]) * drho * dsig * exp(-vh[T]/2) * m_vy[T] / dk_delta
								   - ddelta * drho * dsig * (vz0[T]-dc) / sqrt(1-dc^2)
								   + dsig * sqrt(1-drho^2) * rann(1,1); // log-volatility
					vsampy1[m_k] = (ddelta * (dz01-dc) / sqrt(1-dc^2*ddelta^2) + dk_delta * rann(1,1))
									  * sqrt(dlambda1/dmu_lambda) * exp(vsamph1[m_k]/2); // return
				}
			}
		
			//print counter
			if(m_ncount>0 && !imod(m_k, m_ncount)) println(m_k);
		}		 
	    /*------------ S A M P L I N G   E N D --------------*/
	
		
		/*--- output result ---*/
		
		iBm = iacf = min(1000, nsim/2);	//bandwidth
	
		aspar =	{"phi", "sigma", "rho", "mu"};
		aspar2 = {"$\\phi$", "$\\sigma$", "$\\rho$", "$\\mu$"};
		if(m_fls){
			aspar ~= {"delta"};
			aspar2 ~= {"$\\delta$"};
		}
		if(m_flt){
			aspar ~= {"nu"};
			aspar2 ~= {"$\\nu$"};
		}
		aspar ~= {"h(500)"};
		aspar2 ~= {"$h_{500}$"};
		
		print("\n                       ",
			  "[ESTIMATION RESULT]"
			  "\n-----------------------------------------------------------------"
			  "\nParameter    Mean     Stdev      95%L      95%U   Geweke    Inef."  
			  "\n-----------------------------------------------------------------");
	
		mpost = zeros(npmt, 6);
		for(i=0 ; i<npmt ; i++)
			mpost[i][]
				= meanc(msamp[][i])
				~ sqrt(varc(msamp[][i]))
				~ quantilec(msamp[][i], <.025, .975>)'
				~ fGeweke(msamp[][i],iBm)
				~ fTsvar(msamp[][i], iBm) / varc(msamp[][i]);
	
		mposh = zeros(m_ns, 6);
		for(i=0 ; i<m_ns; i++)
			mposh[i][]
			= meanc(msamph[][i])
			~ sqrt(varc(msamph[][i]))
			~ quantilec(msamph[][i], <.025, .975>)'
			~ fGeweke(msamph[][i],iBm)
			~ fTsvar(msamph[][i], iBm) / varc(msamph[][i]);

		ascf = {"%11.4f", "%10.4f", "%10.4f",
				"%10.4f", "%9.3f", "%9.2f"};
		print("%r", aspar, "%cf", ascf, mpost);
		println("-------------------------------------------------------- RSVAZ.ox");
	
		print("Model: SV-AZ-");
		if(m_fls) print("S");
		if(m_flt) print("T");
		else print("N");
		if(m_iRun) print(" (empirical study: n=", m_ns, ")\n");
		else	   print(" (simulation study: n=", m_ns, ")\n");
		print("Iteration: ", nsim, "\n");

		print("\nAcceptance rate for parameter"
			  "\n phi:      ", "%5.1f", m_cmp/nsim*100, "%");
		if(m_flsigrho)
			print("\n sigma,rho:", "%5.1f", m_cms/nsim*100, "%");
		else{
			print("\n sigma:    ", "%5.1f", m_cms/nsim*100, "%",
				  "\n rho:      ", "%5.1f", m_cmr/nsim*100, "%");
		}
		print("\n delta:    ", "%5.1f", m_cmd/nsim*100, "%");
		vaccptprob = m_cmp/nsim ~ m_cms/nsim ~ m_cmr/nsim ~ m_cmd/nsim;
		if(m_flt){
			if(m_flnuarmh){
				print(
					  "\n nu A-R:   ", "%5.1f", nsim/m_can*100, "%",
					  "\n nu M-H:   ", "%5.1f", m_cmn/nsim*100, "%");
			}
			else print("\n nu:       ", "%5.1f", m_cmn/nsim*100, "%");
					  
			print("\n lambda:   ", "%5.1f", double(m_cml/(nsim*m_ns)*100), "%");
			vaccptprob ~= m_flnuarmh?nsim/m_can:-1 ~ m_cmn/nsim ~ m_cml/(nsim*m_ns);
		}

		if(m_flmultimove){
			print("\nAcceptance rate for multi-move sampler",
				  "\n A-R: ", "%5.1f", double((m_nKm+1)*nsim/m_cah*100), "%",
				  "\n M-H: ", "%5.1f", double(m_cmh/((m_nKm+1)*nsim)*100), "%");
			print("\n");
			vaccptprob ~= (m_nKm+1)*nsim / m_cah ~ m_cmh / ((m_nKm+1)*nsim);
		}
		else{
			print("\nAcceptance rate for single-move sampler"
				  "\n M-H: ", "%5.1f", double(m_cmh/(nsim*m_ns)*100), "%");
			print("\n");
			vaccptprob ~= m_cmh/(nsim*m_ns);
		}

		print("\nCorrelation");
		println("%r", aspar, "%c", aspar, "%8.2f", correlation(msamp));

		/*--- value-at-risk, expected shortfall, and volatility forecast ---*/
		if(m_flf)
		{
			vprob = <.005, .01, .05, .1>;
			vvar = quantilec(vsampy1, vprob)'; // value-at-risk
			mfly = vsampy1 .< vvar;
			ves = vsampy1 ' mfly ./ sumc(mfly); // expected shortfall
			vvolf = meanc(vsamph1) ~ quantilec(vsamph1); // mean and median of one-day-ahead log-volatilities
			print("\n\nOne-day-ahead quantile of returns");
			print("%r", {" Tail", " VaR", " ES"}, 
				  "%cf", {"%10.5f"}, vprob | vvar | ves);
			print("\nOne-day-ahead log-volatility");
			println("%r", {" Mean", " Median"}, "%cf", {"%8.5f"}, vvolf');
			savemat(sprint(m_sdirectory, "/var_",m_sfilename,".csv"), vvar);
			savemat(sprint(m_sdirectory, "/es_",m_sfilename,".csv"), ves);
			savemat(sprint(m_sdirectory, "/volf_",m_sfilename,".csv"), vvolf);
			savemat(sprint(m_sdirectory, "/h1_y_", m_sfilename, ".csv"), vsamph1 ~ vsampy1); //1-day ahead log-volatility and return
		}
	
		if(m_ncount > 0)
		{
			//draw MCMC result
			DrawAdjust(ADJ_AREAMATRIX, 3, npmt);
		 	for(i=0 ; i<npmt ; i++){
		 		DrawCorrelogram(i, msamp[][i]', {}, iacf);
		    	DrawTitle(i, aspar2[i]);
		    	DrawMatrix(i+npmt, thinr(msamp[][i], 500)',
						   {""}, 1, nsim/500);
		    	DrawTitle(i+npmt, aspar2[i]);
		    	DrawDensity(i+2*npmt, msamp[][i]', {""}, 1);
		    	DrawTitle(i+2*npmt, aspar2[i]);
		    } 
		  	SaveDrawWindow(sprint(m_sdirectory, "/samp_",m_sfilename,".eps"));
		  	CloseDrawWindow();
	
			savemat(sprint(m_sdirectory, "/samp_",m_sfilename,".csv"), msamp);
	
			//draw state variable
			DrawMatrix(0, m_vy', "", 1, 1);
			DrawTitle(0, "Data");
			if(m_iRun){
				DrawMatrix(1, (vsamph/nsim)', "", 1, 1);
				DrawTitle(1, "Posterior mean of state");
			}	
			else{
				DrawMatrix(1, (m_vh_t ~ (vsamph/nsim))',{"True", "PosMean"}, 1, 1);
				DrawTitle(1, "True and posterior mean of state");
				DrawMatrix(2, (m_vz0_t ~ (vsampz0/nsim))',{"True", "PosMean"}, 1, 1);
				DrawTitle(2, "True and posterior mean of $z_0$");
				DrawMatrix(3, (m_vlambda_t ~ (vsamplambda/nsim))',{"True", "PosMean"}, 1, 1);
				DrawTitle(3, "True and posterior mean of $\lambda$");
			}
		    SaveDrawWindow(sprint(m_sdirectory, "/state_",m_sfilename,".eps"));
		    CloseDrawWindow();

			//draw p-value of Geweke's CD test and inefficiency factor of h
			DrawMatrix(0, mposh[][4]', "", 1, 1);
			DrawTitle(0, "p-value of Geweke's CD test for state");
			DrawMatrix(1, mposh[][5]', "", 1, 1);
			DrawTitle(1, "Inefficiency factor of state");
			SaveDrawWindow(sprint(m_sdirectory, "/geweke_inef_state_", m_sfilename, ".eps"));
		    CloseDrawWindow();

			//save sample box
			savemat(sprint(m_sdirectory, "/samp_", m_sfilename, ".csv"), msamp);
			savemat(sprint(m_sdirectory, "/samph_", m_sfilename, ".csv"), msamph);

			//draw one-day-ahead prediction result
			if(m_flf){
				mpred = vsamph1 ~ vsampy1;
				DrawAdjust(ADJ_AREAMATRIX, 3, 2);
			 	for(i=0 ; i<2; i++){
			 		DrawCorrelogram(i, mpred[][i]', {}, iacf);
			    	DrawTitle(i, aspred[i]);
			    	DrawMatrix(i+2, thinr(mpred[][i], 500)',
							   {""}, 1, nsim/500);
			    	DrawTitle(i+2, aspred[i]);
			    	DrawDensity(i+4, mpred[][i]', {""}, 1);
			    	DrawTitle(i+4, aspred[i]);
			    } 
			  	SaveDrawWindow(sprint(m_sdirectory, "/prediction_",m_sfilename,".eps"));
			  	CloseDrawWindow();
			}
		}
		savemat(sprint(m_sdirectory, "/posstat_", m_sfilename, ".csv"), mpost);
		savemat(sprint(m_sdirectory, "/posh_", m_sfilename, ".csv"), vsamph/nsim ~ vsampvol/nsim ~ vsampz0/nsim ~ vsamplambda/nsim);
		savemat(sprint(m_sdirectory, "/posstath_", m_sfilename, ".csv"), mposh);
		savemat(sprint(m_sdirectory, "/aprob_", m_sfilename, ".csv"), vaccptprob);
	}	//end of parameter estimation
	else	/*------------------------- ML estimation ---*/
	{
		/*--- set variables ---*/
		ncyc = 5 + 2*m_fls + 2*m_flt;	//# of cycle for reduced run
		m_mmls = zeros(nsim, ncyc);
	 
		/*--- MCMC sampling ---*/
		print("\n[Reduced MCMC]\n  Iteration:");
								  
		for(m_kk=0 ; m_kk<ncyc ; m_kk++)
		{
			print(m_kk);
			/*----------- S A M P L I N G   S T A R T -----------*/
			for(m_k=-nburn ; m_k<nsim ; m_k++)
			{
				/*--- sampling phi ---*/
				if(m_kk < 2)
					dphi = SampPhi(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			
				/*--- sampling sigma and rho ---*/
				if(m_flsigrho && m_kk < 4){ // joint sampling
					[dsig,drho] = SampSigRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
				}
				else{ // separate sampling
					/*--- sampling rho ---*/
					if(m_kk < 3)
						drho = SampRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
	
					/*--- sampling sigma ---*/
					if(m_kk < 4)
						dsig = SampSigma(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
				}
				
				/*--- sampling mu ---*/
				if(m_kk < 5)
					dmu = SampMu(dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			
				/*--- sampling skewness parameter: delta ---*/
				if(m_fls && m_kk < 5 + 2*m_fls) //skew distribution
					ddelta = SampDelta(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
	
				/*--- sampling nu ---*/
				if(m_flt && m_kk < 5 + 2*m_fls + 2*m_flt) // if Student's t
					dnu = SampNu(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
	
				/*--- sampling z0 ---*/
				if(m_fls) // if skewed distribution
					vz0 = SampZ0(dmu, dphi, dsig, drho, ddelta, dnu, vh, vlambda, m_vy);
	
				/*--- sampling h ---*/
				if(m_flmultimove)
					vh = MultiMove(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
				else
					vh = SingleMove(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
	
				/*--- sampling lambda ---*/
				if(m_flt) // if Student's t
					vlambda = SampLambda(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

				if(m_ncount>0 && !imod(m_k, m_ncount)) print(" ",m_k/m_ncount);
			}		 
		    /*------------ S A M P L I N G   E N D --------------*/
			if(m_ncount>0) print("\n  Iteration:");
	
			/*--- fix parameter ---*/
			if(m_kk == 1)
				dphi = m_vposm[0];

			if(m_kk == 3)
			{
				dsig = m_vposm[1];
				drho = m_vposm[2];
			}

			if(m_kk == 4)
				dmu = m_vposm[3];

			if(m_fls && m_kk == 4 + 2*m_fls)
				ddelta = m_vposm[3+m_fls];

			if(m_flt && m_kk == 4 + 2*m_fls + 2*m_flt)
				dnu = m_vposm[3+m_fls+m_flt];
						
		}	/*------------------------------ end of cycle ---*/
	}	//end of ML estimation

	if(m_ncount>0 && !m_flML) println("\n\nTime: ", timespan(time));
}

/*
** SampMu samples parameter mu
*/
SVAZ::SampMu(const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, vy1b, dsigma1sq, dmu1, dmu, dmua;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vy1b = sqrt(dmu_lambda./vlambda) .* (drho*dsig*exp(-vh/2)/dk_delta)
			.* (vy - sqrt(vlambda/dmu_lambda).*exp(vh/2)*ddelta.*(vz0-dc)/sqrt(1-dc^2*ddelta^2));
	dsigma1sq = 1 / (1/m_dsm0^2 + (1 - dphi^2 + (m_ns-1)*(1-dphi)^2/(1-drho^2)) / dsig^2);
	dmu1 = dsigma1sq * (m_dmm0/m_dsm0^2 + ((1-dphi^2)*vh[0] + (1-dphi)/(1-drho^2)*sumc(vh[1:]-dphi*vh[:m_ns-2]-vy1b[:m_ns-2]))/dsig^2);

	if(!m_flML)	//parameter estimation
		dmu = dmu1 + sqrt(dsigma1sq) * rann(1, 1); // mu ~ N(mu1,sigma1^2)
	else	//ML estimation: generate mu for m_kk = 0,1,2,3,4
	{
		if(m_kk < 4)
			dmu = dmu1 + sqrt(dsigma1sq) * rann(1, 1);
		else
		{
			dmu = dmua = m_vposm[3];
			if(m_k>=0)
				m_mmls[m_k][4] = densn((dmua-dmu1)/sqrt(dsigma1sq)) / sqrt(dsigma1sq);
		}		
	}

	return dmu;
}

/*
** SampPhi samples parameter phi by MH algorithm (independent chain) 
*/
SVAZ::SampPhi(const dmu, const dphio, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, vyb, vhb, dsigp, dmup, nite, dphin, dfrac, dphia, dqa;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vyb = sqrt(dmu_lambda./vlambda) .* (drho*dsig*exp(-vh/2)/dk_delta)
			.* (vy - sqrt(vlambda/dmu_lambda).*exp(vh/2)*ddelta.*(vz0-dc)/sqrt(1-dc^2*ddelta^2));
	vhb = vh - dmu;

	/*--- draw candidate ---*/
	dsigp = sqrt(dsig^2*(1-drho^2) / (drho^2*vhb[0]^2 + sumsqrc(vhb[1:m_ns-2])));
	dmup = dsigp^2 * sumc((vhb[1:]-vyb[:m_ns-2]).*vhb[:m_ns-2]) / (dsig^2 * (1-drho^2));
	nite = 0;
	do{
		nite++;
		dphin = dmup + rann(1, 1) * dsigp; // phi ~ N(mu_p, sigma_p^2)
		if(nite>10000){
			print("phi! ");
			dphin = dphio;
		}
	}while(fabs(dphin) >= 1 && nite<=10000);

	/*--- acceptance probability ---*/
	dfrac =	((1+dphin)/(1+dphio))^(m_dap0-0.5) * ((1-dphin)/(1-dphio))^(m_dbp0-0.5);

	if(!m_flML){ //parameter estimation
		if(ranu(1,1) < dfrac){ //accept
			if(m_k>=0) m_cmp++; //acceptance counter
		}
		else dphin = dphio; //reject
	}
	else	//ML estimation: generate phi when m_kk = 0,1
	{
		dphia = m_vposm[0];

		if(m_kk == 0)	//numerator
		{
			if(ranu(1,1) > dfrac)	//reject case
	 			dphin = dphio;
			  
			if(m_k>=0)
			{
				dfrac =	((1+dphia)/(1+dphin))^(m_dap0-0.5) * ((1-dphia)/(1-dphin))^(m_dbp0-0.5); //acceptance probability
				dqa = densn((dphia-dmup)/dsigp)/dsigp / (probn((1-dmup)/dsigp) - probn((-1-dmup)/dsigp)); //candidate density			  
				m_mmls[m_k][0] = min(1, dfrac) * dqa;
			}
		}
		else if(m_k>=0)	//denominator
		{
			dfrac =	((1+dphin)/(1+dphia))^(m_dap0-0.5) * ((1-dphin)/(1-dphia))^(m_dbp0-0.5); //acceptance probability
			m_mmls[m_k][1] = min(1, dfrac);			
			dphin = dphia;
		}
	}
  
	return dphin;
}

/*
** fPosRho computes log posterior distribution of rho
*/
SVAZ::fLogPosRho(const dp, const dfunc, const dsc, const dhe)
{
	decl drho, dmu, dphi, dsig, ddelta, dnu;
	decl dmu_lambda, dc, dk_delta, vyt, vhb;

	// set parameters	
	drho = 2 / (1+exp(-dp)) - 1;
	dmu = m_vpar[0];
	dphi = m_vpar[1];
	dsig = m_vpar[2];
	ddelta = m_vpar[3];
	dnu = m_vpar[4];

	// define variables
	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vyt = m_vy - sqrt(m_vlambda/dmu_lambda).*exp(m_vh/2)*ddelta.*(m_vz0-dc)/sqrt(1-dc^2*ddelta^2);
	vhb = m_vh - dmu;

//	dfunc[0] = log(densbeta((drho+1)/2, m_dar0, m_dbr0)) - 0.5 * (m_ns - 1) * log(1 - drho^2)
//				- 0.5*dmu_lambda/(dk_delta^2*(1-drho^2)) * sumc(exp(-m_vh[:m_ns-2])./m_vlambda[:m_ns-2] .* sqr(vy1t[:m_ns-2]
//					- sqrt(m_vlambda[:m_ns-2]/dmu_lambda).*exp(m_vh[:m_ns-2]/2)*ddelta.*(vhb[1:]-dphi*vhb[:m_ns-2])/dsig));
	dfunc[0] = (m_dar0-m_ns/2-0.5)*log(1+drho) + (m_dbr0-m_ns/2-0.5)*log(1-drho)
				- 0.5*dmu_lambda/(dk_delta^2*(1-drho^2)) * sumc(exp(-m_vh[:m_ns-2])./m_vlambda[:m_ns-2] .* sqr(vyt[:m_ns-2]
					- sqrt(m_vlambda[:m_ns-2]/dmu_lambda)*drho.*exp(m_vh[:m_ns-2]/2)*dk_delta.*(vhb[1:]-dphi*vhb[:m_ns-2])/dsig));
							   
	return 1;
}

/*
** SampRho samples parameter rho by MH algorithm (indepedent chain)
*/
SVAZ::SampRho(const dmu, const dphi, const dsig, const drhoo, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dparo, dparm, dlm, drhom, dsc, dhe, dv, dm, dparn, drhon;
	decl dln, dhn, dlo, dho, dfrac;
	decl drhoa, vpara, dla, dha, dja, dqa;

	m_vpar = dmu | dphi | dsig | ddelta | dnu;
	m_vz0 = vz0;
	m_vh = vh;
	m_vlambda = vlambda;
	m_vy = vy;

	dparo = dparm = log((1+drhoo)/(1-drhoo)); //transformed
	MaxBFGS(fLogPosRho, &dparm, &dlm, 0, 1); //find mode			
	fLogPosRho(dparm, &dlm, 0, 0); //log-posterior density
  	Num1Derivative(fLogPosRho, dparm, &dsc); //score at mode
  	Num2Derivative(fLogPosRho, dparm, &dhe); //hessian at mode

	// draw candidate
	dv = -1 / dhe;
	dm = dparm + dv * dsc;
	dparn = dm + sqrt(dv) * rann(1,1);
	drhon = 2 / (1+exp(-dparn)) - 1;
	
	// compute acceptance probability
	fLogPosRho(dparn, &dln, 0, 0); //log-posterior density at the candidate												    
	dhn = dlm + (dparn-dparm)*dsc - 0.5*(dparn-dparm)^2/dv + log(fabs(-2*drhon/(1-drhon^2))); //log-proposal density at the candidate 
	fLogPosRho(dparo, &dlo, 0, 0); //log-posterior density at the old value
	dho = dlm + (dparo-dparm)*dsc - 0.5*(dparo-dparm)^2/dv + log(fabs(-2*drhoo/(1-drhoo^2))); //log-proposal density at the old value 
	dfrac = exp(dln - dhn - dlo + dho);

	if(!m_flML){ //parameter estimation
		if(ranu(1,1) < dfrac){
			if(m_k>=0) m_cmr++; //acceptance counter
		}
		else drhon = drhoo; //reject
	}
//	else	//ML estimation: generate [sig,rho] when m_kk = 0,1,2,3
//	{
//		if(m_kk < 2)
//		{
//			if(ranu(1,1) > dfrac)	//reject case
//			{
//				dsigetan = dsigetao;
//				drhon = drhoo;
//			}
//		}
//		else
//		{
//			dsigetaa = m_vposm[1];
//			drhoa = m_vposm[2];
//			vpara = log(dsigetaa) | log((1+drhoa)/(1-drhoa));	//transformed
//			fPosSigRho(vpara, &dla, 0, 0);
//	
//			dja = - log(dsigetaa) + log(2) - log(1 - drhoa^2);
//			dha = dlm + (vpara - vparm) ' vsc
//				+ 0.5 * (vpara - vparm) ' mhe * (vpara - vparm)
//				+ dja;
//			dqa = -0.5 * (vpara - vm) ' invert(mS) * (vpara - vm)
//				- log(M_2PI * sqrt(determinant(mS)))
//				+ dja;
//		
//			if(m_kk == 2)	//numerator
//			{
//				if(ranu(1,1) > dfrac)	//reject case
//				{
//		 			dsigetan = dsigetao;
//					drhon = drhoo;
//					dln = dlo;
//					dhn = dho;
//				}
//		
//				if(m_k>=0)
//				{
//					dfrac = exp(dla - dha - dln + dhn);
//					m_mmls[m_k][2] = min(1, dfrac) * exp(dqa);
//				}
//			}
//			else if(m_k>=0)	//denominator
//			{
//				dfrac = exp(dln - dhn - dla + dha);
//				m_mmls[m_k][3] = min(1, dfrac);
//		
//				dsigetan = dsigetaa;
//				drhon = drhoa;
//			}
//		}
//	}
  				  
	return drhon;
}

/*
** SampSigma samples parameter sigma by MH algorithm (indepedent chain)
*/
SVAZ::SampSigma(const dmu, const dphi, const dsigo, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, vhb, dn1, dS1, dsign, vyb, dA, dln, dlo, dfrac;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vhb = vh - dmu;

	/*--- draw candidate ---*/
	dn1 = m_dn0 + m_ns;
	dS1 = m_dS0 + (1-dphi^2)*vhb[0]^2 + sumsqrc(vhb[1:]-dphi*vhb[:m_ns-2])/(1-drho^2);
	dsign = sqrt(1 / rangamma(1, 1, dn1/2, dS1/2)); // sigma^2 ~ InvGamma(n1/2, S1/2)

	/*--- acceptance probability ---*/
	vyb = vy - sqrt(vlambda/dmu_lambda) .* exp(vh/2) * ddelta .* (vz0-dc) / sqrt(1-dc^2*ddelta^2);
	dA = sqrt(dmu_lambda)*drho / (dk_delta*(1-drho^2))
			* sumc(exp(-vh[:m_ns-2]/2)  ./ sqrt(vlambda[:m_ns-2]) .* vyb[:m_ns-2] .* (vhb[1:]-dphi*vhb[:m_ns-2]));
	dln = dA / dsign;
	dlo = dA / dsigo;
	dfrac = exp(dln - dlo);

	if(!m_flML){ //parameter estimation
		if(ranu(1,1) < dfrac){ //accept
			if(m_k>=0) m_cms++; //acceptance counter
		}
		else dsign = dsigo; //reject
	}
  				  
	return dsign;
}

/*
** fPosSigRho computes log posterior distribution
** of sigma and rho 
*/
SVAZ::fLogPosSigRho(const vp, const dfunc, const vsc, const mhe)
{
	decl dsig, drho, dmu, dphi, ddelta, dnu;
	decl dmu_lambda, dc, dk_delta, vyt, vhb;

	dsig = exp(vp[0]);
	drho = 2 / (1+exp(-vp[1])) - 1;
	dmu = m_vpar[0];
	dphi = m_vpar[1];
	ddelta = m_vpar[2];
	dnu = m_vpar[3];

	// define variables
	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vyt = m_vy - sqrt(m_vlambda/dmu_lambda).*exp(m_vh/2)*ddelta.*(m_vz0-dc)/sqrt(1-dc^2*ddelta^2);
	vhb = m_vh - dmu;

	dfunc[0] = (m_dar0-m_ns/2-0.5)*log(1+drho) + (m_dbr0-m_ns/2-0.5)*log(1-drho) - (m_dn0+m_ns+2)*log(dsig) - 0.5*m_dS0/dsig^2
				- 0.5*dmu_lambda/(dk_delta^2*(1-drho^2)) * sumc(exp(-m_vh[:m_ns-2])./m_vlambda[:m_ns-2] .* sqr(vyt[:m_ns-2]
					- sqrt(m_vlambda[:m_ns-2]/dmu_lambda)*drho.*exp(m_vh[:m_ns-2]/2)*dk_delta.*(vhb[1:]-dphi*vhb[:m_ns-2])/dsig))
				- (0.5/dsig^2) * ( sumsqrc(vhb[1:]-dphi*vhb[:m_ns-2]) + (1-dphi^2)*vhb[0]^2 );
							   
	return 1;
}

/*
** SampSigRho samples parameter sigma and rho
** by MH algorithm 
*/
SVAZ::SampSigRho(const dmu, const dphi, const dsigo, const drhoo, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl vparo, vparn, vparm, dlm, vsc, mhe, mS, vm;
	decl dlo, dln, dho, dhn, dsign, drhon, dfrac;
	decl dsiga, drhoa, vpara, dla, dha, dja, dqa, sign;

	m_vpar = dmu | dphi | ddelta | dnu;
	m_vz0 = vz0;
	m_vh = vh;
	m_vlambda = vlambda;
	m_vy = vy;

	vparo = vparm = log(dsigo) | log((1+drhoo)/(1-drhoo)); //transformed

	MaxBFGS(fLogPosSigRho, &vparm, &dlm, 0, 1); //find mode
			
	fLogPosSigRho(vparm, &dlm, 0, 0);
  	Num1Derivative(fLogPosSigRho, vparm, &vsc);
  	Num2Derivative(fLogPosSigRho, vparm, &mhe);
	  
	mS = -invert(mhe);
	vm = vparm + mS * vsc;

	vparn = vm + choleski(mS) * rann(2, 1);	//proposal
				
	fLogPosSigRho(vparo, &dlo, 0, 0);
	fLogPosSigRho(vparn, &dln, 0, 0);
												    
	dsign = exp(vparn[0]);
	drhon = 2 / (1 + exp(-vparn[1])) - 1;
						  
	dho = dlm + (vparo-vparm)'vsc + 0.5*(vparo-vparm)'mhe*(vparo-vparm) - log(dsigo) + log(2) - log(1 - drhoo^2);
	dhn = dlm + (vparn-vparm)'vsc + 0.5*(vparn-vparm)'mhe*(vparn-vparm) - log(dsign) + log(2) - log(1 - drhon^2);

	dfrac = exp(dln - dhn - dlo + dho);

	if(!m_flML)	//parameter estimation
	{
		if(ranu(1,1) < dfrac)
		{
			if(m_k>=0)
				m_cms++;	//acceptance counter
		}
		else	//reject
		{
			dsign = dsigo;
			drhon = drhoo;
		}
	}
	else	//ML estimation: generate [sig,rho] when m_kk = 0,1,2,3
	{
		if(m_kk < 2)
		{
			if(ranu(1,1) > dfrac)	//reject case
			{
				dsign = dsigo;
				drhon = drhoo;
			}
		}
		else
		{
			dsiga = m_vposm[1];
			drhoa = m_vposm[2];
			vpara = log(dsiga) | log((1+drhoa)/(1-drhoa));	//transformed
			fLogPosSigRho(vpara, &dla, 0, 0);
	
			dja = - log(dsiga) + log(2) - log(1 - drhoa^2);
			dha = dlm + (vpara - vparm) ' vsc + 0.5 * (vpara - vparm) ' mhe * (vpara - vparm) + dja;
			dqa = -0.5 * (vpara - vm) ' invert(mS) * (vpara - vm) - 0.5 * log(M_2PI) - 0.5 * logdet(mS,&sign) + dja;
		
			if(m_kk == 2)	//numerator
			{
				if(ranu(1,1) > dfrac)	//reject case
				{
		 			dsign = dsigo;
					drhon = drhoo;
					dln = dlo;
					dhn = dho;
				}
		
				if(m_k>=0)
				{
					dfrac = exp(dla - dha - dln + dhn);
					m_mmls[m_k][2] = min(1, dfrac) * exp(dqa);
				}
			}
			else if(m_k>=0)	//denominator
			{
				dfrac = exp(dln - dhn - dla + dha);
				m_mmls[m_k][3] = min(1, dfrac);
		
				dsign = dsiga;
				drhon = drhoa;
			}
		}
	}
  				  
	return {dsign, drhon};
}

/*
** fLogPosDelta computes log posterior distribution of delta
*/
SVAZ::fLogPosDelta(const dp, const dfunc, const dsc, const dhe)
{
	decl ddelta, dmu, dphi, dsig, drho, dnu;
	decl dmu_lambda, dc, dk_delta, vyt, vhb;

	// set parameters
	ddelta = 2 / (1 + exp(-dp)) - 1;
	dmu = m_vpar[0];
	dphi = m_vpar[1];
	dsig = m_vpar[2];
	drho = m_vpar[3];
	dnu = m_vpar[4];

	// define variables
	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));
	vyt = m_vy - sqrt(m_vlambda/dmu_lambda).*exp(m_vh/2)*ddelta.*(m_vz0-dc)/sqrt(1-dc^2*ddelta^2);
	vhb = m_vh - dmu;

	dfunc[0] = (m_dad0-1)*log(1+ddelta) + (m_dbd0-1)*log(1-ddelta) - m_ns*log(dk_delta)
				- 0.5 * dmu_lambda * exp(-m_vh[m_ns-1]) / (m_vlambda[m_ns-1] * dk_delta^2) * vyt[m_ns-1]^2
				- 0.5 * dmu_lambda/(dk_delta^2*(1-drho^2)) * sumc( exp(-m_vh[:m_ns-2])./m_vlambda[:m_ns-2] .* sqr(vyt[:m_ns-2]
					- sqrt(m_vlambda[:m_ns-2]/dmu_lambda)*drho.*exp(m_vh[:m_ns-2]/2)*dk_delta.*(vhb[1:]-dphi*vhb[:m_ns-2])/dsig) );

	return 1;
}

/*
** SampDelta samples parameter delta
*/
SVAZ::SampDelta(const dmu, const dphi, const dsig, const drho, const ddeltao, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dparo, dparm, dlm, ddeltam, dsc, dhe, dv, dm, dparn, ddeltan;
	decl dln, dhn, dlo, dho, dfrac;
	decl ddeltaa, dpara, dla, dja, dha, dqa;

	m_vpar = dmu | dphi | dsig | drho | dnu;
	m_vz0 = vz0;
	m_vh = vh;
	m_vlambda = vlambda;
	m_vy = vy;

	dparo = dparm = log((1+ddeltao)/(1-ddeltao)); //transformed
	MaxBFGS(fLogPosDelta, &dparm, &dlm, 0, 1); //find mode
	Num1Derivative(fLogPosDelta, dparm, &dsc);
	Num2Derivative(fLogPosDelta, dparm, &dhe);
	dv = -1 / dhe;
	dm = dparm + dv * dsc;
	dparn = dm + sqrt(dv) * rann(1,1);
	ddeltan = 2 / (1 + exp(-dparn)) - 1;

	// compute acceptance probability
	fLogPosDelta(dparn, &dln, 0, 0); //log-posterior density at the candidate
	dhn = dlm + (dparn-dparm)*dsc - 0.5*(dparn-dparm)^2/dv + log(2) - log(1-ddeltan^2); //log-proposal density at the candidate
	fLogPosDelta(dparo, &dlo, 0, 0); //log-posterior density at the old value
	dho = dlm + (dparo-dparm)*dsc - 0.5*(dparo-dparm)^2/dv + log(2) - log(1-ddeltao^2); //log-proposal density at the old value
	dfrac = exp(dln - dhn - dlo + dho);

	if(!m_flML){ //parameter estimation
		if(ranu(1,1) < dfrac){
			if(m_k>=0) m_cmd++; //acceptance counter
		}
		else ddeltan = ddeltao; //reject
	}
	else	//ML estimation: generate delta for m_kk = 0,1,2,3,4 if skew
	{
		if(m_kk < 5)
		{
			if(ranu(1,1) > dfrac)
				ddeltan = ddeltao;	//reject
		}
		else
		{
			ddeltaa = m_vposm[4];
			dpara = log((1+ddeltaa)/(1-ddeltaa));
			fLogPosDelta(dpara, &dla, 0, 0);

			dja = log(2) - log(1-ddeltaa^2);
			dha = dlm + (dpara-dparm)*dsc + 0.5*dhe*(dpara-dparm)^2 + dja;
			dqa = -0.5 * (dpara - dm)^2 / dv - 0.5*log(M_2PI) - 0.5*log(dv) + dja;

			if(m_kk == 5)	//numerator
			{
				if(ranu(1,1) > dfrac)	//reject case
				{
		 			ddeltan = ddeltao;
					dln = dlo;
					dhn = dho;
				}

				if(m_k>=0)
				{
					dfrac = exp(dla - dha - dln + dhn);
					m_mmls[m_k][5] = min(1, dfrac) * exp(dqa);
				}
			}
			else if(m_k>=0)	//denominator
			{
				dfrac = exp(dln - dhn - dla + dha);
				m_mmls[m_k][6] = min(1, dfrac);
		
				ddeltan = ddeltaa;
			}
		}
	}
	
	return ddeltan;
}

/*
** fPosNu computes log posterior distribution of nu
*/
SVAZ::fPosNu(const dp, const dfunc, const dsc, const dhe)
{
	decl dnu, dmu, dphi, dsig, drho, ddelta;
	decl dmu_lambda, dc, dk_delta, vyt, vhb;

	// set parameters
	dnu = exp(dp) + 4;
	dmu = m_vpar[0];
	dphi = m_vpar[1];
	dsig = m_vpar[2];
	drho = m_vpar[3];
	ddelta = m_vpar[4];
	
	// define variables
	dmu_lambda = dnu / (dnu - 2);
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vyt = m_vy - sqrt(m_vlambda/dmu_lambda).*exp(m_vh/2)*ddelta.*(m_vz0-dc)/sqrt(1-dc^2*ddelta^2);
	vhb = m_vh - dmu;

	dfunc[0] = 0.5*m_ns*log(dmu_lambda) + 0.5*dnu*m_ns*log(dnu/2) - m_ns*loggamma(dnu/2)
				- 0.5*dnu*sumc(1 ./ m_vlambda + log(m_vlambda)) + (m_dnn0-1)*log(dnu) - m_dnS0*dnu
				- 0.5 * dmu_lambda * exp(-m_vh[m_ns-1]) / (m_vlambda[m_ns-1] * dk_delta^2) * vyt[m_ns-1]^2
				- 0.5 * dmu_lambda/(dk_delta^2*(1-drho^2)) * sumc(exp(-m_vh[:m_ns-2])./m_vlambda[:m_ns-2] .* sqr(vyt[:m_ns-2]
					- sqrt(m_vlambda[:m_ns-2]/dmu_lambda)*drho.*exp(m_vh[:m_ns-2]/2)*dk_delta.*(vhb[1:]-dphi*vhb[:m_ns-2])/dsig));
				   
	return 1;
}

/*
** SampNu samples parameter nu by A-R M-H algorithm 
*/
SVAZ::SampNu(const dmu, const dphi, const dsig, const drho, const ddelta, const dnuo, const vz0, const vh, const vlambda, const vy)
{
	decl dpm, dpo, dlm, dsc, dhe;
	decl dv, dm, dpn, dnun;
	decl dlo, dln, dho, dhn, dfrac;
	decl dnua, dnuat, dla, dha, dja, dqa;

	m_vpar = dmu | dphi | dsig | drho | ddelta;
	m_vz0 = vz0;
	m_vh = vh;
	m_vlambda = vlambda;
	m_vy = vy;

	dpm = dpo = log(dnuo-4);
	MaxBFGS(fPosNu, &dpm, &dlm, 0, 1); //find mode
	fPosNu(dpm, &dlm, 0, 1);
  	Num1Derivative(fPosNu, dpm, &dsc);
  	Num2Derivative(fPosNu, dpm, &dhe);

	dv = -1 / dhe;
	dm = dpm + dv * dsc;

	dpn = dm + sqrt(dv) * rann(1,1); //generate a candidate
	dnun = exp(dpn) + 4;

	//compute logarithms of posterior and proposal densities at proposed and previous values
	fPosNu(dpn, &dln, 0, 0);												    
	dhn = dlm + (dpn-dpm)*dsc - 0.5*(dpn-dpm)^2/dv - log(dnun-4); //log-proposal density at the candidate 
	fPosNu(dpo, &dlo, 0, 0);
	dho = dlm + (dpo-dpm)*dsc - 0.5*(dpo-dpm)^2/dv - log(dnuo-4); //log-proposal density at the previous value
	
	//compute acceptance probability in M-H part
	if(m_flnuarmh){
		if(dlo < dho) dfrac = 1;
		else if(dln < dhn) dfrac = exp(dho - dlo);
		else dfrac = exp(dln - dhn - dlo + dho);
	}
	else dfrac = exp(dln - dhn - dlo + dho);
						
	if(!m_flML)	//parameter estimation
	{
		if(ranu(1,1) <= dfrac)
		{
			if(m_k>=0) m_cmn++;	//acceptance counter
		}
		else
			dnun = dnuo;	//reject
	}
	else	//ML estimation: generate nu for m_kk = 0,1,2,...,6+2*m_flb
	{
		if(m_kk < 5+2*m_fls)
		{
			if(ranu(1,1) > dfrac)	//reject case
	 			dnun = dnuo;
		}
		else
		{
			dnua = m_vposm[4+m_fls];
			dnuat =	log(dnua);	//transformed
			fPosNu(dnuat, &dla, 0, 0);
				
			dja = - log(dnua-4);
			dha = dlm + (dnuat-dpm)*dsc + 0.5*dhe*(dnuat-dpm)^2 + dja;
			dqa = -0.5 * (dnuat - dm)^2 / dv - 0.5*log(M_2PI) - 0.5*log(dv) + dja - log(1 - probn((log(4) - dm) / sqrt(dv)));

			if(m_kk == 5+2*m_fls)	//numerator
			{
				if(ranu(1,1) > dfrac)	//reject case
				{
		 			dnun = dnuo;
					dln = dlo;
					dhn = dho;
				}

				if(m_k>=0)
				{
					dfrac = exp(dla - dha - dln + dhn);
					m_mmls[m_k][5+2*m_fls] = min(1, dfrac) * exp(dqa);
				}
			}
			else if(m_k>=0)	//denominator
			{
				dfrac = exp(dln - dhn - dla + dha);
				m_mmls[m_k][6+2*m_fls] = min(1, dfrac);
		
				dnun = dnua;
			}
		}
	}

	return dnun;
}

/*
** SampZ0 samples auxiliary variable z0 
*/
SVAZ::SampZ0(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vh, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, vhb, dsigzt, dmuzt, vz0, t;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vhb = vh - dmu;

	dsigzt = sqrt( 1 / ( 1 + ddelta^2 / ((1-ddelta^2)*(1-drho^2)) ) );
	vz0 = zeros(m_ns, 1) - 1;
	for(t=0; t<m_ns-1; t++){
		dmuzt = dsigzt^2 * sqrt(dmu_lambda/vlambda[t]) * exp(-vh[t]/2) * ddelta
					/ (dk_delta^2 * (1-drho^2) * sqrt(1-dc^2*ddelta^2))
					* (vy[t] + sqrt(vlambda[t]/dmu_lambda)*exp(vh[t]/2)*dc*ddelta/sqrt(1-dc^2*ddelta^2)
						- sqrt(vlambda[t]/dmu_lambda)*drho*exp(vh[t]/2)*dk_delta*(vhb[t+1] - dphi*vhb[t]) / dsig);
		while(vz0[t] <= 0){
			vz0[t] = dmuzt + dsigzt * rann(1,1);
		}
	}
	t = m_ns - 1;
	dsigzt = sqrt(1 - ddelta^2);
	dmuzt = dsigzt^2 * sqrt(dmu_lambda/vlambda[t]) * exp(-vh[t]/2) * ddelta
					/ (dk_delta^2 * sqrt(1-dc^2*ddelta^2))
					* (vy[t] + sqrt(vlambda[t]/dmu_lambda)*exp(vh[t]/2)*dc*ddelta/sqrt(1-dc^2*ddelta^2));
	while(vz0[t] <= 0){
		vz0[t] = dmuzt + dsigzt * rann(1,1);
	}
	
	return vz0;		
}

/*
** SingleMove samples log-volatilities h_1, h_2, ..., h_T separately
*/
SVAZ::SingleMove(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vho, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, vxb, vlam1, vlam2, vlamz0, vhn, vmuh, vln, vlo, vac, dsigh1, dsigh2T1, dsighT, t, T;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));
	vlam1 = sqrt(dmu_lambda./vlambda) * drho / (dk_delta * (1-drho^2) * dsig);
	vlam2 = dmu_lambda ./ (vlambda * dk_delta^2 * (1-drho^2));	
	vlamz0 = sqrt(vlambda/dmu_lambda) * ddelta .* (vz0 - dc) / sqrt(1 - dc^2*ddelta^2);

	/*--- Single-move ---*/
	vhn = vho;
	vmuh = vln = vlo = vac = zeros(m_ns, 1);

	// t = 1
	dsigh1 = sqrt(dsig^2 / (1 - dphi^2 + dphi^2/(1-drho^2)));
	vmuh[0]	= dsigh1^2 * (-0.5 + dphi*(vho[1]-(1-dphi)*dmu) / (dsig^2*(1-drho^2)) + (1-dphi^2)*dmu/dsig^2);
	vhn[0] = vmuh[0] + dsigh1 * rann(1,1);
	vln[0] = -0.5 * vlam2[0] * exp(-vhn[0]) * (vy[0] - vlamz0[0] * exp(vhn[0]/2))^2
				+ vlam1[0] * exp(-vhn[0]/2) * (vy[0] - vlamz0[0] * exp(vhn[0]/2)) * (vhn[1]-dmu - dphi*(vhn[0]-dmu));
	vlo[0] = -0.5 * vlam2[0] * exp(-vho[0]) * (vy[0] - vlamz0[0] * exp(vho[0]/2))^2
				+ vlam1[0] * exp(-vho[0]/2) * (vy[0] - vlamz0[0] * exp(vho[0]/2)) * (vhn[1]-dmu - dphi*(vho[0]-dmu));
	vac[0] = ranu(1,1) < exp(vln[0] - vlo[0]);
	vhn[0] = vac[0] ? vhn[0] : vho[0];

	// t = 2,...,T-1
	dsigh2T1 = sqrt(dsig^2 * (1-drho^2) / (1+dphi^2));
	for(t=1; t<=m_ns-2; t++){
		vmuh[t] = dsigh2T1^2 * ( -0.5 + (dphi*(vhn[t+1]+vhn[t-1])+(1-dphi)^2*dmu)/(dsig^2*(1-drho^2))
									+ vlam1[t-1] * exp(-vhn[t-1]/2) * (vy[t-1] - vlamz0[t-1]*exp(vhn[t-1]/2)) );
		vhn[t] = vmuh[t] + dsigh2T1 * rann(1,1);
		vln[t] = -0.5 * vlam2[t] * exp(-vhn[t]) * (vy[t] - vlamz0[t] * exp(vhn[t]/2))^2
					+ vlam1[t] * exp(-vhn[t]/2) * (vy[t] - vlamz0[t] * exp(vhn[t]/2)) * (vhn[t+1]-dmu - dphi*(vhn[t]-dmu));
		vlo[t] = -0.5 * vlam2[t] * exp(-vho[t]) * (vy[t] - vlamz0[t] * exp(vho[t]/2))^2
					+ vlam1[t] * exp(-vho[t]/2) * (vy[t] - vlamz0[t] * exp(vho[t]/2)) * (vhn[t+1]-dmu - dphi*(vho[t]-dmu));
		vac[t] = ranu(1,1) < exp(vln[t] - vlo[t]);
		vhn[t] = vac[t] ? vhn[t] : vho[t];
	}

	// t = T
	dsighT = sqrt(dsig^2 * (1-drho^2));
	T = m_ns - 1;
	vmuh[T] = dsighT^2 * ( -0.5 + ((1-dphi)*dmu + dphi*vhn[T-1])/(dsig^2*(1-drho^2))
								+ vlam1[T-1] * exp(-vhn[T-1]/2) * (vy[T-1] - vlamz0[T-1]*exp(vhn[T-1]/2)) );
	vhn[T] = vmuh[T] + dsighT * rann(1,1);
	vln[T] = -0.5 * dmu_lambda * exp(-vhn[T]) / (vlambda[T]*dk_delta^2) * (vy[T] - vlamz0[T]*exp(vhn[T]/2))^2;
	vlo[T] = -0.5 * dmu_lambda * exp(-vho[T]) / (vlambda[T]*dk_delta^2) * (vy[T] - vlamz0[T]*exp(vho[T]/2))^2;
	vac[T] = ranu(1,1) < exp(vln[T] - vlo[T]);
	vhn[T] = vac[T] ? vhn[T] : vho[T];
	
	if(!m_flML && m_k>=0) m_cmh += sumc(vac);

	return vhn;
}

/*
** MultiMove samples log-volatilities h_1, h_2, ..., h_T for each block
*/
SVAZ::MultiMove(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambda, const vy)
{
	decl dmu_lambda, dc, dk_delta, dgam, vz0b;
	decl nite, vk, ir, id, ird, fli, vh_s, vetao, vetaoi;
	decl vyi, vho, vhn, dhrd1, vhh, i, j, t;
	decl vmu, vsig2, dlfhrd1, vder, vder1, vdel, vA, vB;
	decl vlambdai, vz0bi, vDh, vKh, vJh, vbh, vgn;
	decl vgh, vyh, dLh, mBh, mQh;
	decl vD, vZ, mG, vH, vH0, da, da0, dP, dH2, dV, du0;
	decl ve, vDinv, vK, vL, mJ, vu, dr, dU;
	decl deps, dC, dCinv;
	decl vmun, vsign2, dlfhrdn1, vmuo, vsigo2, dlfhrdo1;
	decl vetan, dposn, dposo, dpron, dproo, dfrac;
	decl car;

	// Define parameters and variable
	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));
	dgam = exp(0.5*dmu) / sqrt((1-dc^2*ddelta^2)*dmu_lambda);
	vz0b = vz0 - dc;

	// Setup
	nite = 5; //# of iteration for h_hat
	vh_s = vh - dmu; //sample box
	vk = 0 ~ floor(m_ns*(range(1,m_nKm)+ranu(1,m_nKm))/(m_nKm+2)) ~ m_ns-1; //stochastic knots
	vetao = ((vh_s[1:] - dphi*vh_s[:m_ns-2]) / dsig) | 0; //eta_old
	
	for(i=0 ; i<m_nKm+1 ; i++)
	{
		ir  = vk[i];
		id  = vk[i+1] - vk[i]; 
		ird = vk[i+1];

		if(id<2){if(!m_flML && m_k>=0) m_cah++;}
		else
		{
			fli = i < m_nKm;	// fli = 1 if it is not the last block		 
		
			vyi = m_vy[ir:ird]; //r:r+d, ((id+1), 2) matrix
			vlambdai = vlambda[ir:ird];
			vz0bi = vz0b[ir:ird];

			if(i)
			{
				vetaoi = vetao[ir:ird];
				vho = vh_s[ir] | zeros(id, 1); //current(old) point, r:r+d
				vhn	= vho; //new point, r:r+d
		
				for(t=0 ; t<id ; t++)
					vho[t+1] = dphi * vho[t] + vetaoi[t]*dsig;  //update r+1:r+d	

				if(fli) dhrd1 = dphi * vho[id] + vetaoi[id]*dsig;
			}
			else	//initial case
			{
				ir = -1;
				id++;
				vyi = 0 | vyi;	//add y_0 = 0 = unconditional means
				vlambdai = dmu_lambda | vlambdai;
				vz0bi = 0 | vz0bi;
				vho = 0 | vh_s[:ird];	//add h_0
				vhn = zeros(id+1, 1);
				dhrd1 = vh_s[ird+1];
				vetaoi = vh_s[0]*sqrt(1-dphi^2)/dsig | vetao[:ird];
			}

			/* finding mode	& draw candidate */
			for(j=0 ; j<=nite ; j++)
			{
				vhh = j ? vhn : vho; //h_hat, r:r+d
			   
				// vmu, vsig2, vder, vder1 are t=(r:r+d),(id+1)*1 vector
				vmu = vyi[:id-1] - ddelta*vz0bi[:id-1] .* sqrt(vlambdai[:id-1]) .* exp(vhh[:id-1]/2) * dgam	// y_{1,t} - mu_t (t<n)
					- drho * sqrt((1-ddelta^2)*vlambdai[:id-1]) .* exp(vhh[:id-1]/2) * dgam .* (vhh[1:] - dphi*vhh[:id-1])/dsig;
				vsig2 = (1-drho^2) * (1-ddelta^2) * vlambdai .* exp(vhh) * dgam^2;	// sigma_t^2 (t<n)
			
				vder = -0.5 * ddelta * vz0bi[:id-1] .* sqrt(vlambdai[:id-1]) .* exp(vhh[:id-1]/2) * dgam	//dmu_t/dh_t (t<n)
					 - drho * sqrt((1-ddelta^2)*vlambdai[:id-1]) .* exp(vhh[:id-1]/2) * dgam .* ((vhh[1:]-dphi*vhh[:id-1])/2-dphi)/dsig;

				vder1 = - drho * sqrt((1-ddelta^2)*vlambdai[:id-1]) .* exp(vhh[:id-1]/2) * dgam / dsig;
				vder1 = i == 0 ? (0|0|vder1[1:]) : (0|vder1);	//dmu_{t-1}/dh_t

				if(fli)	// if not the last block (s+m < n)
				{
					vmu |= vyi[id] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhh[id]/2) * dgam
						 - drho * sqrt((1-ddelta^2)*vlambdai[id]) * exp(vhh[id]/2) * dgam * (dhrd1 - dphi*vhh[id])/dsig;
					vder |= -0.5 * ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhh[id]/2) * dgam
						  - drho * sqrt((1-ddelta^2)*vlambdai[id]) * exp(vhh[id]/2) * dgam * ((dhrd1-dphi*vhh[id])/2-dphi) / dsig;
					dlfhrd1 = -0.5 * (dhrd1 - dphi * vhh[id])^2 / dsig^2;	// additional term at t = s+m < n
				}
				else{
					vmu |= vyi[id] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhh[id]/2) * dgam;
					vsig2[id] = (1-ddelta^2) * vlambdai[id] * exp(vhh[id]) * dgam^2;
					vder |= -0.5 * ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhh[id]/2) * dgam;
					dlfhrd1 = 0;
				}
				if(!i) vmu[0] = 0;
	
				//vdel, vA, vB are t=(r+1:r+d), id*1 vector			
				vdel = -0.5	+ 0.5 * vmu[1:].^2 ./ vsig2[1:] - vmu[1:] ./ vsig2[1:] .* vder[1:] - vmu[:id-1] ./ vsig2[:id-1] .* vder1[1:];
				vA = 0.5 + vder[1:].^2 ./ vsig2[1:] + vder1[1:].^2 ./ vsig2[:id-1];
				vB = 0 | (vder[1:id-1] .* vder1[2:] ./ vsig2[1:id-1]);
			
				if(fli)	// additional term for d_t at t = s+m < n
				{
					vdel[id-1] += dphi * (dhrd1 - dphi * vhh[id]) / dsig^2;
					vA[id-1] += dphi^2 / dsig^2;
				}
				  
 
				vDh = vJh = vKh = vbh = vZ = zeros(id, 1); //r+1:r+d
				
				vDh[0] = vA[0];
				vKh[0] = sqrt(vDh[0]);
				vbh[0] = vdel[0];
				for(t=1 ; t<id ; t++)
				{
					vDh[t] = vA[t] - vB[t]^2 / vDh[t-1];
					vKh[t] = sqrt(vDh[t]);
					vJh[t] = vB[t] / vKh[t-1];
					vbh[t] = vdel[t] - vJh[t] * vbh[t-1] / vKh[t-1];
				}
				vJh |= 0;	//r+1:r+d+1
			
				vgh = vhh[1:] + ((vJh[1:id-1].*vhh[2:]./vKh[:id-2])|0);
				vyh = vgh + vbh ./ vDh;		//r+1:r+d
			
				//simulation smoother
			
				vZ = 1 + dphi * vJh[1:] ./ vKh;
				mG = (1 ~ vJh[1:]*dsig) ./ vKh;
				vH = 0 ~ dsig;
			
				vH0 = i ? (0 ~ dsig) : (0 ~ dsig/sqrt(1-dphi^2));
				da = da0 = i ? dphi * vh_s[ir] : 0;
				dP = vH0 * vH0';
			
				ve = vDinv = vK = vL = vu = zeros(id, 1);
				mJ = zeros(id, 2);
							
				for(t=0 ; t<id ; t++)		//Kalman filter
				{
					ve[t] = vyh[t] - vZ[t] * da;
					vDinv[t] = 1 / (vZ[t]^2*dP + mG[t][]*mG[t][]');
					vK[t] = (dphi*dP*vZ[t] + vH*mG[t][]') * vDinv[t];
					vL[t] = dphi - vK[t] * vZ[t];
					mJ[t][] = vH - vK[t] * mG[t][];
			
					da = dphi * da + vK[t] * ve[t];
					dP = dphi * dP * vL[t] + vH * mJ[t][]';
				}

				if(j < nite)	//finding mode
				{
			
					dr = dU = 0;
				
					for(t=id-1 ; t>=0 ; t--)	//Simulation smoother
					{
						dC = vH * (unit(2) - mG[t][] ' vDinv[t] * mG[t][] - mJ[t][] ' dU * mJ[t][]) * vH';
						deps = 0;
						vu[t] = vH * (mG[t][] ' vDinv[t] * ve[t] + mJ[t][] ' dr) + deps;
						dV = vH * (mG[t][] ' vDinv[t] * vZ[t] + mJ[t][] ' dU * vL[t]);
						dCinv = 1 / dC;
						dr = vZ[t] * vDinv[t] * ve[t] + vL[t] * dr - dV * dCinv * deps;
						dU = vZ[t]^2 * vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv;
					}

					du0 = vH0 * vH0 ' dr;
				
					vhn[1] = da0 + du0;			//t=r+1
					for(t=0 ; t<id-1 ; t++)
						vhn[t+2] = dphi * vhn[t+1] + vu[t];
			 	}	//end finding mode
		   
				else	//draw candidate
				{
					car = 0;
					do
					{
						car++;
						dr = dU = 0;
														 
						for(t=id-1 ; t>=0 ; t--)	//Simulation smoother
						{
							dC = vH * (unit(2) - mG[t][] ' vDinv[t] * mG[t][] - mJ[t][] ' dU * mJ[t][]) * vH';
							deps = sqrt(dC) * rann(1, 1); 
							vu[t] = vH * (mG[t][] ' vDinv[t] * ve[t] + mJ[t][] ' dr) + deps;
							dV = vH * (mG[t][] ' vDinv[t] * vZ[t] + mJ[t][] ' dU * vL[t]);
							dCinv = 1 / dC;
							dr = vZ[t] * vDinv[t] * ve[t] + vL[t] * dr - dV * dCinv * deps;
							dU = vZ[t]^2 * vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv;
						}
					
						dH2 = vH0 * vH0';
						dC = dH2 - dH2^2 * dU;
						du0 = dH2 * dr + sqrt(dC) * rann(1, 1);
					
						vhn[1] = da0 + du0;			//t=r+1
						for(t=0 ; t<id-1 ; t++)
							vhn[t+2] = dphi * vhn[t+1] + vu[t];	//propose new h
									 
						//AR step
						dLh = sumc(-0.5 * vhh - 0.5 * vmu.^2 ./ vsig2) + dlfhrd1;
					
						mBh = (0 ~ diag(vB[1:])) | 0;
						mQh = diag(vA) + mBh + mBh';
									
						dpron = dLh + vdel ' (vhn - vhh)[1:] - 0.5 * (vhn - vhh)[1:] ' mQh * (vhn - vhh)[1:];
					
						vmun = vyi[:id-1] - ddelta * vz0bi[:id-1] .* sqrt(vlambdai[:id-1]) .* exp(vhn[:id-1]/2) * dgam
							 - drho * sqrt((1-ddelta^2)*vlambdai[:id-1]) .* exp(vhn[:id-1]/2) * dgam .* (vhn[1:]-dphi*vhn[:id-1])/dsig;
						vsign2 = (1 - drho^2) * (1 - ddelta^2) * vlambdai .* exp(vhn) * dgam^2;
				
						if(fli)
						{
							vmun |= vyi[id] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhn[id]/2) * dgam
								  - drho * sqrt((1 - ddelta^2) * vlambdai[id]) * exp(vhn[id]/2) * dgam * (dhrd1 - dphi*vhn[id]) / dsig;
							dlfhrdn1 = -0.5 * (dhrd1 - dphi * vhn[id])^2 / dsig^2;
						}
						else
						{
							vmun |= vyi[id] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vhn[id]/2) * dgam;
							vsign2[id] = (1 - ddelta^2) * vlambdai[id] * exp(vhn[id]) * dgam^2;
							dlfhrdn1 = 0;
						}
						if(!i) vmun[0] = 0;
				
						dposn = sumc(-0.5*vhn - 0.5 * vmun.^2 ./ vsign2) + dlfhrdn1;
				
						if(!m_flML && m_k>=0) m_cah++;

						if(car>=100)
							print("h! ");

					}while(ranu(1, 1) >= exp(dposn - dpron) && car<100);

				}	//end draw candidate
  						   
			}	//end finding mode & draw candidate	

			//MH step
			dproo = dLh + vdel ' (vho - vhh)[1:] - 0.5 * (vho - vhh)[1:] ' mQh * (vho - vhh)[1:];
									  
			vmuo = vyi[:id-1] - ddelta * vz0bi[:id-1] .* sqrt(vlambdai[:id-1]) .* exp(vho[:id-1]/2) * dgam
				 - drho * sqrt((1-ddelta^2) * vlambdai[:id-1]) .* exp(vho[:id-1]/2) * dgam .* (vho[1:] - dphi * vho[:id-1]) / dsig;
			vsigo2 = (1 - drho^2) * (1 - ddelta^2) * vlambdai .* exp(vho) * dgam^2;
			if(fli)
			{
				vmuo |= vyi[id] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vho[id]/2) * dgam
					  - drho * sqrt((1 - ddelta^2) * vlambdai[id]) * exp(vho[id]/2) * dgam * (dhrd1 - dphi * vho[id]) / dsig;
				dlfhrdo1 = -0.5 * (dhrd1 - dphi * vho[id])^2 / dsig^2;
			}
			else{
				vmuo |= vyi[id][0] - ddelta * vz0bi[id] * sqrt(vlambdai[id]) * exp(vho[id]/2) * dgam;
				vsigo2[id] = (1 - ddelta^2) * vlambdai[id] * exp(vho[id]) * dgam^2;
				dlfhrdo1 = 0;
			}
			if(!i) vmuo[0] = 0;
		 
			dposo = sumc(-0.5 * vho	- 0.5 * vmuo.^2 ./ vsigo2) + dlfhrdo1;

			dfrac = exp(  dposn + min(dposo|dproo) - dposo - min(dposn|dpron) );
							
			if(ranu(1, 1) < dfrac){
				if(!m_flML && m_k>=0) m_cmh++;
			}					   
			else
				vhn[1:] = vho[1:];
	
			vh_s[ir+1:ird] = vhn[1:];
		
		}
	}	//end sampling all blocks		

	return vh_s + dmu;
}

/*
** SampLambda samples parameter lambda by MH algorithm 
*/
SVAZ::SampLambda(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu, const vz0, const vh, const vlambdao, const vy)
{
	decl dmu_lambda, dc, dk_delta, vhb, da, vb, vlambdan, t, vln, vlo, vac;

	dmu_lambda = dnu / (dnu - 2);
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));	
	vhb = vh - dmu;

	/*--- draw candidate ---*/
	da = dnu + 1;
	vb = dnu + dmu_lambda * exp(-vh[:m_ns-2]) .* sqr(vy[:m_ns-2]) / (dk_delta^2 * (1-drho^2)); // t=1,...,T-1
	vb |= dnu + dmu_lambda * exp(-vh[m_ns-1]) * sqr(vy[m_ns-1]) / dk_delta^2; // t=T
	vlambdan = zeros(m_ns, 1);
	for(t=0; t<m_ns; t++){
		vlambdan[t] = 1 / rangamma(1, 1, da/2, vb[t]/2); // lambda_t ~ InvGamma(a/2,b/2)
	}
	
	/*--- acceptance probability ---*/
	vln = sqrt(dmu_lambda./vlambdan[:m_ns-2]) .* exp(-vh[:m_ns-2]/2) .* vy[:m_ns-2] / (dk_delta^2 * (1-drho^2)) // t=1,...,T-1
			.* ( ddelta .* (vz0[:m_ns-2] - dc) / sqrt(1 - dc^2 * ddelta^2) + drho * dk_delta * (vhb[1:] - dphi * vhb[:m_ns-2])/dsig );
	vln |= sqrt(dmu_lambda/vlambdan[m_ns-1]) * exp(-vh[m_ns-1]/2) * vy[m_ns-1] / dk_delta^2 // t=T
			* ddelta * (vz0[m_ns-1] - dc) / sqrt(1 - dc^2 * ddelta^2);
	vlo = sqrt(dmu_lambda./vlambdao[:m_ns-2]) .* exp(-vh[:m_ns-2]/2) .* vy[:m_ns-2] / (dk_delta^2 * (1-drho^2)) // t=1,...,T-1
			.* ( ddelta .* (vz0[:m_ns-2] - dc) / sqrt(1 - dc^2 * ddelta^2) + drho * dk_delta * (vhb[1:] - dphi * vhb[:m_ns-2])/dsig );
	vlo |= sqrt(dmu_lambda/vlambdao[m_ns-1]) * exp(-vh[m_ns-1]/2) * vy[m_ns-1] / dk_delta^2 // t=T
			* ddelta * (vz0[m_ns-1] - dc) / sqrt(1 - dc^2 * ddelta^2);
	vac = ranu(m_ns, 1) .< exp(vln - vlo);
	vlambdan = vac .? vlambdan .: vlambdao;

	if(!m_flML && m_k>=0) m_cml += sumc(vac);

	return vlambdan;
}

/*
** fLogPri evaluates log-prior density
*/
SVAZ::fLogPri(const dmu, const dphi, const dsig, const drho, const ddelta, const dnu)
{
	decl dlogpri;

	dlogpri
	= -0.5*log(2*M_PI) - log(m_dsm0) - 0.5*(dmu-m_dmm0)^2/m_dsm0^2
	+ log(densbeta((dphi+1)/2, m_dap0, m_dbp0)) - log(2)
	+ log(densbeta((drho+1)/2, m_dar0, m_dbr0)) - log(2)
	+ 0.5*m_dn0*log(m_dS0/2) - loggamma(m_dn0/2) - (m_dn0+2)*log(dsig) - 0.5*m_dS0/dsig^2;
	
	if(m_flt) dlogpri += log(densgamma(dnu, m_dnn0, m_dnS0));	
	if(m_fls) dlogpri += log(densbeta((ddelta+1)/2, m_dad0, m_dbd0)) - log(2);
	
	return dlogpri;
}

/*
** ML computes marginal likelihood 
*/
SVAZ::ML(const nI, const nR, const nsim)
{
	decl msamp, vnusamp, vnu, vf, dnum;
	decl vlik, dpos, dposse, iBm, i;
	decl dvm, vfl, mlsh, mOm, mOmsum, vtaua;
	decl dphi, dsig, drho, dmu, ddelta, dnu, dpri;
	decl mml;
	decl times = timer();

	m_vposm = loadmat(sprint(m_sdirectory, "/posstat_",m_sfilename,".csv"))[][0];
		//posterior mean of phi, sigma, rho, mu, delta (if skew), nu (if t), xi, sigma_u

	iBm = min(500, nsim/2);			//bandwidth

	/*--- likelihood part ---*/
	vlik = PFilter(nI, nR);

	/*--- posterior part ---*/

	m_flML = 1;
	MCMC(nsim, 0.1*nsim);		//reduced run
	m_flML = 0;

	dvm = meanc(m_mmls);
	vfl = <1; -1; 1; -1; 1>;	// phi, (sigma,rho), mu
	if(m_fls)
		vfl |= <1; -1>;	// delta
	if(m_flt)
		vfl |= <1; -1>;	// nu
	dpos = log(dvm) * vfl;
						 
	mlsh = m_mmls - dvm;
	mOmsum = zeros(5+2*m_fls+2*m_flt, 5+2*m_fls+2*m_flt);
	for(i=0 ; i<iBm ; i++){
		mOm = mlsh[i:][] ' mlsh[:nsim-i-1][] / nsim;
		mOmsum += (1 - i/(iBm+1)) * (mOm + mOm' * (i>0));
	}

	mOmsum /= nsim;
	vtaua = vfl ./ dvm';
	dposse = sqrt(vtaua ' mOmsum * vtaua);


	/*--- prior part ---*/
	dphi = m_vposm[0];
	dsig = m_vposm[1];
	drho = m_vposm[2];
	dmu = m_vposm[3];
	ddelta = m_fls ? m_vposm[4] : 0;
	dnu = m_flt ? m_vposm[4+m_fls] : M_INF_POS;
	dpri = fLogPri(dmu, dphi, dsig, drho, ddelta, dnu);


	/*--- output ---*/
	print("\n\n========================================"
		  "\n     Marginal likelihood estimation"
		  "\n========================================"
		  "\n                                SVAZ.ox"
		  "\nModel: SV-");
	if(m_fls)
		print("S");
	if(m_flt)
		print("T");
	else
		print("N");
		 
	print("\n\n    ESTIMATION RESULT    (S.E.)"
		  "\n--------------------------------");
	print("%r",{" Likelihood", " Prior", " Posterior",
			    " log-ML"},
		  "%cf", {"%10.3f", "%8.3f"},
		    (meanc(vlik) ~ sqrt(varc(vlik)/nR))
		  | (dpri ~ 0)
		  | (dpos ~ dposse)
		  | ((meanc(vlik)+dpri-dpos)
		  	  ~ sqrt(varc(vlik)/nR+dposse^2)));
	print("--------------------------------"
		  "\n[Particle filter]\nParticle", "%7.0f", nI,
		  "\nRepricants", "%5.0f", nR,
		  "\n\n[Reduced MCMC]\nSample", "%9.0f", nsim);

	mml = (meanc(vlik) ~ sqrt(varc(vlik)/nR))
		  | (dpri ~ 0)
		  | (dpos ~ dposse)
		  | ((meanc(vlik)+dpri-dpos) ~ sqrt(varc(vlik)/nR+dposse^2));

	savemat(sprint(m_sdirectory, "/ml_",m_sfilename,".csv"), mml);

	println("\n\nTime: ", timespan(times));
}

/*
** PFilter estimates likelihood by particle filter
*/
SVAZ::PFilter(const nI, const nR)
{
	decl dphi, dsig, drho, dmu, ddelta, dnu;
	decl dmu_lambda, dc, dk_delta, vlambda, vz0;
	decl vf, vlik, vh, vmuy, vsigy, vw, vfhat, vmu, vgj, vuc, vuci;
	decl vlambdai, vz0i, vind, cuci, cuc, vhi, vfhati, vgji, i, j, t;

	dphi = m_vposm[0];
	dsig = m_vposm[1];
	drho = m_vposm[2];
	dmu = m_vposm[3];
	ddelta = m_fls ? m_vposm[4] : 0;
	dnu = m_flt ? m_vposm[4+m_fls] : M_INF_POS;

	if(m_flt) dmu_lambda = dnu / (dnu - 2);
	else dmu_lambda = 1;
	dc = sqrt(2 / M_PI);
	dk_delta = sqrt((1-ddelta^2) / (1-dc^2*ddelta^2));

	vlambda = ones(nI, 1);
	vz0 = ones(nI, 1) * dc;

	vf = zeros(m_ns, 1);
	vlik = zeros(nR, 1);
	
	if(nR > 1)
		print("\n[Particle filter] Iteration:");

	/*---------- Iteration START -----------*/
	for(j=0 ; j<nR ; j++)
	{
		if(nR > 1)
			print(j);

		/* Procedure 1 */
		// generate h_1^i, lambda_1^i, z_01^i from the unconditional distributions 
		vh = dmu + rann(nI, 1) * dsig / sqrt(1 - dphi^2); // generate h_1
		if(m_flt)
			vlambda = 1 ./ rangamma(nI, 1, dnu/2, dnu/2);
		if(m_fls){
			for(i=0; i<nI; i++){
				vz0[i] = rann(1,1);
				while(vz0[i] <= 0){
					vz0[i] = rann(1,1);
				}
			}
		}

		// (a) compute w_i = f(x_1, y_1 | h_1^i, lambda_1^i, z_01^i, theta) and record their mean 
		vmuy = ddelta*(vz0-dc)/sqrt(1-dc^2*ddelta^2) .* sqrt(vlambda/dmu_lambda) .* exp(vh/2);
		vsigy = dk_delta * sqrt(vlambda/dmu_lambda) .* exp(vh/2);
		vw = sqrt(0.5/M_PI) .* exp(-0.5 * sqr((m_vy[0] - vmuy) ./ vsigy)) ./ vsigy;
		vf[0] = meanc(vw);

		// (b) compute fhat(h_1^i, lambda_1^i, z_01^i | x_1, y_1, theta)
		vfhat = vw / sumc(vw);
		
		/*---------- Filtering START -----------*/
		for(t=0 ; t<m_ns-1 ; t++)
		{
			/* Procedure 2 */
			vmu = dmu + dphi * (vh - dmu)
				  + drho*dsig*(sqrt(dmu_lambda./vlambda)*m_vy[t]./(dk_delta.*exp(vh/2)) - ddelta*(vz0-dc)/sqrt(1-ddelta^2));
			vgj = sqrt(0.5/M_PI) .* exp(-vmu/2 - 0.5*sqr(m_vy[t+1])./(dk_delta^2 .* exp(vmu))) ./ dk_delta .*vfhat;
			vgj /= sumc(vgj);	// 1st stage weight generated from g(k|Y_{t+1}) =~ pi^k * f(y_{t+1}|mu_{t+1}^k)
			vuc = ranmultinomial(nI, vgj);	// generate the number of successes for each of nI particles
			vuci = vecindex(vuc);	// return index with at least one success
		
			vind = <>;	// vector of index in vuci: each index appears at the number of occurance in vuc
			for(i=0 ; i<rows(vuci) ; i++)
			{
				cuci = vuci[i];					// example: nI = 3, vuci = (0, 2)
				cuc = vuc[cuci];				// vuc = (1, 0, 3, 0)
				vind |= cuci * ones(cuc, 1);	// vind = (0, 2, 2, 2)
			}
		
			vhi = vh[vind];			// extract particles for the corresponding index
			vfhati = vfhat[vind];
			vgji = vgj[vind];
		
			vlambdai = vlambda[vind];
			vz0i = vz0[vind];

			vmu = dmu + dphi * (vhi - dmu)
				  + drho*dsig*(sqrt(dmu_lambda./vlambdai)*m_vy[t]./(dk_delta.*exp(vhi/2)) - ddelta*(vz0i-dc)/sqrt(1-ddelta^2));
			vh = vmu + rann(nI, 1) * dsig * sqrt(1-drho^2);
			
			if(m_flt)
				vlambda = 1 ./ rangamma(nI, 1, dnu/2, dnu/2);
			if(m_fls){
				for(i=0; i<nI; i++){
					vz0[i] = rann(1,1);
					while(vz0[i] <= 0){
						vz0[i] = rann(1,1);
					}
				}
			}
				
			// (a) compute w_i = f(x_t+1,y_t+1|*) fhat(h_t^i,lambda_t^i,z_0t^i|*) / g(h_t^i,lambda_t^i,z_0t^i |*) and record their mean 
			vmuy = ddelta*(vz0-dc)/sqrt(1-dc^2*ddelta^2) .* sqrt(vlambda/dmu_lambda) .* exp(vh/2);
			vsigy = dk_delta * sqrt(vlambda/dmu_lambda) .* exp(vh/2);
			vw = sqrt(0.5/M_PI) .* exp(-0.5 * sqr((m_vy[t+1] - vmuy) ./ vsigy)) ./ vsigy .* vfhati ./ vgji;
			vf[t+1] = meanc(vw);

			// (b) compute fhat(h_t+1^i, lambda_t+1^i, z_0,t+1^i | x_t, y_t, theta)
			vfhat = vw / sumc(vw);
			
		}
		/*---------- Filtering END -----------*/
	  		   
		vlik[j] = sumc(log(vf));
	}
	/*----------- Iteration END ------------*/

  
	return vlik;
}

/*
** PPA generate returns and log-RVs from posteriors 
*/
SVAZ::PPA(const nsim)
{
	decl dc, dmu_lambda, ns, msampy;
	decl dphi, dsig, drho, dmu, ddelta, dnu, dxi;
	decl mstate, vh, vz0, vlambda;
	decl veta, vz, veps;

	// setup
	m_flML = 0;
	m_cmp = m_cms = m_cmr = m_cmd = m_can = m_cmn = m_cml = m_cmh = 0; //acceptance counter
	dc = sqrt(2 / M_PI);
	dmu_lambda = 1;
	ns = m_ns;
	msampy = zeros(ns, nsim);
	
	/*--- Set priors ---*/
	if(ismatrix(m_vpri)){
		m_dap0 = m_vpri[0];		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = m_vpri[1];
		m_dn0  = m_vpri[2];		//sigma^2 ~ InvGamma(n0/2, S0/2)
		m_dS0  = m_vpri[3];
		m_dar0 = m_vpri[4];			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = m_vpri[5];
		m_dmm0 = m_vpri[6];		//mu ~ N(mm0, sm0^2)
		m_dsm0 = m_vpri[7];	
		m_dad0 = m_vpri[8];		//(delta+1)/2 ~ Beta(ad0, bd0)
		m_dbd0 = m_vpri[9];
		m_dnn0  = m_vpri[10];		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = m_vpri[11];
	}
	else{
		m_dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = 1.5;
		m_dn0  = 5;		//sigma^2 ~ InvGamma(n0/2, S0/2)
		m_dS0  = 0.2;
		m_dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = 2;
		m_dmm0 = dmu;		//mu ~ N(mm0, sm0^2)
		m_dsm0 = 1;	
		m_dad0 = 1;		//(delta+1)/2 ~ Beta(ad0, bd0)
		m_dbd0 = 1;
		m_dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = 0.5;
	}

	// posterior means of parameters
	m_vposm = loadmat(sprint(m_sdirectory, "/posstat_",m_sfilename,".csv"))[][0]; //posterior means
	dphi = m_vposm[0];
	dsig = m_vposm[1];
	drho = m_vposm[2];
	dmu = m_vposm[3];
	ddelta = m_fls ? m_vposm[4] : 0;
	dnu = m_flt ? m_vposm[4+m_fls] : M_INF_POS;

	// posterior means of state variables
	mstate = loadmat(sprint(m_sdirectory, "/posh_",m_sfilename,".csv"));
	vh = mstate[][0];
	vz0 = mstate[][2];
	vlambda = mstate[][3];
	if(!m_fls)
		vz0 = ones(m_ns,1) * dc;

	// sampling from posteriors
	for(m_k=-0.1*nsim; m_k<nsim; m_k++)
	{
		// sampling parameters and state variables
		dmu = SampMu(dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		dphi = SampPhi(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		if(m_flsigrho)
			[dsig,drho] = SampSigRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		else{
			drho = SampRho(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			dsig = SampSigma(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		}
		if(m_fls) //skew distribution
			ddelta = SampDelta(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		if(m_flt){ // if Student's t
			dnu = SampNu(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
			dmu_lambda = dnu / (dnu-2);
		}
		if(m_fls) // if skewed distribution
			vz0 = SampZ0(dmu, dphi, dsig, drho, ddelta, dnu, vh, vlambda, m_vy);
		vh = SingleMove(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);
		if(m_flt) // if Student's t
			vlambda = SampLambda(dmu, dphi, dsig, drho, ddelta, dnu, vz0, vh, vlambda, m_vy);

		// sampling returns and log-RVs
		if(m_k >= 0){
			veta = vh[1:] - dmu - dphi*(vh[0:ns-2]-dmu) | 0;
			vz = drho*veta/dsig + sqrt(1-drho^2)*rann(ns,1);
			veps = (ddelta.*(vz0-dc) + sqrt(1-ddelta^2).*vz) ./ sqrt(1 - dc^2*ddelta^2);
			msampy[][m_k] = veps .* sqrt(vlambda/dmu_lambda) .* exp(vh/2);
		}
	}

	// save returns and log-RVs
	savemat(sprint(m_sdirectory, "/ppa_sampy_",m_sfilename,".csv"), msampy);
}

/*
**  fTsvar & fGeweke return the diagnostics for convergence
*/
SVAZ::fTsvar(const vx, const iBm)
{
	decl ns, vsp;

	ns = rows(vx);
	vsp = periodogram(vx, iBm, 2, 1) / ns;
	
	return M_2PI * vsp[0];
}

SVAZ::fGeweke(const vx, const iBm)
{
	decl ns, n1, n2, vx1, vx2, dx1bar, dx2bar;
	decl dvar1, dvar2, dz;

	ns = rows(vx);
	n1 = floor(0.1*ns);
	n2 = floor(0.5*ns);
	vx1 = vx[:n1-1];
	vx2 = vx[ns-n2:];
	dx1bar = meanc(vx1);
	dx2bar = meanc(vx2);
	dvar1 = fTsvar(vx1, iBm);
	dvar2 = fTsvar(vx2, iBm);
	dz = (dx1bar - dx2bar) / sqrt(dvar1 / n1 + dvar2 / n2);
	
	return 2 * tailn(fabs(dz));
}
