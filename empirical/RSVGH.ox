/*----------------------------------------------------------
                        RSVGH.ox
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
**  MCMC estimation for realized stochastic volatility model
**  with generalized hyperbolic skew t-distribution
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<RSVGH.h>

RSVGH::RSVGH(const mdata)	//vdata=0 means simulation study
{												   
	m_iRun = ismatrix(mdata) ? 1 : 0;
						//0:simulation, 1:empirical
	m_my = m_myorg = mdata;		  

	if(m_iRun){
		m_ns = rows(mdata);
		m_nKm = int(m_ns/25)-1;	//# of blocks in multi-move sampler
	}

	m_flt = 1;			//model specification: t
	m_flb = 1;			//model specification: skew
	m_flbc = 1;			//bias correction: constant and slope
	m_vini = 0;			//initial value
	m_ncount = 1000;	//counter of iteration
	m_flML = 0;			//flag for ML estimation
	m_flf  = 0;			//flag for VaR and ES backtest
}

/*
** SetFilename set main name for output files
*/
RSVGH::SetFilename(const sfilename)
{
	m_sfilename = sfilename;
}

/*
** SetDirectory set directory for output files
*/
RSVGH::SetDirectory(const sdirectory)
{
	m_sdirectory = sdirectory;
}

/*
** SetDist set distribution of return
**   iDist:  1(RSVskt) / 0(RSVt) / -1(RSVn)
*/
RSVGH::SetDist(const iDist)
{
	m_flt = (iDist>=0);
	m_flb = (iDist==1);
}

/*
** SetBC set bias correction
**   iBC:  0(only constant) / 1(constant and slope)
*/
RSVGH::SetBC(const iBC)
{
	m_flbc = iBC;
}

/*
** SetMultiMove set sampling method for latent volatilities h
*/
RSVGH::SetMultiMove(const flmultimove)
{
	m_flmultimove = flmultimove; //0:single-move; 1:multi-move
}

/*
** SetRanseed set ranseed for random generater
*/
RSVGH::SetRanseed(const irs)
{
	ranseed(irs);
}


/*
** SetBlocks set # of blocks in multi-move sampler
*/
RSVGH::SetBlocks(const nKm)
{
	m_nKm = nKm-1;
}

/*
** SetCount set counter of iteration
*/
RSVGH::SetCount(const ncount)
{
	m_ncount = ncount;
}

/*
** SetForecast set volatility, VaR and ES forecasts
**   iForecast:  1(on) / 0(off)
*/
RSVGH::SetForecast(const iForecast)
{
	m_flf = iForecast;
}

/*
** SetPar set true parameter for simulation
*/
RSVGH::SimObs(const ns, const dphi, const dsigeta, const drho,
			 const dmu, const dbeta_a, const dnu_a, const dxi, const dpsi_a, const dsigu)
{
	decl mSigsq, my, vh, me, vz, dmuz, dsigz, vzb, dbeta, dnu, dpsi, db, i;
	
	m_ns = ns;

	mSigsq = choleski( (1 ~ 0 ~ drho*dsigeta) | (0 ~ dsigu^2 ~ 0) | (drho*dsigeta ~ 0 ~ dsigeta^2));
	
	my = zeros(ns, 2);
	vh = zeros(ns, 1);
	me = mSigsq * rann(3,ns);
	
	vh[0] = dmu + rann(1, 1) * dsigeta / sqrt(1-dphi^2);
	for(i=0 ; i<ns-1; i++)
	{
		vh[i+1] = dmu + dphi * (vh[i] - dmu) + me[2][i];
	}

	if(!m_flt)	// normal
	{
		dnu = M_INF_POS;
		vz = ones(ns, 1);
		dmuz = 1;
		dsigz = 0;
	}
	else	// (skewed) t
	{
		dnu = dnu_a;
		vz = 1 ./ rangamma(ns, 1, dnu/2, dnu/2);
		dmuz = dnu / (dnu - 2);
		dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	}
	vzb = vz - dmuz;

	dbeta = m_flb ? dbeta_a : 0;
	dpsi = m_flbc ? dpsi_a : 1;

	db = dbeta^2*dsigz^2+dmuz;

	my[][0] = (dbeta*vzb + sqrt(vz).*me[0][]')/sqrt(db) .* exp(vh/2);
	my[][1] = dxi + dpsi*vh + me[1][]';

	m_my = m_myorg = my;
	m_vh_t = vh;
	m_vz_t = vz;	
	
	print("\n>>> Data is simualted. [RSV");
	if(m_flbc)
		print("C");
	if(m_flt)
	{
		if(m_flb) print("sk");
		print("t");
	}
	else
		print("n");
	print(" model]"
		  "\nTrue value:\n phi    = ", dphi,
		  "\n sigeta = ", dsigeta, "\n rho    = ", drho,
		  "\n mu     = ", dmu,  "\n beta   = ", dbeta,
		  "\n nu     = ", dnu, "\n xi     = ", dxi, "\n psi    = ", dpsi, "\n sigu   = ", dsigu,
		  "\n mu_z   = ", dmuz,
		  "\n sig_z  = ", dsigz, "\n b      = ", db);
	println("");

	savemat(sprint(m_sdirectory, "/simdata_", m_sfilename, ".csv"), my~vh~vz);
	savemat(sprint(m_sdirectory, "/simpar_", m_sfilename, ".csv"), dphi|dsigeta|drho|dmu|dbeta|dnu|dxi|dpsi|dsigu);

	m_nKm = int(ns/25);
}

/*
** SetIni set initial value for MCMC sampling
*/
RSVGH::SetIni(const dphi, const dsigeta, const drho,
			 const dmu, const dbeta, const dnu, const dxi, const dpsi, const dsigu, const vh)
{
	m_vini = dphi ~ dsigeta ~ drho ~ dmu ~ dbeta ~ dnu ~ dxi ~ dpsi ~ dsigu;
	m_vhini = vh;
}

/*
** SetPri set hyper parameters for priors
*/
RSVGH::SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0,
	const db0, const dsb02, const dnn0, const dnS0, const dmx0, const dsx0, const dmp0, const dsp0, const dsun0,
	const dsuS0, const flprint)
{
	m_vpri = dap0 ~ dbp0 ~ dn0 ~ dS0 ~ dar0 ~ dbr0 ~ dmm0 ~ dsm0 ~ db0 ~ dsb02 ~ dnn0 ~ dnS0
			 ~ dmx0 ~ dsx0 ~ dmp0 ~ dsp0 ~ dsun0 ~ dsuS0;
	if(flprint)
	{
		print("\n");
		print("Prior: (phi+1)/2      ~ Beta(", dap0, ", ", dbp0, ")\n");
		print("       sigma_eta^(-2) ~ Gamma(", dn0, ", ", dS0, ")\n");
		print("       (rho+1)/2      ~ Beta(", dar0, ", ", dbr0, ")\n");
		print("       mu             ~ N(", dmm0, ", ", dsm0^2, ")\n");
		if(m_flt)
		{
			if(m_flb)
				print("       beta           ~ N(", db0, ", ", dsb02, ")\n");
			print("       nu             ~ Gamma(", dnn0, ", ", dnS0, ")\n");
			print("       xi             ~ N(", dmx0, ", ", dsx0^2, ")\n");
			if(m_flbc)
				print("       psi            ~ N(", dmp0, ", ", dsp0^2, ")\n");			
			print("       sigma_u^(-2)   ~ Gamma(", dsun0, ", ", dsuS0, ")\n");
			print("%c", {"Mean", "Stdev"},
				  "%r", {"       phi", "       sigma_eta^2", "       rho", "       nu", "       sigma_u^2"}, "%cf", {"%9.2f"},
				  fBetaMoment(dap0,dbp0,1)*2-1 ~ sqrt(fBetaMoment(dap0,dbp0,2))*2 |
				  fInvGammaMoment(dn0,1/dS0,1) ~ sqrt(fInvGammaMoment(dn0,1/dS0,2)) |
				  fBetaMoment(dar0,dbr0,1)*2-1 ~ sqrt(fBetaMoment(dar0,dbr0,2))*2 |
				  fGammaMoment(dnn0,dnS0,1) ~ sqrt(fGammaMoment(dnn0,dnS0,2)) |
				  fInvGammaMoment(dsun0,1/dsuS0,1) ~ sqrt(fInvGammaMoment(dsun0,1/dsuS0,2)));
		}
		else
		{
			print("       xi           ~ N(", dmx0, ", ", dsx0^2, ")\n");
			if(m_flbc)
				print("       psi            ~ N(", dmp0, ", ", dsp0^2, ")\n");			
			print("       sigma_u^(-2) ~ Gamma(", dsun0, ", ", dsuS0, ")\n");
			print("%c", {"Mean", "Stdev"},
				  "%r", {"       phi", "       sigma_eta^2", "       rho", "       sigma_u^2"}, "%cf", {"%9.2f"},
				  fBetaMoment(dap0,dbp0,1)*2-1 ~ sqrt(fBetaMoment(dap0,dbp0,2))*2 |
				  fInvGammaMoment(dn0,1/dS0,1) ~ sqrt(fInvGammaMoment(dn0,1/dS0,2)) |
				  fBetaMoment(dar0,dbr0,1)*2-1 ~ sqrt(fBetaMoment(dar0,dbr0,2))*2 |
				  fInvGammaMoment(dsun0,1/dsuS0,1) ~ sqrt(fInvGammaMoment(dsun0,1/dsuS0,2)));
		}
	}
}

/*
** DrawData draw graph of data (and state if m_iRun=0)
*/
RSVGH::DrawData()
{
	Draw(0, m_myorg[][0]');
	DrawTitle(0, "Daily return");
	Draw(1, m_myorg[][1]');
	DrawTitle(1, "Log realized volatility");
	if(!m_iRun){
		Draw(2, m_vh_t');
		DrawTitle(2, "Simulated log volatility ($h_t$)");
	}
	DrawDensity(2+!m_iRun, m_myorg[][0]', 1, 1, 0, 1, 1);
	DrawTitle(2+!m_iRun, "Density of daily return");	
	DrawTitle(3+!m_iRun, "QQ-plot to normal density");	
	DrawDensity(4+!m_iRun, m_myorg[][1]', 1, 1, 0, 1, 1);
	DrawTitle(4+!m_iRun, "Density of log realized volatility");	
	DrawTitle(5+!m_iRun, "QQ-plot to normal density");	
  	SaveDrawWindow(sprint(m_sdirectory, "/datadraw_",m_sfilename,".eps"));
  	CloseDrawWindow();
}


/*
** MCMC computes MCMC estimation
*/

RSVGH::MCMC(const nsim, const nburn)
{
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu, i, ncyc;
	decl vh, dmuz, dsigz, db, vz, vzb, vy1b0, npmt, msamp, vsamph, vsampz, vsampvol, msamph;
	decl vsamp, iBm, iacf, aspar, aspar2, mpost, mposh, ascf, vaccptprob;
	decl msampy, msamph1, dh1, dz1, vprob, vvar, mfly, ves, vvolf;
	decl time;

	if(m_ncount>0) time = timer();
	
	/*--- set initial ---*/

	if(!m_flML)
	{
		if(ismatrix(m_vini))
		{
			dphi     = m_vini[0];
			dsigeta  = m_vini[1];
			drho     = m_vini[2];
			dmu      = m_vini[3];
			dbeta    = m_vini[4];
			dnu      = m_vini[5];
			dxi   	 = m_vini[6];
			dpsi	 = m_vini[7];
			dsigu    = m_vini[8];
		}
		else
		{
			dphi = 0.95;
			dsigeta = 0.2;
			drho = -0.5;
			dmu = 0.5;
			dbeta = 0;
			dnu = 10;
			dxi = 0;
			dpsi = 1;
			dsigu = 0.1;
		}
	}
	else
	{
		dphi    = m_vposm[0];
		dsigeta = m_vposm[1];
		drho    = m_vposm[2];
		dmu     = m_vposm[3];
		dbeta   = m_flb ? m_vposm[4] : 0;
		dnu     = m_flt ? m_vposm[4+m_flb] : M_INF_POS;
		dxi     = m_vposm[4+m_flt+m_flb];
		dpsi    = m_flbc ? m_vposm[4+m_flt+m_flb+1] : 1;
		dsigu   = m_vposm[5+m_flt+m_flb+m_flbc];
//		println(dbeta ~ dnu ~ dpsi);
	}

	/*--- prior ---*/
	if(ismatrix(m_vpri))
	{
		m_dap0 = m_vpri[0];		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = m_vpri[1];
		m_dn0  = m_vpri[2];		//sigma^{-2} ~ Gamma(n0, S0)
		m_dS0  = m_vpri[3];
		m_dar0 = m_vpri[4];			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = m_vpri[5];
		m_dmm0 = m_vpri[6];		//mu ~ N(mm0, sm0^2)
		m_dsm0 = m_vpri[7];	
		m_db0   = m_vpri[8];		//beta ~ N(b0, sigb20)
		m_dsb02 = m_vpri[9];
		m_dnn0  = m_vpri[10];		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = m_vpri[11];
		m_dmx0  = m_vpri[12];		//xi ~ N(mx0, sx0^2)
		m_dsx0  = m_vpri[13];
		m_dmp0  = m_vpri[14];		//psi ~ N(mp0, sp0^2)
		m_dsp0  = m_vpri[15];
		m_dsun0	= m_vpri[16];		//sigma_u^{-2} ~ Gamma(sun0, suS0)
		m_dsuS0 = m_vpri[17];
	}
	else
	{
		m_dap0 = 20;		//(phi+1)/2 ~ Beta(ap0, bp0)
		m_dbp0 = 1.5;
		m_dn0  = 2.5;		//sigma^{-2} ~ Gamma(n0, S0)
		m_dS0  = 0.025;
		m_dar0 = 1;			//(rho+1)/2 ~ Beta(ar0, br0)
		m_dbr0 = 2;
		m_dmm0 = dmu;		//mu ~ N(mm0, sm0^2)
		m_dsm0 = 1;	
		m_db0   = 0;		//beta ~ N(b0, sigb20)
		m_dsb02 = 1;
		m_dnn0  = 5;		//nu ~ Gamma(nn0, nS0)
		m_dnS0  = 0.5;
		m_dmx0  = 0;		//xi ~ N(mx0, sx0^2)
		m_dsx0  = 1;
		m_dmp0  = 1;		//psi ~ N(mp0, sp0^2)
		m_dsp0  = 1;
		m_dsun0	= 2.5;		//sigma_u^{-2} ~ Gamma(sun0, suS0)
		m_dsuS0 = 0.1;
	}
	
	if(!m_flb) dbeta = 0;
	if(!m_flt) dnu = M_INF_POS;
	if(!m_flbc) dpsi = 1;
								
	if(!m_flML)
	{
		if (ismatrix(m_vhini)) vh = m_vhini;
		else vh = dmu + rann(m_ns, 1) * dsigeta;
	}
	else vh = loadmat(sprint(m_sdirectory, "/posh_",m_sfilename,".csv"))[][0];

	if(m_flt)
	{
		vz = 1 ./ rangamma(m_ns, 1, dnu/2, dnu/2);
		dmuz = dnu / (dnu - 2);
		dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	}
	else
	{
		vz = ones(m_ns, 1);
		dmuz = 1;
		dsigz = 0;
	}
	db = dbeta^2*dsigz^2 + dmuz;
	vzb = vz - dmuz;

	if(!m_iRun)
	{
		vh = m_vh_t;
		vz = m_vz_t;
		if(m_flt)
		{
			dmuz = dnu / (dnu - 2);
			dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
		}
		else
		{
			dmuz = 1;
			dsigz = 0;
		}
		db = dbeta^2*dsigz^2 + dmuz;
		vzb = vz - dmuz;
	}

	/*--- set sample option ---*/
//	nburn = 0.1 * nsim;	//burn-in period
	npmt = 6+m_flt+m_flb+m_flbc;	//# of parameter

	
	if(!m_flML)	/*--------- parameter estimation --------*/
	{
		msamp = zeros(nsim, npmt);	//sample box
		vsamph = vsampvol = vsampz = zeros(m_ns, 1);
		msamph = zeros(nsim, m_ns);
		msampy = zeros(nsim, 3);
		msamph1 = zeros(nsim, 2);
		m_cmp = m_cms = m_cmn = m_cmz = m_cmb = m_cmh = m_cah = 0; //acceptance counter
		m_cnu1 = m_cnu2 = 0; //counter of breaks in loop for sampling nu
	
		/*--- MCMC sampling ---*/
		if(m_ncount > 0) println("\n\nIteration:");
	
		/*----------- S A M P L I N G   S T A R T -----------*/
		for(m_k=-nburn ; m_k<nsim ; m_k++)
		{
			vy1b0 = (sqrt(db)*m_my[][0].*exp(-vh/2) - dbeta*vzb) ./ sqrt(vz);
		
			/*--- sampling phi ---*/
			dphi = SampPhi(drho*dsigeta*vy1b0, dphi, dsigeta, drho, vh-dmu);
		
			/*--- sampling sigma and rho ---*/
			[dsigeta, drho] = SampSigRho(dphi, dsigeta, drho, vh-dmu, vy1b0);
		  
			/*--- sampling mu ---*/
			dmu = SampMu(dphi, dsigeta, drho, dbeta, drho*dsigeta*vy1b0, vh);
		
			if(m_flt)	// if t
			{
				/*--- sampling beta ---*/
				if(m_flb)
					dbeta = SampBeta(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz, vzb);
	
				/*--- sampling nu ---*/
				dnu = SampNu(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz);
				dmuz = dnu / (dnu - 2);
				dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
				db = dbeta^2*dsigz^2 + dmuz;
		
				/*--- sampling z ---*/
				vz = SampZ(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz);
				vzb = vz - dmuz;
			}
		
			/*--- sampling xi ---*/
			dxi = SampXi(dpsi, dsigu, vh);
		
			/*--- sampling psi ---*/
			if(m_flbc)
				dpsi = SampPsi(dxi, dsigu, vh);

			/*--- sampling sigma_u ---*/
			dsigu = SampSigu(dxi, dpsi, vh);
								
			/*--- sampling h ---*/
			if(m_flmultimove)
				vh = MultiMove(dphi, dsigeta, drho, exp(dmu/2)/sqrt(db), dbeta, dnu, dxi+dpsi*dmu, dpsi, dsigu, vh-dmu, vz, vzb) + dmu;
			else
				vh = SingleMove(dmu, dphi, dsigeta, drho, dbeta, dnu, dxi, dsigu, vh, vz, m_my);
	
			if(m_k >= 0)
			{
				/*--- storing sample ---*/
				vsamp = dphi ~ dsigeta ~ drho ~ dmu;
				if(m_flb)
					vsamp ~= dbeta;
				if(m_flt)
					vsamp ~= dnu;
				vsamp ~= dxi;
				if(m_flbc)
					vsamp ~= dpsi;
				vsamp ~= dsigu;
				
				msamp[m_k][] = vsamp;
				vsamph += vh;
				vsampvol += exp(vh/2);
				vsampz += vz;
				msamph[m_k][] = vh';
			
				if(m_flf)
				{
					/*--- VaR and Expected Shortfall ---*/
					if(m_flt)
						dz1 = 1 / rangamma(1, 1, dnu/2, dnu/2);
					else
						dz1 = 1;
	
					dh1 = dmu + dphi * (vh[m_ns-1] - dmu)
						+ drho * dsigeta * (sqrt(db)*m_my[m_ns-1][0] * exp(-vh[m_ns-1]/2) - dbeta * (vz[m_ns-1] - dmuz)) / sqrt(vz[m_ns-1])
						+ sqrt(1-drho^2) * dsigeta * rann(1,1);
					msampy[m_k][0] = (dbeta * (dz1-dmuz) + sqrt(dz1) .* rann(1, 1)) / sqrt(db) * exp(dh1/2);	// daily return y_{1,n+1}
					msampy[m_k][1] = dxi + dpsi*dh1 + dsigu*rann(1,1);		// log-realized volatility y_{2,n+1}
					msampy[m_k][2] = exp(msampy[m_k][1]);	// realized volatility RV_{n+1}
					msamph1[m_k][] = dh1 ~ exp(dh1);	// predicted log-voalitlity and volatility
				}
			}
		
			//print counter
			if(m_ncount>0 && !imod(m_k, m_ncount)) println(m_k);
		}		 
	    /*------------ S A M P L I N G   E N D --------------*/
	
		
		/*--- output result ---*/
		
		iBm = iacf = min(1000, nsim/2);	//bandwidth
	
		aspar =	{"phi", "sigma_eta", "rho", "mu"};
		aspar2 = {"$\\phi$", "$\\sigma_\\eta$", "$\\rho$", "$\\mu$"};
		if(m_flt)
		{
			if(m_flb){
			aspar ~= {"beta"};
			aspar2 ~= {"$\\beta$"};
			}
			aspar ~= {"nu"};
			aspar2 ~= {"$\\nu$"};
		}
		aspar ~= {"xi"};
		aspar2 ~= {"$\\xi$"};
		if(m_flbc)
		{
			aspar ~= {"psi"};
			aspar2 ~= {"$\\psi$"};
		}
		aspar ~= {"sigma_u"};
		aspar2 ~= {"$\\sigma_u$"};
		aspar ~= {"h(500)"};
		aspar2 ~= {"$h_{500}$"};
		
		print("\n                       ",
			  "[ESTIMATION RESULT]"
			  "\n------------------------------------------------------------------------"
			  "\nParameter           Mean     Stdev     95%L      95%U   Geweke     Inef."  
			  "\n------------------------------------------------------------------------");
	
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
		print("%r", aspar, "%cf", ascf, mpost|mposh[499][]);
		print("--------------------------------------------------------------- RSVGH.ox\n");
	
		print("Model: RSV");
		if(m_flbc) print("C");
		if(m_flt)
		{
			if(m_flb) print("sk");
			print("t");
		}
		else
			print("n");
		if(m_iRun) print(" (empirical study: n=", m_ns, ")\n");
		else	   print(" (simulation study: n=", m_ns, ")\n");
		print("Iteration: ", nsim, "\n");

		print("\nAcceptance rate for parameter"
			  "\n phi:    ", "%5.1f", m_cmp/nsim*100, "%",
			  "\n sig,rho:", "%5.1f", m_cms/nsim*100, "%");
		vaccptprob = m_cmp/nsim ~ m_cms/nsim;
		if(m_flt)
		{
			print(
				  "\n nu:     ", "%5.1f", m_cmn/nsim*100, "%",
				  " (Break1:", "%5.1f", m_cnu1/nsim*100, "%, "
				  " Break2:", "%5.1f", m_cnu2/nsim*100, "%)"
				  "\n z:      ", "%5.1f", double(m_cmz/(nsim*m_ns)*100), "%");
			vaccptprob ~= m_cmn/nsim ~ m_cmz/(nsim*m_ns);
		}
		if(m_flb)
		{
			print("\n beta:   ", "%5.1f", double(m_cmb/nsim*100), "%");
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
	
		if(m_flf)
		{
			vprob = <.005, .01, .05, .1>;
			vvar = quantilec(msampy[][0], vprob)';
			mfly = msampy[][0] .< vvar;
			ves = msampy[][0] ' mfly ./ sumc(mfly);
			vvolf = meanc(msamph1 ~ msampy[][1:2]);	// (h, exp(h), y_2, exp(y_2))
			print("\n\n1-day-ahead VaR and ES");
			print("%r", {"tail:", " VaR", " ES"}, 
				  "%cf", {"%10.5f"}, vprob | vvar | ves);
			print("\n1-day-ahead volatility");
			print("%c", {"h", "exp(h)", "y_2", "exp(y_2)"}, 
				  "%cf", {"%10.5f"}, vvolf);
			savemat(sprint(m_sdirectory, "/var_",m_sfilename,".csv"), vvar);
			savemat(sprint(m_sdirectory, "/es_",m_sfilename,".csv"), ves);
			savemat(sprint(m_sdirectory, "/volf_",m_sfilename,".csv"), vvolf);
			savemat(sprint(m_sdirectory, "/h1_y_", m_sfilename, ".csv"), msamph1[][0] ~ msampy[][0:1]); //1-day ahead log-volatility, return, log-RV
		}
	
		if(m_ncount > 0)
		{
			//draw MCMC result
			DrawAdjust(ADJ_AREAMATRIX, 3, npmt);
		 	for(i=0 ; i<npmt ; i++){
		 		DrawCorrelogram(i, msamp[][i]', {}, iacf);
		    	DrawTitle(i, aspar2[i]);
		    	DrawMatrix(i+npmt, thinr(msamp[][i], 500)', {""}, 1, nsim/500);
		    	DrawTitle(i+npmt, aspar2[i]);
		    	DrawDensity(i+2*npmt, msamp[][i]', {""}, 1);
		    	DrawTitle(i+2*npmt, aspar2[i]);
		    } 
		  	SaveDrawWindow(sprint(m_sdirectory, "/samp_",m_sfilename,".eps"));
		  	CloseDrawWindow();

	
			//draw state variable
			DrawMatrix(0, m_my[][0]', "", 1, 1);
			DrawTitle(0, "Data: y_1");
			DrawMatrix(1, m_my[][1]', "", 1, 1);
			DrawTitle(1, "Data: y_2");
			if(m_iRun){
			DrawMatrix(2, mposh[][0]', "", 1, 1);
			DrawTitle(2, "Posterior mean of state");
			}	
			else{
			DrawMatrix(2, (m_vh_t ~ mposh[][0])', {"True", "PosMean"}, 1, 1);
			DrawTitle(2, "True and posterior mean of state");
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
		}
		savemat(sprint(m_sdirectory, "/posstat_", m_sfilename, ".csv"), mpost);
		savemat(sprint(m_sdirectory, "/posh_", m_sfilename, ".csv"), vsamph/nsim ~ vsampvol/nsim ~ vsampz/nsim);
		savemat(sprint(m_sdirectory, "/posstath_", m_sfilename, ".csv"), mposh);
		savemat(sprint(m_sdirectory, "/aprob_", m_sfilename, ".csv"), vaccptprob);

	}	//end of parameter estimation
	else	/*------------------------- ML estimation ---*/
	{
		/*--- set variables ---*/
		ncyc = 7+2*m_flt+2*m_flb+m_flbc;	//# of cycle for reduced run
		m_mmls = zeros(nsim, ncyc);
	 
		/*--- MCMC sampling ---*/
		print("\n[Reduced MCMC]\n  Iteration:");
								  
		for(m_kk=0 ; m_kk<ncyc ; m_kk++)
		{
			print(m_kk);
			/*----------- S A M P L I N G   S T A R T -----------*/
			for(m_k=-nburn ; m_k<nsim ; m_k++)
			{
				vy1b0 = (sqrt(db)*m_my[][0].*exp(-vh/2) - dbeta*vzb) ./ sqrt(vz);
			
				/*--- sampling phi ---*/
				if(m_kk < 2)
					dphi = SampPhi(drho*dsigeta*vy1b0, dphi, dsigeta, drho, vh-dmu);
			
				/*--- sampling sigma and rho ---*/
				if(m_kk < 4)
					[dsigeta, drho] = SampSigRho(dphi, dsigeta, drho, vh-dmu, vy1b0);
			  
				/*--- sampling mu ---*/
				if(m_kk < 5)
					dmu = SampMu(dphi, dsigeta, drho, dbeta, drho*dsigeta*vy1b0, vh);
			
				if(m_flt)	// if t
				{
					/*--- sampling beta ---*/
					if(m_flb && m_kk < 7)
						dbeta = SampBeta(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz, vzb);
		
					/*--- sampling nu ---*/
					if(m_kk < 7+2*m_flb)
					{
						dnu = SampNu(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz);
						dmuz = dnu / (dnu - 2);
						dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
						db = dbeta^2*dsigz^2 + dmuz;
					}
			
					/*--- sampling z ---*/
					vz = SampZ(dphi, dsigeta, drho, dmu, dbeta, dnu, vh, vz);
					vzb = vz - dnu / (dnu - 2);
				}
			
				/*--- sampling xi ---*/
				if(m_kk < 6+2*m_flt+2*m_flb)
					dxi = SampXi(dpsi, dsigu, vh);
			
				/*--- sampling psi ---*/
				if(m_flbc && m_kk < 7+2*m_flt+2*m_flb)
					dpsi = SampPsi(dxi, dsigu, vh);

				/*--- sampling sigma_u ---*/
				if(m_kk < 7+2*m_flt+2*m_flb+m_flbc)
					dsigu = SampSigu(dxi, dpsi, vh);
									
				/*--- sampling h ---*/
				if(m_flmultimove)
					vh = MultiMove(dphi, dsigeta, drho, exp(dmu/2)/sqrt(db), dbeta, dnu, dxi+dpsi*dmu, dpsi, dsigu, vh-dmu, vz, vzb) + dmu;
				else
						vh = SingleMove(dmu, dphi, dsigeta, drho, dbeta, dnu, dxi, dsigu, vh, vz, m_my);

				if(m_ncount>0 && !imod(m_k, m_ncount)) print(" ",m_k/m_ncount);
			}		 
		    /*------------ S A M P L I N G   E N D --------------*/
			if(m_ncount>0) print("\n  Iteration:");
	
			/*--- fix parameter ---*/
			if(m_kk == 0)
				dphi = m_vposm[0];
			else if(m_kk == 2)
			{
				dsigeta = m_vposm[1];
				drho = m_vposm[2];
			}
			else if(m_kk == 4)
				dmu = m_vposm[3];
		
			if(m_flt)
			{
				if(m_flb && m_kk == 5)
					dbeta = m_vposm[4];
				if(m_kk == 5+2*m_flb)
					dnu = m_vposm[4+m_flb];
			}
	
			if(m_kk == 5+2*m_flt+2*m_flb)
				dxi = m_vposm[4+m_flt+m_flb];
			if(m_flbc && m_kk == 6+2*m_flt+2*m_flb)
				dpsi = m_vposm[5+m_flt+m_flb];
			if(m_kk == 6+2*m_flt+2*m_flb+m_flbc)
				dsigu = m_vposm[5+m_flt+m_flb+m_flbc];
						
		}	/*------------------------------ end of cycle ---*/
	}	//end of ML estimation

	if(m_ncount>0 && !m_flML) println("\n\nTime: ", timespan(time));
}


/*
** SampPhi samples parameter phi by MH algorithm 
*/
RSVGH::SampPhi(const vy1b, const dphio, const dsigeta, const drho, const vhb)
{
	decl vmu, dsum, dmup, dsigp, dphin, dfrac, dphia, dqa, nite;

	/*--- draw candidate ---*/

	dsum = drho^2 * vhb[0]^2 + sumsqrc(vhb[1:m_ns-2]);

	dmup = sumc((vhb[1:]-vy1b[:m_ns-2]).*vhb[:m_ns-2])/dsum;
	dsigp = sqrt(dsigeta^2*(1-drho^2) / dsum);

	nite = 0;
	do
	{
		nite++;
		dphin = dmup + rann(1, 1) * dsigp;
		if(nite>10000)
		{
			print("phi! ");
			dphin = dphio;
		}
	}
	while(fabs(dphin) >= 1 && nite<=10000);
	

	/*--- accept candidate ---*/

	dfrac = densbeta((dphin+1)/2, m_dap0, m_dbp0)
		  / densbeta((dphio+1)/2, m_dap0, m_dbp0)
		  * sqrt(1-dphin^2) / sqrt(1-dphio^2);


	if(!m_flML)	//parameter estimation
	{
		if(ranu(1,1) < dfrac)
		{
			if(m_k>=0)
				m_cmp++;	//acceptance counter
		}
		else
			dphin = dphio;		//reject
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
				dfrac = densbeta((dphia+1)/2, m_dap0, m_dbp0)
					  / densbeta((dphin+1)/2, m_dap0, m_dbp0)
					  * sqrt(1-dphia^2) / sqrt(1-dphin^2);
	
				dqa = densn((dphia-dmup)/dsigp)/dsigp / (probn((1-dmup)/dsigp) - probn((-1-dmup)/dsigp));
			  
				m_mmls[m_k][0] = min(1, dfrac) * dqa;
			}
		}
		else if(m_k>=0)	//denominator
		{
			dfrac = densbeta((dphin+1)/2, m_dap0, m_dbp0)
				  / densbeta((dphia+1)/2, m_dap0, m_dbp0)
				  * sqrt(1-dphin^2) / sqrt(1-dphia^2);
	
			m_mmls[m_k][1] = min(1, dfrac);
			
			dphin = dphia;
		}
	}
  
	return dphin;
}

/*
** fPosSigRho computes log posterior distribution
** of sigma and rho 
*/
RSVGH::fPosSigRho(const vp, const dfunc, const vsc,
				 const mhe)
{
	decl dsig, drho;

	dsig = exp(vp[0]);
	drho = 2 / (1+exp(-vp[1])) - 1;

	dfunc[0]
	  = log(1/dsig^4 * densgamma(1/dsig^2, m_dn0, m_dS0))
	  + log(densbeta((drho+1)/2, m_dar0, m_dbr0))
	  - m_ns * log(dsig)
	  - 0.5 * (m_ns - 1) * log(1 - drho^2)
	  - 0.5 * sumsqrc(m_vh[1:] - m_dphi*m_vh[:m_ns-2] - dsig*drho*m_vy1b0[:m_ns-2]) / (dsig^2*(1 - drho^2))
	  - 0.5 * (1 - m_dphi^2) * m_vh[0]^2 / dsig^2;
							   
	return 1;
}

/*
** SampSigRho samples parameter sigma and rho
** by MH algorithm 
*/
RSVGH::SampSigRho(const dphi, const dsigetao, const drhoo, const vhb, const vy1b0)
{
	decl vparo, vparn, vparm, dlm, vsc, mhe, mS, vm;
	decl dlo, dln, dho, dhn, dsigetan, drhon, dfrac;
	decl dsigetaa, drhoa, vpara, dla, dha, dja, dqa;

	m_vh = vhb;
	m_dphi = dphi;
	m_vy1b0 = vy1b0;

	vparo = vparm
		  = log(dsigetao)
		  | log((1+drhoo)/(1-drhoo));	//transformed

	MaxBFGS(fPosSigRho, &vparm, &dlm, 0, 1);//find mode
			
	fPosSigRho(vparm, &dlm, 0, 0);
  	Num1Derivative(fPosSigRho, vparm, &vsc);
  	Num2Derivative(fPosSigRho, vparm, &mhe);
	  
	mS = -invert(mhe);
	vm = vparm + mS * vsc;

	vparn = vm + choleski(mS) * rann(2, 1);	//proposal
				
	fPosSigRho(vparo, &dlo, 0, 0);
	fPosSigRho(vparn, &dln, 0, 0);
												    
	dsigetan = exp(vparn[0]);
	drhon = 2 / (1 + exp(-vparn[1])) - 1;
						  
	dho = dlm + (vparo-vparm)'vsc
		+ 0.5*(vparo-vparm)'mhe*(vparo-vparm)
		- log(dsigetao)
		+ log(2) - log(1 - drhoo^2);
	dhn = dlm + (vparn-vparm)'vsc
		+ 0.5*(vparn-vparm)'mhe*(vparn-vparm)
		- log(dsigetan)
		+ log(2) - log(1 - drhon^2);

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
			dsigetan = dsigetao;
			drhon = drhoo;
		}
	}
	else	//ML estimation: generate [sig,rho] when m_kk = 0,1,2,3
	{
		if(m_kk < 2)
		{
			if(ranu(1,1) > dfrac)	//reject case
			{
				dsigetan = dsigetao;
				drhon = drhoo;
			}
		}
		else
		{
			dsigetaa = m_vposm[1];
			drhoa = m_vposm[2];
			vpara = log(dsigetaa) | log((1+drhoa)/(1-drhoa));	//transformed
			fPosSigRho(vpara, &dla, 0, 0);
	
			dja = - log(dsigetaa) + log(2) - log(1 - drhoa^2);
			dha = dlm + (vpara - vparm) ' vsc
				+ 0.5 * (vpara - vparm) ' mhe * (vpara - vparm)
				+ dja;
			dqa = -0.5 * (vpara - vm) ' invert(mS) * (vpara - vm)
				- log(M_2PI * sqrt(determinant(mS)))
				+ dja;
		
			if(m_kk == 2)	//numerator
			{
				if(ranu(1,1) > dfrac)	//reject case
				{
		 			dsigetan = dsigetao;
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
		
				dsigetan = dsigetaa;
				drhon = drhoa;
			}
		}
	}
  				  
	return {dsigetan, drhon};
}
  
/*
** SampMu samples parameter mu
*/
RSVGH::SampMu(const dphi, const dsigeta, const drho, const dbeta, const vy1b, const vh)
{
	decl vdif, dsigm2, dmuh, dmu, dmua;

	vdif = vh[1:] - dphi*vh[:m_ns-2] - vy1b[:m_ns-2];
	
	dsigm2 = 1/ (1 / m_dsm0^2
				 + ((1-drho^2)*(1-dphi^2)
					+ (m_ns-1)*(1-dphi)^2)
					/ (dsigeta^2*(1-drho^2)));
	dmuh = dsigm2
		 *(m_dmm0 / m_dsm0^2
		   + ((1-drho^2)*(1-dphi^2)*vh[0]
			  + (1-dphi)*sumc(vdif))/(dsigeta^2*(1-drho^2)));

	if(!m_flML)	//parameter estimation
		dmu = dmuh + sqrt(dsigm2) * rann(1, 1);
	else	//ML estimation: generate mu for m_kk = 0,1,2,3,4
	{
		if(m_kk < 4)
			dmu = dmuh + sqrt(dsigm2) * rann(1, 1);
		else
		{
			dmu = dmua = m_vposm[3];
			if(m_k>=0)
				m_mmls[m_k][4]
				  = densn((dmua-dmuh)/sqrt(dsigm2))
				  / sqrt(dsigm2);
		}		
	}

	return dmu;
}

/*
** fPosBeta computes log posterior distribution of beta
*/
RSVGH::fPosBeta(const dp, const dfunc, const dsc, const dhe)
{
	decl dbeta, dphi, dsigeta, drho, dmu, dnu, vh, vz, dmuz, dsigz, db, vzb, vhbd, vy1b;

	dbeta = dp;
	dphi = m_dphi;
	dsigeta = m_dsigeta;
	drho = m_drho;
	dmu = m_dmu;
	dnu = m_dnu;
	vh = m_vh;
	vz = m_vz;

	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = dbeta^2*dsigz^2+dmuz;
	vzb = vz - dmuz;
	vhbd = (vh[1:] - dmu) - dphi * (vh[:m_ns-2] - dmu);
	vy1b = (sqrt(db)*m_my[][0].*exp(-vh/2) - dbeta*vzb) ./ sqrt(vz);
	
	dfunc[0]
	  = - 0.5*(dbeta-m_db0)^2/m_dsb02
	  	+ 0.5*m_ns*log(db)
		- 0.5*sumsqrc(vy1b)
		- 0.5*sumsqrc(vhbd-drho*dsigeta*vy1b[:m_ns-2])/((1-drho^2)*dsigeta^2);
				   
	return 1;
}

/*
** SampBeta samples parameter beta
*/
RSVGH::SampBeta(const dphi, const dsigeta, const drho, const dmu, const dbetao, const dnu, const vh, const vz, const vzb)
{
	decl dbetam, dlm, dsc, dhe, dS, dm, dbetan, dlo, dln, dho, dhn, dfrac, dbetaa, dla, dha, dqa;

	m_dphi = dphi;
	m_dsigeta = dsigeta;
	m_drho = drho;
	m_dmu = dmu;
	m_dnu = dnu;
	m_vh = vh;
	m_vz = vz;

	dbetam = dbetao; //initial value

	MaxBFGS(fPosBeta, &dbetam, &dlm, 0, 1); //find mode
			
	fPosBeta(dbetam, &dlm, 0, 0);
  	Num1Derivative(fPosBeta, dbetam, &dsc);
  	Num2Derivative(fPosBeta, dbetam, &dhe);
	  
	dS = -1 / dhe;
	dm = dbetam + dS * dsc;

	dbetan = dm + rann(1,1)*sqrt(dS);	//proposal
	  
	fPosBeta(dbetao, &dlo, 0, 0);
	fPosBeta(dbetan, &dln, 0, 0);
												    
	dho = dlm + (dbetao-dbetam)*dsc + 0.5*dhe*(dbetao-dbetam)^2;
	dhn = dlm + (dbetan-dbetam)*dsc + 0.5*dhe*(dbetan-dbetam)^2;

	dfrac = exp(dln - dhn - dlo + dho);
		
	if(!m_flML)	//parameter estimation
	{
		if(ranu(1,1) <= dfrac)
		{
			if(m_k>=0) m_cmb++;	//acceptance counter
		}
		else
			dbetan = dbetao;	//reject
	}
	else	//ML estimation: generate beta for m_kk = 0,1,2,3,4 if skew
	{
		if(m_kk < 5)
		{
			if(ranu(1,1) > dfrac)
				dbetan = dbetao;	//reject
		}
		else
		{
			dbetaa = m_vposm[4];
			fPosBeta(dbetaa, &dla, 0, 0);
				
			dha = dlm + (dbetaa-dbetam)*dsc + 0.5*dhe*(dbetaa-dbetam)^2;
			dqa = -0.5 * (dbetaa - dm)^2 / dS - 0.5 * log(M_2PI * sqrt(dS));

			if(m_kk == 5)	//numerator
			{
				if(ranu(1,1) > dfrac)	//reject case
				{
		 			dbetan = dbetao;
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
		
				dbetan = dbetaa;
			}
		}
	}
	
	return dbetan;
}

/*
** fPosNut computes log posterior distribution of p=log(nu-4)
*/
RSVGH::fPosNut(const dp, const dfunc, const dsc, const dhe)
{
	decl dnu, dmuz, dsigz, db, vzb, vy1b;

	dnu = exp(dp)+4;
	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = m_dbeta^2*dsigz^2+dmuz;

	vzb = m_vz - dmuz;
	vy1b = (sqrt(db)*m_my[][0].*exp(-m_vh/2) - m_dbeta*vzb) ./ sqrt(m_vz);
	
	dfunc[0]
	  =	dp //Jacobian = |dnu/dp| = exp(p) --> log(Jacobian) = p
	  + log(densgamma(dnu, m_dnn0, m_dnS0))
	  - 0.5 * sumsqrc(vy1b)
	  - 0.5 * sumsqrc(m_vhbd - m_dsigeta*m_drho*vy1b[:m_ns-2])/(m_dsigeta^2 * (1-m_drho^2))
	  + 0.5 * m_ns * dnu * log(dnu / 2)
	  - m_ns * loggamma(dnu / 2)
	  - 0.5 * dnu * sumc(log(m_vz) + 1 ./ m_vz)
	  + 0.5 * m_ns * log(db);
				   
	return 1;
}

/*
** fPosNu computes log posterior distribution of nu
*/
RSVGH::fPosNu(const dp, const dfunc, const dsc, const dhe)
{
	decl dnu, dmuz, dsigz, db, vzb, vy1b;

	dnu = exp(dp);
	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = m_dbeta^2*dsigz^2+dmuz;

	vzb = m_vz - dmuz;
	vy1b = (sqrt(db)*m_my[][0].*exp(-m_vh/2) - m_dbeta*vzb) ./ sqrt(m_vz);
	
	dfunc[0]
	  = log(densgamma(dnu, m_dnn0, m_dnS0))
	  - 0.5 * sumsqrc(vy1b)
	  - 0.5 * sumsqrc(m_vhbd - m_dsigeta*m_drho*vy1b[:m_ns-2])/(m_dsigeta^2 * (1-m_drho^2))
	  + 0.5 * m_ns * dnu * log(dnu / 2)
	  - m_ns * loggamma(dnu / 2)
	  - 0.5 * dnu * sumc(log(m_vz) + 1 ./ m_vz)
	  + 0.5 * m_ns * log(db);
				   
	return 1;
}

/*
** SampNu samples parameter nu by MH algorithm 
*/
RSVGH::SampNu(const dphi, const dsigeta, const drho, const dmu, const dbeta, const dnuo, const vh, const vz)
{
	decl dnuot, dnun, dnunt, dnum, dlm, dsc, dhe, dS, dm;
	decl dfrac, dlo, dln, dho, dhn;
	decl dnua, dnuat, dla, dha, dja, dqa;
	decl dp, nite;

	m_dsigeta = dsigeta;
	m_drho = drho;
	m_dbeta = dbeta;
	m_vh = vh;
	m_vhbd = (vh[1:] - dmu) - dphi * (vh[:m_ns-2] - dmu);
	m_vz = vz;

	dnuot = dnum = log(dnuo);//transformed

	dp = log(dnuo-4); //nu>4
	MaxBFGS(fPosNut, &dp, &dlm, 0, 1); //find mode

	dnum = log(exp(dp)+4);
	fPosNu(dnum, &dlm, 0, 0);
  	Num1Derivative(fPosNu, dnum, &dsc);
  	Num2Derivative(fPosNu, dnum, &dhe);
	  
	dS = -1 / dhe;
	dm = dnum + dS * dsc;

//	if(probn((log(4)-dm)/sqrt(dS))>0.99)
//	{
//		if(m_k>=0 && !m_flML) m_cnu1++;
////		print("nu1! ");
////		print("\nIteration: ","%5.0f",m_k,"\n  Mode=","%6.1f",dm," with SD=","%6.1f",sqrt(dS)," is too small");
//		dnum = log(4.1);
//		fPosNu(dnum, &dlm, 0, 0);
//	  	Num1Derivative(fPosNu, dnum, &dsc);
//	  	Num2Derivative(fPosNu, dnum, &dhe);
//		dS = -1 / dhe;
//		dm = dnum + dS * dsc;
//	}
	nite = 0;
	do{
		nite++;
		dnunt = dm + rann(1, 1) * sqrt(dS);	//proposal
		if(nite>100)
		{
			if(m_k>=0 && !m_flML) m_cnu2++;
//			print("\nStep: ","%5.0f",nite," dm=","%6.1f",dm," dS=",dS,"%6.1f","dnunt=","%6.1f",dnunt);
//			print("nu2! ");
			dnunt = log(4) + sqrt(dS)*fabs(rann(1,1));
		}
	}while(dnunt <= log(4) && nite<=100);
	  
	fPosNu(dnuot, &dlo, 0, 0);
	fPosNu(dnunt, &dln, 0, 0);
												    
	dnun = exp(dnunt);
	
	dho = dlm + (dnuot-dnum)*dsc + 0.5*dhe*(dnuot-dnum)^2
		- log(dnuo);
	dhn = dlm + (dnunt-dnum)*dsc + 0.5*dhe*(dnunt-dnum)^2
		- log(dnun);


	dfrac = exp(dln - dhn - dlo + dho);
						
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
		if(m_kk < 5+2*m_flb)
		{
			if(ranu(1,1) > dfrac)	//reject case
	 			dnun = dnuo;
		}
		else
		{
			dnua = m_vposm[4+m_flb];
			dnuat =	log(dnua);	//transformed
			fPosNu(dnuat, &dla, 0, 0);
				
			dja = - log(dnua);
			dha = dlm + (dnuat-dnum)*dsc + 0.5*dhe*(dnuat-dnum)^2
				+ dja;
			dqa = -0.5 * (dnuat - dm)^2 / dS
				- 0.5 * log(M_2PI * sqrt(dS))
				+ dja
				- log(1 - probn((log(4) - dm) / sqrt(dS)));

			if(m_kk == 5+2*m_flb)	//numerator
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
					m_mmls[m_k][5+2*m_flb] = min(1, dfrac) * exp(dqa);
				}
			}
			else if(m_k>=0)	//denominator
			{
				dfrac = exp(dln - dhn - dla + dha);
				m_mmls[m_k][6+2*m_flb] = min(1, dfrac);
		
				dnun = dnua;
			}
		}
	}

	return dnun;
}

/*
** SampZ samples parameter lambda by MH algorithm 
*/
RSVGH::SampZ(const dphi, const dsigeta, const drho, const dmu, const dbeta, const dnu, const vh, const vzo)
{
	decl dmuz, dsigz, db, vzn, vzbo, vzbn, vy1bo, vy1bn, vlo, vln, vac, vhbd;

	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = dbeta^2*dsigz^2+dmuz;
	vzn  = 1 ./ rangamma(m_ns, 1, (dnu+1)/2, dnu/2); //candidate

	vhbd = (vh[1:] - dmu) - dphi * (vh[:m_ns-2] - dmu);

	vzbo = vzo - dmuz;
	vzbn = vzn - dmuz;
	vy1bo = (sqrt(db)*m_my[][0].*exp(-vh/2) - dbeta*vzbo) ./ sqrt(vzo);
	vy1bn = (sqrt(db)*m_my[][0].*exp(-vh/2) - dbeta*vzbn) ./ sqrt(vzn);
						   
	vlo = -0.5 * vy1bo.^2
		- 0.5 * (((vhbd - dsigeta*drho*vy1bo[:m_ns-2]).^2) | 0) / (dsigeta^2 * (1 - drho^2));
	vln = -0.5 * vy1bn.^2
		- 0.5 * (((vhbd - dsigeta*drho*vy1bn[:m_ns-2]).^2) | 0) / (dsigeta^2 * (1 - drho^2));
 						   
	vac = ranu(m_ns, 1) .< exp(vln - vlo); 
	vzn = vac .? vzn .: vzo;

	if(!m_flML && m_k>=0) m_cmz += sumc(vac);
						
	return vzn;
}

/*
** SampXi samples parameter xi
*/
RSVGH::SampXi(const dpsi, const dsigu, const vh)
{
	decl dden, vdif, dmx, dsx, dxi, dxia;

	dden = m_ns*m_dsx0^2 + dsigu^2;
	vdif = m_my[][1] - dpsi*vh;

	dmx = (m_dsx0^2 * sumc(vdif) + m_dmx0*dsigu^2) / dden;
	dsx = sqrt(m_dsx0^2 * dsigu^2 / dden);

	if(!m_flML)	//parameter estimation
		dxi = dmx + dsx * rann(1,1);	// xi|. ~ N(mx, sx^2)
	else	//ML estimation: generate xi for m_kk = 0,1,...,5+2*m_flt+m_flb
	{
		if(m_kk < 5+2*m_flt+2*m_flb)
			dxi = dmx + dsx * rann(1,1);	// xi|. ~ N(mx, sx^2)
		else
		{
			dxi = dxia = m_vposm[4+m_flt+m_flb];
			if(m_k>=0)
				m_mmls[m_k][5+2*m_flt+2*m_flb] = densn((dxia-dmx)/dsx) / dsx;
		}		
	}
	
	return dxi;
}

/*
** SampXi samples parameter xi
*/
RSVGH::SampPsi(const dxi, const dsigu, const vh)
{
	decl dden, vdif, dmp, dsp, dpsi, dpsia;

	dden = m_dsp0^2*sumsqrc(vh) + dsigu^2;
	vdif = (m_my[][1] - dxi) .* vh;

	dmp = (m_dsp0^2 * sumc(vdif) + m_dmp0*dsigu^2) / dden;
	dsp = sqrt(m_dsp0^2 * dsigu^2 / dden);

	if(!m_flML)	//parameter estimation
		dpsi = dmp + dsp * rann(1,1);	// psi|. ~ N(mvp, svp^2)
	else	//ML estimation: generate xi for m_kk = 0,1,...,5+2*m_flt+m_flb
	{
		if(m_kk < 6+2*m_flt+2*m_flb)
			dpsi = dmp + dsp * rann(1,1);	// psi|. ~ N(mvp, svp^2)
		else
		{
			dpsi = dpsia = m_vposm[5+m_flt+m_flb];
			if(m_k>=0)
				m_mmls[m_k][6+2*m_flt+2*m_flb] = densn((dpsia-dmp)/dsp) / dsp;
		}		
	}
	
	return dpsi;
}

/*
** SampSigu samples parameter sigma_u
*/
RSVGH::SampSigu(const dxi, const dpsi, const vh)
{
	decl dn, dd, dsigu, dsigua;

	dn = 0.5*m_ns + m_dsun0;
	dd = m_dsuS0 + 0.5 * sumsqrc(m_my[][1]-dxi-dpsi*vh);

	dsigu = sqrt(1 / rangamma(1, 1, dn, dd));

	if(!m_flML)	//parameter estimation
		dsigu = sqrt(1 / rangamma(1, 1, dn, dd));
	else	//ML estimation: generate sigma_u for m_kk = 0,1,...,6+2*m_flt+m_flb
	{
		if(m_kk < 6+2*m_flt+2*m_flb+m_flbc)
			dsigu = sqrt(1 / rangamma(1, 1, dn, dd));
		else
		{
			dsigu = dsigua = m_vposm[5+m_flt+m_flb+m_flbc];
			if(m_k>=0)
				m_mmls[m_k][6+2*m_flt+2*m_flb+m_flbc] = densgamma(1/dsigua^2, dn, dd) / dsigua^4;
		}		
	}

	return dsigu;
}

/*
** SingleMove samples log-volatilities h_1, h_2, ..., h_T separately
*/
RSVGH::SingleMove(const dmu, const dphi, const dsigeta, const drho, const dbeta, const dnu, const dxi, const dsigu,
					const vho, const vz, const my)
{
	decl dmu_z, dsig_z, dbz, vzb, vxb, vmuh, vhn, vyn, vln, vyo, vlo, vac, dsigh1, dsigh2T1, dsighT, t, T;

	if(m_flt){
		dmu_z = dnu / (dnu - 2);
		dsig_z = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	}
	else{
		dmu_z = 1;
		dsig_z = 0;
	}
	dbz = sqrt(dbeta^2 * dsig_z^2 + dmu_z);
	vzb = vz - dmu_z;
	vxb = my[][1] - dxi;

	/*--- Single-move ---*/
	vhn = vho;
	vmuh = vyn = vyo = vln = vlo = vac = zeros(m_ns, 1);

	// t = 1
	dsigh1 = sqrt( 1 / ( (1 - dphi^2 + dphi^2/(1-drho^2))/dsigeta^2 + 1/dsigu^2 ) );
	vmuh[0]	= dsigh1^2 * (-0.5 + dphi*(vho[1]-(1-dphi)*dmu) / (dsigeta^2*(1-drho^2)) + (1-dphi^2)*dmu/dsigeta^2 + vxb[0] /dsigu^2);
	vhn[0] = vmuh[0] + dsigh1 * rann(1,1);
	vyn[0] = (dbz * my[0][0] * exp(-vhn[0]/2) - dbeta * vzb[0]) / sqrt(vz[0]);	
	vln[0] = -0.5 * vyn[0]^2 / (1 - drho^2) + drho * vyn[0] * (vhn[1]-dmu - dphi*(vhn[0]-dmu)) / ((1-drho^2) * dsigeta);
	vyo[0] = (dbz * my[0][0] * exp(-vho[0]/2) - dbeta * vzb[0]) / sqrt(vz[0]);	
	vlo[0] = -0.5 * vyo[0]^2 / (1 - drho^2) + drho * vyo[0] * (vho[1]-dmu - dphi*(vho[0]-dmu)) / ((1-drho^2) * dsigeta);
	vac[0] = ranu(1,1) < exp(vln[0] - vlo[0]);
	vhn[0] = vac[0] ? vhn[0] : vho[0];

	// t = 2,...,T-1
	dsigh2T1 = sqrt( 1 / ( (1+dphi^2)/(1-drho^2) / dsigeta^2 + 1/dsigu^2 ) );
	for(t=1; t<=m_ns-2; t++){
		vmuh[t] = dsigh2T1^2 * ( -0.5 + (dphi*(vhn[t+1]+vhn[t-1])+(1-dphi)^2*dmu)/(dsigeta^2*(1-drho^2)) + vxb[t]/dsigu^2
									+ drho * dsigeta * vyn[t-1] / ((1-drho^2)*dsigeta^2) );
		vhn[t] = vmuh[t] + dsigh2T1 * rann(1,1);
		vyn[t] = (dbz * my[t][0] * exp(-vhn[t]/2) - dbeta * vzb[t]) / sqrt(vz[t]);	
		vln[t] = -0.5 * vyn[t]^2 / (1 - drho^2) + drho * vyn[t] * (vhn[t+1]-dmu - dphi*(vhn[t]-dmu)) / ((1-drho^2) * dsigeta);
		vyo[t] = (dbz * my[t][0] * exp(-vho[t]/2) - dbeta * vzb[t]) / sqrt(vz[t]);	
		vlo[t] = -0.5 * vyo[t]^2 / (1 - drho^2) + drho * vyo[t] * (vho[t+1]-dmu - dphi*(vho[t]-dmu)) / ((1-drho^2) * dsigeta);
		vac[t] = ranu(1,1) < exp(vln[t] - vlo[t]);
		vhn[t] = vac[t] ? vhn[t] : vho[t];
	}

	// t = T
	dsighT = sqrt( 1 / ( 1 / (dsigeta^2*(1-drho^2)) + 1/dsigu^2 ) );
	T = m_ns - 1;
	vmuh[T] = dsighT^2 * ( -0.5 + ((1-dphi)*dmu + dphi*vhn[T-1] + drho*dsigeta*vyn[T-1])/(dsigeta^2*(1-drho^2)) + vxb[T]/dsigu^2 );
	vhn[T] = vmuh[T] + dsighT * rann(1,1);
	vyn[T] = (dbz * my[T][0] * exp(-vhn[T]/2) - dbeta * vzb[T]) / sqrt(vz[T]);
	vln[T] = -0.5 * vyn[T]^2;
	vyo[T] = (dbz * my[T][0] * exp(-vho[T]/2) - dbeta * vzb[T]) / sqrt(vz[T]);	
	vlo[T] = -0.5 * vyo[T]^2;
 	vac[T] = ranu(1,1) < exp(vln[T] - vlo[T]);
	vhn[T] = vac[T] ? vhn[T] : vho[T];
	
	if(!m_flML && m_k>=0) m_cmh += sumc(vac);

	return vhn;
}

/*
**  MultiMove implements multi-move sampler
**  by Omori and Watanabe(2008)
*/
RSVGH::MultiMove(const dphi, const dsigeta, const drho, const dgam, const dbeta, const dnu, const dxb, const dpsi, const dsigu,
			   const vh, const vz, const vzb)
{
	decl dmu_z, nite, vk, ir, id, ird, fli, vh_s, vetao, vetaoi;
	decl myi, vho, vhn, dhrd1, vhh, i, j, t;
	decl vmu, vsig2, dlfhrd1, vder, vder1, vdel, vA, vB;
	decl vzi, vzbi, vDh, vKh, vJh, vbh, vgn;
	decl vgh, vyh, dLh, mBh, mQh;
	decl vD, vZ, mG, vH, vH0, da, da0, dP, dH2, dV, du0;
	decl ve, vDinv, vK, vL, mJ, vu, dr, dU;
	decl deps, dC, dCinv;
	decl vmun, vsign2, dlfhrdn1, vmuo, vsigo2, dlfhrdo1;
	decl vetan, dposn, dposo, dpron, dproo, dfrac;
	decl car;

	if(m_flt) dmu_z = dnu / (dnu - 2);
	else dmu_z = 1;
	
	nite = 5;	//# of iteration for h_hat
			
	vh_s = vh;	//sample box

	vk = 0	  						//stochastic knots
	   ~ floor(m_ns * (range(1, m_nKm)
	   				   +ranu(1, m_nKm))/(m_nKm+2))
	   ~ m_ns-1;

	vetao = ((vh[1:] - dphi*vh[:m_ns-2]) / dsigeta) | 0;  //eta_old
	
	for(i=0 ; i<m_nKm+1 ; i++)
	{
		ir  = vk[i];
		id  = vk[i+1] - vk[i]; 
		ird = vk[i+1];

		if(id<2){if(!m_flML && m_k>=0) m_cah++;}
		else
		{
			fli = i < m_nKm;	// fli = 1 if it is not the last block		 
		
			myi = m_my[ir:ird][];	//r:r+d, (id+1)*2 vector
			vzi  = vz[ir:ird];
			vzbi = vzb[ir:ird];

			if(i)
			{
				vetaoi = vetao[ir:ird];
		
				vho			//current(old) point
				 = vhn		//new point, r:r+d
				 = vh_s[ir] | zeros(id, 1);
		
				for(t=0 ; t<id ; t++)
					vho[t+1] = dphi * vho[t] + vetaoi[t]*dsigeta;  //update r+1:r+d	

				if(fli) dhrd1 = dphi * vho[id] + vetaoi[id]*dsigeta;
			}
			else	//initial case
			{
				ir = -1;
				id++;
				myi = (0 ~ dxb) | myi;	//add y_0 = (y_{1,0}, y_{2,0}) = (0, c) = unconditional means (c = xi + psi*mu)
				vzi =  dmu_z | vzi;
				vzbi = 0 | vzbi;
				vho = 0 | vh_s[:ird];	//add h_0
				vhn = zeros(id+1, 1);
				dhrd1 = vh_s[ird+1];
				vetaoi = vh[0]*sqrt(1-dphi^2)/dsigeta | vetao[:ird];
			}

			/* finding mode	& draw candidate */
			for(j=0 ; j<=nite ; j++)
			{
				vhh = j ? vhn : vho;		//h_hat, r:r+d
			   
				vmu = myi[:id-1][0] - dbeta * vzbi[:id-1] .* exp(vhh[:id-1]/2) * dgam	// y_{1,t} - mu_t (t<n)
					- drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam .* (vhh[1:] - dphi * vhh[:id-1]) / dsigeta;
				vsig2 = (1 - drho^2) * vzi .* exp(vhh) * dgam^2;	// sigma_t^2 (t<n)
			
				vder = -0.5 * dbeta * vzbi[:id-1] .* exp(vhh[:id-1]/2) * dgam	//dmu_t/dh_t (t<n)
					 - drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam .* ((vhh[1:] - dphi * vhh[:id-1])/2 - dphi) / dsigeta;

				vder1 = - drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam / dsigeta;
				vder1 = i == 0 ? (0|0|vder1[1:]) : (0|vder1);	//dmu_{t-1}/dh_t

				if(fli)	// if not the last block (s+m < n)
				{
					vmu |= myi[id][0] - dbeta * vzbi[id] * exp(vhh[id]/2) * dgam
						 - drho * sqrt(vzi[id])	* exp(vhh[id]/2) * dgam * (dhrd1 - dphi * vhh[id]) / dsigeta;
					vder |= -0.5 * dbeta * vzbi[id] * exp(vhh[id]/2) * dgam
						  - drho * sqrt(vzi[id]) * exp(vhh[id]/2) * dgam * ((dhrd1 - dphi * vhh[id])/2 - dphi) / dsigeta;
					dlfhrd1 = -0.5 * (dhrd1 - dphi * vhh[id])^2 / dsigeta^2;	// additional term for d_t at t = s+m < n (kappa in my note)
				}
				else{
					vmu |= myi[id][0] - dbeta * vzbi[id] * exp(vhh[id]/2) * dgam;
					vsig2[id] = vzi[id] * exp(vhh[id]) * dgam^2;
					vder |= -0.5 * dbeta * vzbi[id] * exp(vhh[id]/2) * dgam;
					dlfhrd1 = 0;
				}
				if(!i) vmu[0] = 0;
	
				//vmu, vsig2, vder, vder1 are t=(r:r+d),(id+1)*1 vector
			
				vdel = -0.5									//delta
					 + 0.5 * vmu[1:].^2 ./ vsig2[1:]
					 - vmu[1:] ./ vsig2[1:] .* vder[1:]
					 - vmu[:id-1] ./ vsig2[:id-1] .* vder1[1:]
					 + dpsi*(myi[1:][1] - dxb - dpsi*vhh[1:]) ./ (dsigu^2);	// additional term due y_2
				vA = 0.5 + vder[1:].^2 ./ vsig2[1:]
						 + vder1[1:].^2 ./ vsig2[:id-1]
						 + ones(id,1)*(dpsi^2)./(dsigu^2);	// additional term due to y_2
				vB = 0
				   | (vder[1:id-1] .* vder1[2:] ./ vsig2[1:id-1]);
			
				if(fli)	// additional term for d_t at t = s+m < n (kappa in my note)
				{
					vdel[id-1] += dphi * (dhrd1 - dphi * vhh[id]) / dsigeta^2;
					vA[id-1] += dphi^2 / dsigeta^2;
//					vB[id-1] += dphi / dsigeta^2; //omitted in Takahashi et al. (2016 IJF) 
				}
				  
				//vdel, vA, vB are t=(r+1:r+d), id*1 vector
 
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
				mG = (1 ~ vJh[1:]*dsigeta) ./ vKh;
				vH = 0 ~ dsigeta;
			
				vH0 = i ? (0 ~ dsigeta) : (0 ~ dsigeta/sqrt(1-dphi^2));
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
						dC = vH * (unit(2) - mG[t][] ' vDinv[t] * mG[t][]
								   - mJ[t][] ' dU * mJ[t][]) * vH';
						deps = 0;
						vu[t] = vH * (mG[t][] ' vDinv[t] * ve[t]
									  + mJ[t][] ' dr)
							  + deps;
						dV = vH * (mG[t][] ' vDinv[t] * vZ[t]
								   + mJ[t][] ' dU * vL[t]);
						dCinv = 1 / dC;
						dr = vZ[t] * vDinv[t] * ve[t] + vL[t] * dr
						   - dV * dCinv * deps;
						dU = vZ[t]^2 * vDinv[t] + vL[t]^2 * dU
						   + dV^2 * dCinv;
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
							dC = vH * (unit(2) - mG[t][] ' vDinv[t] * mG[t][]
									   - mJ[t][] ' dU * mJ[t][]) * vH';
							deps = sqrt(dC) * rann(1, 1); 
							vu[t] = vH * (mG[t][] ' vDinv[t] * ve[t]
										  + mJ[t][] ' dr)
								  + deps;
							dV = vH * (mG[t][] ' vDinv[t] * vZ[t]
									   + mJ[t][] ' dU * vL[t]);
							dCinv = 1 / dC;
							dr = vZ[t] * vDinv[t] * ve[t] + vL[t] * dr
							   - dV * dCinv * deps;
							dU = vZ[t]^2 * vDinv[t] + vL[t]^2 * dU
							   + dV^2 * dCinv;
					
						}
					
						dH2 = vH0 * vH0';
						dC = dH2 - dH2^2 * dU;
						du0 = dH2 * dr + sqrt(dC) * rann(1, 1);
					
						vhn[1] = da0 + du0;			//t=r+1
						for(t=0 ; t<id-1 ; t++)
							vhn[t+2] = dphi * vhn[t+1] + vu[t];	//propose new h
									 
						//AR step
						dLh = sumc(-0.5 * vhh - 0.5 * vmu.^2 ./ vsig2
							  - 0.5*(myi[][1] - dxb - dpsi*vhh).^2 ./ (dsigu^2)) + dlfhrd1;
					
						mBh = (0 ~ diag(vB[1:])) | 0;
						mQh = diag(vA) + mBh + mBh';
									
						dpron = dLh + vdel ' (vhn - vhh)[1:]
						 	  - 0.5 * (vhn - vhh)[1:] ' mQh * (vhn - vhh)[1:];
					
						vmun = myi[:id-1][0] - dbeta * vzbi[:id-1] .* exp(vhn[:id-1]/2) * dgam
							 - drho * sqrt(vzi[:id-1]) .* exp(vhn[:id-1]/2) * dgam .* (vhn[1:] - dphi * vhn[:id-1]) / dsigeta;
						vsign2 = (1 - drho^2) * vzi .* exp(vhn) * dgam^2;
				
						if(fli)
						{
							vmun |= myi[id][0] - dbeta * vzbi[id] * exp(vhn[id]/2) * dgam
								  - drho * sqrt(vzi[id]) * exp(vhn[id]/2) * dgam * (dhrd1 - dphi * vhn[id]) / dsigeta;
							dlfhrdn1 = -0.5 * (dhrd1 - dphi * vhn[id])^2 / dsigeta^2;
						}
						else
						{
							vmun |= myi[id][0] - dbeta * vzbi[id] * exp(vhn[id]/2) * dgam;
							vsign2[id] = vzi[id] * exp(vhn[id]) * dgam^2;
							dlfhrdn1 = 0;
						}
						if(!i) vmun[0] = 0;
				
						dposn = sumc(-0.5 * vhn	- 0.5 * vmun.^2 ./ vsign2
								- 0.5*(myi[][1] - dxb - dpsi*vhn).^2 ./ (dsigu^2)) + dlfhrdn1;
				
						if(!m_flML && m_k>=0) m_cah++;

						if(car>=100)
							print("h! ");

					}while(ranu(1, 1) >= exp(dposn - dpron) && car<100);

				}	//end draw candidate
  						   
			}	//end finding mode & draw candidate	

			//MH step
			dproo = dLh + vdel ' (vho - vhh)[1:]
			 	  - 0.5 * (vho - vhh)[1:] ' mQh * (vho - vhh)[1:];
									  
			vmuo = myi[:id-1][0] - dbeta * vzbi[:id-1] .* exp(vho[:id-1]/2) * dgam
				 - drho * sqrt(vzi[:id-1]) .* exp(vho[:id-1]/2) * dgam .* (vho[1:] - dphi * vho[:id-1]) / dsigeta;
			vsigo2 = (1 - drho^2) * vzi .* exp(vho) * dgam^2;
			if(fli)
			{
				vmuo |= myi[id][0] - dbeta * vzbi[id] * exp(vho[id]/2) * dgam
					  - drho * sqrt(vzi[id]) * exp(vho[id]/2) * dgam * (dhrd1 - dphi * vho[id]) / dsigeta;
				dlfhrdo1 = -0.5 * (dhrd1 - dphi * vho[id])^2 / dsigeta^2;
			}
			else{
				vmuo |= myi[id][0] - dbeta * vzbi[id] * exp(vho[id]/2) * dgam;
				vsigo2[id] = vzi[id] * exp(vho[id]) * dgam^2;
				dlfhrdo1 = 0;
			}
			if(!i) vmuo[0] = 0;
		 
			dposo = sumc(-0.5 * vho	- 0.5 * vmuo.^2 ./ vsigo2 - 0.5*(myi[][1] - dxb - dpsi*vho).^2 ./ (dsigu^2)) + dlfhrdo1;

			dfrac = exp(  dposn + min(dposo|dproo) - dposo - min(dposn|dpron) );
							
			if(ranu(1, 1) < dfrac){
				if(!m_flML && m_k>=0) m_cmh++;
			}					   
			else
				vhn[1:] = vho[1:];
	
			vh_s[ir+1:ird] = vhn[1:];
		
		}
	}	//end sampling all blocks		

	return vh_s;
}


/*
** ML computes marginal likelihood 
*/
RSVGH::ML(const nI, const nR, const nsim)
{
	decl msamp, vnusamp, vnu, vf, dnum;
	decl vlik, dpos, dposse, dpri, iBm, i;
	decl dvm, vfl, mlsh, mOm, mOmsum, vtaua;
	decl mml;
	decl dmu, dsigeta, drho, dphi, dbeta, dnu, dxi, dsigu, dpsi;
	decl times = timer();

	m_vposm = loadmat(sprint(m_sdirectory, "/posstat_",m_sfilename,".csv"))[][0];
		//posterior mean of phi, sigma_eta, rho, mu, beta (if skew t), nu (if t), xi, psi(if bc), sigma_u

	iBm = min(500, nsim/2);			//bandwidth

	/*--- likelihood part ---*/
	
	vlik = PFilter(nI, nR);

	
	/*--- posterior part ---*/

	m_flML = 1;
	MCMC(nsim, 0.1*nsim);		//reduced run
	m_flML = 0;

	dvm = meanc(m_mmls);
	vfl = <1; -1; 1; -1; 1>;	// phi, (sigma_eta,rho), mu
	if(m_flt)
	{
		if(m_flb)
			vfl |= <1; -1>;	// beta

		vfl |= <1; -1>;	// nu
	}
	vfl |= <1; 1>;	// xi, sigma_u
	if(m_flbc)
		vfl |= <1>;	// psi
	dpos = log(dvm) * vfl;
						 
	mlsh = m_mmls - dvm;
	mOmsum = zeros(7+2*m_flt+2*m_flb+m_flbc, 7+2*m_flt+2*m_flb+m_flbc);
	for(i=0 ; i<iBm ; i++){
		mOm = mlsh[i:][] ' mlsh[:nsim-i-1][] / nsim;
		mOmsum += (1 - i/(iBm+1)) * (mOm + mOm' * (i>0));
	}

	mOmsum /= nsim;
	vtaua = vfl ./ dvm';
	dposse = sqrt(vtaua ' mOmsum * vtaua);


	/*--- prior part ---*/
	dphi = m_vposm[0];
	dsigeta = m_vposm[1];
	drho = m_vposm[2];
	dmu = m_vposm[3];
	dxi = m_vposm[4+m_flt+m_flb];
	dsigu = m_vposm[5+m_flt+m_flb+m_flbc];

	dpri
	= -0.5*log(2*M_PI) - log(m_dsm0) - 0.5*(dmu-m_dmm0)^2/m_dsm0^2
	+ log(densbeta((dphi+1)/2, m_dap0, m_dbp0)) - log(2)
	+ log(densbeta((drho+1)/2, m_dar0, m_dbr0)) - log(2)
	+ 0.5*m_dn0*log(m_dS0/2) - loggamma(m_dn0/2) - (m_dn0+2)*log(dsigeta) - 0.5*m_dS0/dsigeta^2
	+ -0.5*log(2*M_PI) - log(m_dsx0) - 0.5*(dxi-m_dmx0)^2/m_dsx0^2
	+ 0.5*m_dsun0*log(m_dsuS0/2) - loggamma(m_dsun0/2) - (m_dsun0+2)*log(dsigu) - 0.5*m_dsuS0/dsigu^2;

	if(m_flt)
	{
		dnu = m_vposm[4+m_flb];
		dpri += log(densgamma(dnu, m_dnn0, m_dnS0));

		if(m_flb){
			dbeta = m_vposm[4];
			dpri += -0.5*log(2*M_PI) - log(sqrt(m_dsb02)) - 0.5*(dbeta - m_db0)^2/m_dsb02;
		}
	}

	if(m_flbc)
	{
		dpsi = m_vposm[5+m_flt+m_flb];
		dpri += -0.5*log(2*M_PI) - log(m_dsp0) - 0.5*(dpsi - m_dmp0)^2/m_dsp0^2;
	}

//	dpri = log(0.5*densbeta((m_vposm[0]+1)/2, m_dap0, m_dbp0))
//		 + log(1/m_vposm[1]^4 * densgamma(1/m_vposm[1]^2, m_dn0, m_dS0))
//		 + log(0.5*densbeta((m_vposm[2]+1)/2, m_dar0, m_dbr0))
//		 + log(densn((m_vposm[3]-m_dmm0)/m_dsm0) / m_dsm0);	// phi, sigma_eta, rho, mu
//	if(m_flt)
//	{
//		if(m_flb)
//			dpri +=	log(densn((m_vposm[4]-m_db0)/sqrt(m_dsb02)) / sqrt(m_dsb02));	// beta
//
//		dpri += log(densgamma(m_vposm[4+m_flb], m_dnn0, m_dnS0) / (1 - probgamma(4, m_dnn0, m_dnS0)));	// nu
//	}
//	dpri +=	log(densn((m_vposm[4+m_flt+m_flb]-m_dmx0)/m_dsx0) / m_dsx0);	// xi
//	if(m_flbc)
//		dpri += log(densn((m_vposm[5+m_flt+m_flb]-m_dmp0)/m_dsp0) / m_dsp0);	// psi
//	dpri += log(1/m_vposm[5+m_flt+m_flb+m_flbc]^4 * densgamma(1/m_vposm[5+m_flt+m_flb+m_flbc]^2, m_dsun0, m_dsuS0));
//		// sigma_u


	/*--- output ---*/

	print("\n\n========================================"
		  "\n     Marginal likelihood estimation"
		  "\n========================================"
		  "\n                                rsvgh.ox"
		  "\nModel: RSV");
	if(m_flbc) print("C");
	if(m_flt)
	{
		if(m_flb) print("SK");
		print("t model");
	}
	else
		print("n model");
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
RSVGH::PFilter(const nI, const nR)
{
	decl dphi, dsigeta, drho, dmu, dbeta, dnu, dxi, dpsi, dsigu, dmuz, dsigz, db, vzb;
	decl vf, vlik, vz, vh, vw, vfhat, vmu, vgj, vuc, vuci;
	decl vzi, vind, cuci, cuc, vhi, vfhati, vgji, i, j, t;

	dphi     = m_vposm[0];
	dsigeta  = m_vposm[1];
	drho     = m_vposm[2];
	dmu      = m_vposm[3];
	dbeta    = m_flb ? m_vposm[4] : 0;
	dnu      = m_flt ? m_vposm[4+m_flb] : M_INF_POS;
	dxi      = m_vposm[4+m_flb+m_flt];
	dpsi     = m_flbc ? m_vposm[5+m_flb+m_flt] : 1;
	dsigu    = m_vposm[5+m_flb+m_flt+m_flbc];

	if(m_flt)
	{
		dmuz = dnu / (dnu - 2);
		dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
		db = dbeta^2*dsigz^2 + dmuz;
	}
	else
	{
		dmuz = 1;
		dsigz = 0;
		db = 1;
		vz = ones(nI, 1);
	}
	
	vf = zeros(m_ns, 1);
	vlik = zeros(nR, 1);
	
	if(nR > 1)
		print("\n[Particle filter] Iteration:");

	/*---------- Iteration START -----------*/
	for(j=0 ; j<nR ; j++)
	{
		if(nR > 1)
			print(j);
	
		if(m_flt)
			vz = 1 ./ rangamma(nI, 1, dnu/2, dnu/2);

		vh = dmu + rann(nI, 1) * dsigeta / sqrt(1 - dphi^2); // generate h_1
	
		vw = sqrt(db)*densn((sqrt(db)*m_my[0][0] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2))) ./ (sqrt(vz) .* exp(vh/2))
				.* densn((m_my[0][1] - dxi - dpsi*vh)./dsigu) ./ dsigu;
	
		vf[0] = meanc(vw);
		vfhat = vw / sumc(vw);
		
		/*---------- Filtering START -----------*/
		for(t=0 ; t<m_ns-1 ; t++)
		{
			vmu = dmu + dphi * (vh - dmu) + drho * dsigeta * (sqrt(db)*m_my[t][0] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2));
			vgj = sqrt(db)*densn(sqrt(db)*m_my[t+1][0]./(sqrt(dmuz)*exp(vmu/2))) ./ (sqrt(dmuz)*exp(vmu/2))
					.* densn((m_my[t+1][1] - dxi - dpsi*vmu)./dsigu) ./ dsigu
					.* vfhat;
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
		
			vzi = vz[vind];

			vmu = dmu + dphi * (vhi - dmu) + drho * dsigeta * (sqrt(db)*m_my[t][0] - dbeta*(vzi-dmuz).*exp(vhi/2)) ./ (sqrt(vzi) .* exp(vhi/2));
			vh = vmu + rann(nI, 1) * dsigeta * sqrt(1-drho^2);
			
			if(m_flt)
				vz = 1 ./ rangamma(nI, 1, dnu/2, dnu/2);
				
			vw = sqrt(db)*densn((sqrt(db)*m_my[t+1][0] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2))) ./ (sqrt(vz) .* exp(vh/2))
					.* densn((m_my[t+1][1] - dxi - dpsi*vh)./dsigu) ./ dsigu
					.* vfhati ./ vgji;
		
			vf[t+1] = meanc(vw);
			vfhat = vw / sumc(vw);
			
		}
		/*---------- Filtering END -----------*/
	  		   
		vlik[j] = sumc(log(vf));
	}
	/*----------- Iteration END ------------*/

  
	return vlik;
}

/*
**  fTsvar & fGeweke return the diagnostics for convergence
*/
RSVGH::fTsvar(const vx, const iBm)
{
	decl ns, vsp;

	ns = rows(vx);
	vsp = periodogram(vx, iBm, 2, 1) / ns;
	
	return M_2PI * vsp[0];
}

RSVGH::fGeweke(const vx, const iBm)
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

RSVGH::fBetaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(da+db));
	else if(cm==2) return(da*db/((da+db)^2*(da+db+1)));
}

RSVGH::fGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/db);
	else if(cm==2) return(da/db^2);
}

RSVGH::fInvGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(db-1));
	else if(cm==2) return(db^2/(sqr(da-1)*(da-2)));
}
