/*----------------------------------------------------------
                        SVGH.ox
----------------------------------------------------------*/
											 
/*
**  CLASS FILE
**
**  MCMC estimation for stochastic volatility model
**  with generalized hyperbolic distribution
*/

#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<SVGH.h>

SVGH::SVGH(const vdata)	//vdata=0 means simulation study
{												   
	m_iRun = ismatrix(vdata) ? 1 : 0;
						//0:simulation, 1:empirical
	m_vy = m_vyorg = vdata;		  

	if(m_iRun){
		m_ns = rows(vdata);
		m_nKm = int(m_ns/25)-1; //# of blocks in multi-move sampler
	}
						
	m_flt = 1;			//model specification: t
	m_flb = 1;			//model specification: Skew
	m_vini = 0;			//initial value
	m_ncount = 1000;	//counter of iteration
	m_flML = 0;			//flag for ML estimation
	m_flf  = 0;			//flag for VaR and ES backtest
}

/*
** SetFilename set main name for output files
*/
SVGH::SetFilename(const sfilename)
{
	m_sfilename = sfilename;
}

/*
** SetDirectory set directory for output files
*/
SVGH::SetDirectory(const sdirectory)
{
	m_sdirectory = sdirectory;
}

/*
** SetModel set model specification
**   iModel:  1(SVSKt) / 0(SVt) / -1(SVn)
*/
SVGH::SetDist(const iDist)
{
	m_flt = (iDist>=0);
	m_flb = (iDist==1);
}

/*
** SetRanseed set ranseed for random generater
*/
SVGH::SetRanseed(const irs)
{
	ranseed(irs);
}


/*
** SetBlocks set # of blocks in multi-move sampler
*/
SVGH::SetBlocks(const nKm)
{
	m_nKm = nKm-1;
}

/*
** SetCount set counter of iteration
*/
SVGH::SetCount(const ncount)
{
	m_ncount = ncount;
}

/*
** SetForecast set volatility, VaR and ES forecasts
**   iForecast:  1(on) / 0(off)
*/
SVGH::SetForecast(const iForecast)
{
	m_flf = iForecast;
}

/*
** SetPar set true parameter for simulation
*/
SVGH::SimObs(const ns, const dphi, const dsig, const drho, const dmu, const dbeta_a, const dnu_a)
{
	decl mSigsq, me, vy, vh, vz, dmuz, dsigz, vzb, db, dbeta, dnu, i;
	
	m_ns = ns;

	mSigsq
	  = choleski( (1 ~ drho*dsig) | (drho*dsig ~ dsig^2));

	vy = vh = zeros(ns, 1);
	me = mSigsq * rann(2,ns);
	
	vh[0] = dmu + rann(1, 1) * dsig / sqrt(1-dphi^2);
	for(i=0 ; i<ns-1; i++)
	{
		vh[i+1] = dmu + dphi * (vh[i] - dmu) + me[1][i];
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
	db = dbeta^2*dsigz^2+dmuz;

	vy = (dbeta*vzb + sqrt(vz).*me[0][]')/sqrt(db) .* exp(vh/2);

	m_vy = m_vyorg = vy;
	m_vh_t = vh;
	m_vz_t = vz;
	
	print("\n>>> Data is simualted. [SV");
	if(m_flt)
	{
		if(m_flb) print("SK");
		print("t");
	}
	else
		print("n");
	print(" model]"
		  "\nTrue value:\n phi   = ", dphi,
		  "\n sigma = ", dsig, "\n rho   = ", drho,
		  "\n mu    = ", dmu, "\n beta  = ", dbeta,
		  "\n nu    = ", dnu, "\n mu_z  = ", dmuz,
		  "\n sig_z = ", dsigz, "\n b     = ", db);

	savemat(sprint(m_sdirectory, "/simdata_", m_sfilename, ".csv"), vy~vh~vz);
	savemat(sprint(m_sdirectory, "/simpar_", m_sfilename, ".csv"), dphi|dsig|drho|dmu|dbeta|dnu);

	m_nKm = int(ns/25);
}

/*
** SetIni set initial value for MCMC sampling
*/
SVGH::SetIni(const dphi, const dsig, const drho, const dmu, const dbeta, const dnu, const vh)
{
	m_vini = dphi ~ dsig ~ drho ~ dmu ~ dbeta ~ dnu;
	m_vhini = vh;
}

/*
** SetPri set hyper parameters for priors
*/
SVGH::SetPri(const dap0, const dbp0, const dn0, const dS0, const dar0, const dbr0, const dmm0, const dsm0,
	const db0, const dsb02, const dnn0, const dnS0, const flprint)
{
	m_vpri = dap0 ~ dbp0 ~ dn0 ~ dS0 ~ dar0 ~ dbr0 ~ dmm0 ~ dsm0 ~ db0 ~ dsb02 ~ dnn0 ~ dnS0;
	if(flprint)
	{
		print("\n");
		print("Prior: (phi+1)/2  ~ Beta(", dap0, ", ", dbp0, ")\n");
		print("       sigma^(-2) ~ Gamma(", dn0, ", ", dS0, ")\n");
		print("       (rho+1)/2  ~ Beta(", dar0, ", ", dbr0, ")\n");
		print("       mu         ~ N(", dmm0, ", ", dsm0^2, ")\n");
		if(m_flt)
		{
			if(m_flb) print("       beta       ~ N(", db0, ", ", dsb02, ")\n");
			print("       nu         ~ Gamma(", dnn0, ", ", dnS0, ")\n");
			print("%c", {"Mean", "Stdev"}, "%r", {"       phi", "       sigma^2", "       rho", "       nu"}, "%cf", {"%7.2f"},
				  fBetaMoment(dap0,dbp0,1)*2-1 ~ sqrt(fBetaMoment(dap0,dbp0,2))*2 |
				  fInvGammaMoment(dn0,1/dS0,1) ~ sqrt(fInvGammaMoment(dn0,1/dS0,2)) |
				  fBetaMoment(dar0,dbr0,1)*2-1 ~ sqrt(fBetaMoment(dar0,dbr0,2))*2 |
				  fGammaMoment(dnn0,dnS0,1) ~ sqrt(fGammaMoment(dnn0,dnS0,2)));
		}
		else
		{
			print("%c", {"Mean", "Stdev"}, "%r", {"       phi", "       sigma^2", "       rho"}, "%cf", {"%7.2f"},
				  fBetaMoment(dap0,dbp0,1)*2-1 ~ sqrt(fBetaMoment(dap0,dbp0,2))*2 |
				  fInvGammaMoment(dn0,1/dS0,1) ~ sqrt(fInvGammaMoment(dn0,1/dS0,2)) |
				  fBetaMoment(dar0,dbr0,1)*2-1 ~ sqrt(fBetaMoment(dar0,dbr0,2))*2 );
		}
	}
}

/*
** DrawData draw graph of data (and state if m_iRun=0)
*/
SVGH::DrawData()
{
	Draw(0, m_vyorg');
	DrawTitle(0, "Data");
	if(!m_iRun){
		Draw(1, m_vh_t');
		DrawTitle(1, "Simulated log volatility ($h_t$)");
	}
	DrawDensity(1+!m_iRun, m_vyorg', 1, 1, 0, 1, 1);
	DrawTitle(1+!m_iRun, "Density of Data");	
	DrawTitle(2+!m_iRun, "QQ-plot to normal density");	
  	SaveDrawWindow(sprint(m_sdirectory, "/datadraw_",m_sfilename,".eps"));
  	CloseDrawWindow();
}

/*
** MCMC computes MCMC estimation
*/
SVGH::MCMC(const nsim, const nburn)
{
	decl dphi, dsig, drho, dmu, dbeta, dnu, i, ncyc;
	decl vh, vz, dmuz, dsigz, db, vzb, vyb0, npmt, msamp, vsamph, msamph;
	decl vsamp, iBm, iacf, aspar, aspar2, mpost, mposh, ascf, vaccptprob;
	decl vsampy, msamph1, dh1, dz1, vprob, vvar, mfly, ves, vvolf;
	decl vhb, vmu, vloglik, vloglikpf, vlogpri;
	decl time = timer();
	
	/*--- set initial ---*/
	if(!m_flML){
	  if(ismatrix(m_vini)){
		dphi  = m_vini[0];
		dsig  = m_vini[1];
		drho  = m_vini[2];
		dmu   = m_vini[3];
		dbeta = m_vini[4];
		dnu   = m_vini[5];
	  }
	  else{
		dphi  = 0.97;
		dsig  = 0.2;
		drho  = -0.3;
		dmu   = -10;
		dbeta = -0.3;
		dnu   = 15;
	  }
	}
	else{
		dphi  = m_vposm[0];
		dsig  = m_vposm[1];
		drho  = m_vposm[2];
		dmu   = m_vposm[3];
		dbeta = m_flb ? m_vposm[4] : 0;
		dnu   = m_flt ? m_vposm[4+m_flb] : M_INF_POS;
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
	}
		
	if(!m_flb) dbeta = 0;
	if(!m_flt) dnu = M_INF_POS;
								
	if(!m_flML)
	{
		if (ismatrix(m_vhini)) vh = m_vhini;
		else vh = dmu + rann(m_ns, 1) * dsig;
	}
	else{
		vh = loadmat(sprint(m_sdirectory, "/posh_",m_sfilename,".csv"))[][0];
		println("\n Posterior mean of state");
		println(vh[0:10]);
	}

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

	if(!m_iRun){
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
	npmt = 4 + m_flt + m_flb;			//# of parameters: phi, sigma, rho, mu, beta (if SKt), nu (if not normal)

	
  if(!m_flML){	/*--------- parameter estimation --------*/
  
	msamp = zeros(nsim, npmt);	//sample box
	vsamph = zeros(m_ns, 1);
	msamph = zeros(nsim, m_ns);
	vloglik = zeros(nsim, 1);
	vloglikpf = zeros(nsim, 1);
	vlogpri = zeros(nsim, 1);
	vsampy = zeros(nsim, 1);
	msamph1 = zeros(nsim, 2);
	m_cmp = m_cms = m_cmn = m_cmz = m_cmb = m_cmh = m_cah = 0; //acceptance counter
	m_cnu2 = 0; //counter of breaks in loop for sampling nu
							
	/*--- MCMC sampling ---*/
	if(m_ncount > 0) println("\n\nIteration:");

	/*----------- S A M P L I N G   S T A R T -----------*/
	for(m_k=-nburn ; m_k<nsim ; m_k++)
	{
		vyb0 = (sqrt(db)*m_vy.*exp(-vh/2) - dbeta*vzb) ./ sqrt(vz);

		/*--- sampling phi ---*/
		dphi = SampPhi(drho*dsig*vyb0, dphi, dsig, drho, vh-dmu);

		/*--- sampling sigma and rho ---*/
		[dsig, drho] = SampSigRho(dphi, dsig, drho, vh-dmu, vyb0);
  
		/*--- sampling mu ---*/
		dmu = SampMu(dphi, dsig, drho, dbeta, drho*dsig*vyb0, vh);

		if(m_flt)	// if t
		{
			/*--- sampling beta ---*/
			if(m_flb)
				dbeta = SampBeta(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz, vzb);

			/*--- sampling nu ---*/
			dnu = SampNu(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz);
			dmuz = dnu / (dnu - 2);
			dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
			db = dbeta^2*dsigz^2 + dmuz;
	
			/*--- sampling z ---*/
			vz = SampZ(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz);
			vzb = vz - dmuz;
		}

		/*--- sampling h ---*/
		vh = MultiMove(dphi, dsig, drho, exp(dmu/2)/sqrt(db), dbeta, dnu, vh-dmu, vz, vzb) + dmu;

		if(m_k >= 0)
		{
			/*--- storing sample ---*/
			if(m_flt)
				vsamp = m_flb
					  ? dphi ~ dsig ~ drho ~ dmu ~ dbeta ~ dnu
					  : dphi ~ dsig ~ drho ~ dmu ~ dnu;
			else
				vsamp = dphi ~ dsig ~ drho ~ dmu;
	
			msamp[m_k][] = vsamp;
			vsamph += vh;
			msamph[m_k][] = vh';

			if(m_flf)
			{
				/*--- VaR and Expected Shortfall ---*/
				if(m_flt)
					dz1 = 1 / rangamma(1, 1, dnu/2, dnu/2);
				else
					dz1 = 1;

				dh1 = dmu + dphi * (vh[m_ns-1] - dmu) + drho * dsig * (sqrt(db)*m_vy[m_ns-1] * exp(-vh[m_ns-1]/2)
					- dbeta * (vz[m_ns-1] - dmuz)) / sqrt(vz[m_ns-1])  // 2026-09-16 fix: was y * (-beta(z-muz)e^{h/2}) e^{h/2}; now matches RSVGH.ox
					+ sqrt(1-drho^2) * dsig * rann(1,1);
				vsampy[m_k] = (dbeta * (dz1-dmuz) + sqrt(dz1) .* rann(1,1)) / sqrt(db) * exp(dh1/2);
				msamph1[m_k][] = dh1 ~ exp(dh1);	// predicted log-voalitlity and volatility
			}
		}

		//print counter
		if(!imod(m_k, m_ncount) && m_ncount>0) println(m_k);
	}		 
    /*------------ S A M P L I N G   E N D --------------*/

	
	/*--- output result ---*/
	iBm = iacf = min(1000, nsim/2);	//bandwidth

	aspar =	{"phi", "sigma", "rho", "mu"};
	aspar2 = {"$\\phi$", "$\\sigma$", "$\\rho$", "$\\mu$"};
	if(m_flt)
	{
		if(m_flb){
		aspar ~= {"beta"};
		aspar2 ~= {"$\\beta$"};
		}
		aspar ~= {"nu"};
		aspar2 ~= {"$\\nu$"};
	}
	aspar ~= {"h(500)"};
	aspar2 ~= {"$h_{500}$"};
	
	print("\n                       ",
		  "[ESTIMATION RESULT]"
		  "\n--------------------------------"
		  "----------------------------------"
		  "\nParameter    Mean     Stdev     "
		  " 95%L      95%U   Geweke     Inef."  
		  "\n--------------------------------"
		  "----------------------------------");

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

	print("---------------------------------"
		  "------------------------ SVGH.ox\n");

	print("Model: SV");
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
		print("\n nu:     ", "%5.1f", m_cmn/nsim*100, "%",
			  "\n z:      ", "%5.1f", double(m_cmz/(nsim*m_ns)*100), "%");
		vaccptprob ~= m_cmn/nsim ~ m_cmz/(nsim*m_ns);
	}
	if(m_flb)
	{
		print("\n beta:   ", "%5.1f", double(m_cmb/nsim*100), "%");
	}
	print("\nAcceptance rate for multi-move sampler",
		  "\n AR step: ", "%5.1f",	  
		  (m_nKm+1)*nsim / m_cah * 100, "%",
		  "\n MH step: ", "%5.1f",
		  m_cmh / ((m_nKm+1)*nsim) * 100, "%");
	vaccptprob ~= (m_nKm+1)*nsim / m_cah ~ m_cmh / ((m_nKm+1)*nsim);
	
	if(m_flf){
		vprob = <.005, .01, .05, .1>;
		vvar = quantilec(vsampy, vprob)';
		mfly = vsampy .< vvar;
		ves = vsampy ' mfly ./ sumc(mfly); 
		vvolf = meanc(msamph1);
		print("\n\n1-day-ahead VaR and ES");
		print("%r", {"tail:", " VaR", " ES"}, 
			  "%cf", {"%10.5f"}, vprob | vvar | ves);
		print("\n\n1-day-ahead volatility forecast");
		print("%c", {"h", "exp(h)"}, 
			  "%cf", {"%10.5f"}, vvolf);
		savemat(sprint(m_sdirectory, "/var_",m_sfilename,".csv"), vvar);
		savemat(sprint(m_sdirectory, "/es_",m_sfilename,".csv"), ves);
		savemat(sprint(m_sdirectory, "/volf_",m_sfilename,".csv"), vvolf);
		savemat(sprint(m_sdirectory, "/h1_y_", m_sfilename, ".csv"), msamph1[][0] ~ vsampy); //1-day ahead log-volatility and return
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
		DrawMatrix(1, (m_vh_t ~ (vsamph/nsim))',
				   {"True", "PosMean"}, 1, 1);
		DrawTitle(1, "True and posterior mean of state");
		}
	    SaveDrawWindow(sprint(m_sdirectory, "/state_",m_sfilename,".eps"));
	    CloseDrawWindow();
	}

	savemat(sprint(m_sdirectory, "/posstat_", m_sfilename, ".csv"), mpost);
	savemat(sprint(m_sdirectory, "/posh_", m_sfilename, ".csv"), mposh);
	savemat(sprint(m_sdirectory, "/posstath_", m_sfilename, ".csv"), mposh);
	savemat(sprint(m_sdirectory, "/aprob_", m_sfilename,".csv"), vaccptprob);
  }			//end of parameter estimation

  else{		/*------------------------- ML estimation ---*/
	
	/*--- set variables ---*/
	m_mmls = zeros(nsim, 5+2*m_flt+2*m_flb);	// numerator and denominator for phi, (sigma, rho), nu (if t), beta (if skew) with MH and density of mu
	ncyc = 5+2*m_flt+2*m_flb;	//# of cycle for reduced run
 
	/*--- MCMC sampling ---*/
	
	print("\n[Reduced MCMC]    Iteration:");
							  
	/*----------- S A M P L I N G   S T A R T -----------*/
	for(m_kk=0 ; m_kk<ncyc ; m_kk++)
	{
		for(m_k=-nburn ; m_k<nsim ; m_k++)
		{
			vyb0 = (sqrt(db)*m_vy - dbeta * vzb .* exp(vh/2)) .* exp(-vh/2) ./ sqrt(vz);
	
			/*--- sampling phi ---*/
			if(m_kk < 2)
				dphi = SampPhi(drho*dsig*vyb0, dphi, dsig, drho, vh-dmu);
	
			/*--- sampling sigma and rho ---*/
			if(m_kk < 4)
				[dsig, drho] = SampSigRho(dphi, dsig, drho, vh-dmu, vyb0);
	
			/*--- sampling mu ---*/
			if(m_kk < 5)
				dmu = SampMu(dphi, dsig, drho, dbeta, drho*dsig*vyb0, vh);

			if(m_flt)
			{
				/*--- sampling beta ---*/
				if(m_flb && m_kk < 7)
					dbeta = SampBeta(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz, vzb);
			
				/*--- sampling nu ---*/
				if(m_kk < 7+2*m_flb)
				{
					dnu = SampNu(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz);
					dmuz = dnu / (dnu - 2);
					dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
					db = dbeta^2*dsigz^2 + dmuz;
				}
			
				/*--- sampling z ---*/
				vz  = SampZ(dphi, dsig, drho, dmu, dbeta, dnu, vh, vz);
				vzb = vz - dnu / (dnu - 2);
			}
			else
			{
				dbeta = 0;
				dnu = M_INF_POS;
				vz = ones(m_ns, 1);
				dmuz = 1;
				dsigz = 0;
				vzb = vz - 1;
				db = 1;
			}

			/*--- sampling h ---*/
			vh = MultiMove(dphi, dsig, drho, exp(dmu/2)/sqrt(db), dbeta, dnu, vh-dmu, vz, vzb) + dmu;
		}		 
	    /*------------ S A M P L I N G   E N D --------------*/

		/*--- fix parameter ---*/
		if(m_kk == 0)
			dphi = m_vposm[0];
		else if(m_kk == 2)
		{
			dsig = m_vposm[1];
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
		
		print(m_kk);
	}	/*------------------------------ end of cycle ---*/
	
  }		//end of ML estimation
  
	if(m_ncount>0 && !m_flML) println("\n\nTime: ", timespan(time));
}


/*
** SampPhi samples parameter phi by MH algorithm 
*/
SVGH::SampPhi(const vyb, const dphio, const dsig, const drho, const vhb)
{
	decl vmu, dsum, dmup, dsigp, dphin, dfrac, dphia, dqa, nite;

	/*--- draw candidate ---*/
	dsum = drho^2 * vhb[0]^2 + sumsqrc(vhb[1:m_ns-2]);

	dmup = sumc((vhb[1:]-vyb[:m_ns-2]).*vhb[:m_ns-2])/dsum;
	dsigp = sqrt(dsig^2*(1-drho^2) / dsum);

	nite = 0;
	do{
		nite++;
	  dphin = dmup + rann(1, 1) * dsigp;
		if(nite>10000){
			print("phi! ");
			dphin = dphio;
		}
	}while(fabs(dphin) >= 1 && nite<=10000);
	

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
SVGH::fPosSigRho(const vp, const dfunc, const vsc, const mhe)
{
	decl dsig, drho;

	dsig = exp(vp[0]);
	drho = 2 / (1+exp(-vp[1])) - 1;

	dfunc[0]
	  = log(1/dsig^4 * densgamma(1/dsig^2, m_dn0, m_dS0))
	  + log(densbeta((drho+1)/2, m_dar0, m_dbr0))
	  - m_ns * log(dsig)
	  - 0.5 * (m_ns - 1) * log(1 - drho^2)
	  - 0.5 * sumsqrc(m_vh[1:] - m_dphi*m_vh[:m_ns-2] - dsig*drho*m_vyb0[:m_ns-2]) / (dsig^2*(1-drho^2))
	  - 0.5 * (1 - m_dphi^2) * m_vh[0]^2 / dsig^2;
							   
	return 1;
}

/*
** SampSigRho samples parameter sigma and rho
** by MH algorithm 
*/
SVGH::SampSigRho(const dphi, const dsigo, const drhoo, const vhb, const vyb0)
{
	decl vparo, vparn, vparm, dlm, vsc, mhe, mS, vm;
	decl dlo, dln, dho, dhn, dsign, drhon, dfrac;
	decl dsiga, drhoa, vpara, dla, dha, dja, dqa;

	m_vh = vhb;
	m_dphi = dphi;
	m_vyb0 = vyb0;

	vparo = vparm = log(dsigo) | log((1+drhoo)/(1-drhoo));	//transformed

	MaxBFGS(fPosSigRho, &vparm, &dlm, 0, 1);//find mode
			
	fPosSigRho(vparm, &dlm, 0, 0);
  	Num1Derivative(fPosSigRho, vparm, &vsc);
  	Num2Derivative(fPosSigRho, vparm, &mhe);
	  
	mS = -invert(mhe);
	vm = vparm + mS * vsc;

	vparn = vm + choleski(mS) * rann(2, 1);	//proposal
				
	fPosSigRho(vparo, &dlo, 0, 0);
	fPosSigRho(vparn, &dln, 0, 0);
												    
	dsign = exp(vparn[0]);
	drhon = 2 / (1 + exp(-vparn[1])) - 1;
						  
	dho = dlm + (vparo-vparm)'vsc
		+ 0.5*(vparo-vparm)'mhe*(vparo-vparm)
		- log(dsigo)
		+ log(2) - log(1 - drhoo^2);
	dhn = dlm + (vparn-vparm)'vsc
		+ 0.5*(vparn-vparm)'mhe*(vparn-vparm)
		- log(dsign)
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
			fPosSigRho(vpara, &dla, 0, 0);
	
			dja = - log(dsiga) + log(2) - log(1 - drhoa^2);
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
** SampMu samples parameter mu
*/
SVGH::SampMu(const dphi, const dsig, const drho, const dbeta, const vyb, const vh)
{
	decl vdif, dsigm2, dmuh, dmu, dmua;

	vdif = vh[1:] - dphi*vh[:m_ns-2] - vyb[:m_ns-2];
	
	dsigm2 = 1/ (1 / m_dsm0^2
				 + ((1-drho^2)*(1-dphi^2)
					+ (m_ns-1)*(1-dphi)^2)
					/ (dsig^2*(1-drho^2)));
	dmuh = dsigm2
		 *(m_dmm0 / m_dsm0^2
		   + ((1-drho^2)*(1-dphi^2)*vh[0]
			  + (1-dphi)*sumc(vdif))/(dsig^2*(1-drho^2)));
	
	if(!m_flML)	//parameter estimation
		dmu = dmuh + sqrt(dsigm2) * rann(1, 1);
	else	//ML estimation: generate mu for m_k = 0,1,2,3,4
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
SVGH::fPosBeta(const dp, const dfunc, const dsc, const dhe)
{
	decl dbeta, dphi, dsig, drho, dmu, dnu, vh, vz, dmuz, dsigz, db, vzb, vhbd, vy1b;

	dbeta = dp;
	dphi = m_dphi;
	dsig = m_dsig;
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
	vy1b = (sqrt(db)*m_vy.*exp(-vh/2) - dbeta*vzb) ./ sqrt(vz);
	
	dfunc[0]
	  = - 0.5*(dbeta-m_db0)^2/m_dsb02
	  	+ 0.5*m_ns*log(db)
		- 0.5*sumsqrc(vy1b)
		- 0.5*sumsqrc(vhbd-drho*dsig*vy1b[:m_ns-2])/((1-drho^2)*dsig^2);
				   
	return 1;
}

/*
** SampBeta samples parameter beta
*/
SVGH::SampBeta(const dphi, const dsig, const drho, const dmu, const dbetao, const dnu, const vh, const vz, const vzb)
{
	decl dbetam, dlm, dsc, dhe, dS, dm, dbetan, dlo, dln, dho, dhn, dfrac, dbetaa, dla, dha, dqa;

	m_dphi = dphi;
	m_dsig = dsig;
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
	else	//ML estimation: generate beta for m_kk = 0,1,2,3,4,5 if skew
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
** fPosNu computes log posterior distribution of nu
*/
SVGH::fPosNu(const dp, const dfunc, const dsc, const dhe)
{
	decl dnu, dmuz, dsigz, db, vzb, vy1b;

	dnu = exp(dp);
	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = m_dbeta^2*dsigz^2+dmuz;

	vzb = m_vz - dmuz;
	vy1b = (sqrt(db)*m_vy.*exp(-m_vh/2) - m_dbeta*vzb) ./ sqrt(m_vz);
	
	dfunc[0]
	  = log(densgamma(dnu, m_dnn0, m_dnS0))
	  - 0.5 * sumsqrc(vy1b)
	  - 0.5 * sumsqrc(m_vhbd - m_dsig*m_drho*vy1b[:m_ns-2])/(m_dsig^2*(1-m_drho^2))
	  + 0.5 * m_ns * dnu * log(dnu / 2)
	  - m_ns * loggamma(dnu/2)
	  - 0.5 * dnu * sumc(log(m_vz) + 1 ./ m_vz)
	  + 0.5 * m_ns * log(db);
				   
	return 1;
}

/*
** SampNu samples parameter nu by MH algorithm 
*/
SVGH::SampNu(const dphi, const dsig, const drho, const dmu, const dbeta, const dnuo, const vh, const vz)
{
	decl dnuot, dnun, dnunt, dnum, dlm, dsc, dhe, dS, dm, nite;
	decl dfrac, dlo, dln, dho, dhn;
	decl dnua, dnuat, dla, dha, dja, dqa;

	m_dsig = dsig;
	m_drho = drho;
	m_dbeta = dbeta;
	m_vh = vh;
	m_vhbd = (vh[1:] - dmu) - dphi * (vh[:m_ns-2] - dmu);
	m_vz = vz;

	dnuot = dnum = log(dnuo);//transformed
				
	MaxBFGS(fPosNu, &dnum, &dlm, 0, 1); //find mode
			
	fPosNu(dnum, &dlm, 0, 0);
  	Num1Derivative(fPosNu, dnum, &dsc);
  	Num2Derivative(fPosNu, dnum, &dhe);
	  
	dS = -1 / dhe;
	dm = dnum + dS * dsc;

	nite = 0;
	do{
		nite++;
	  dnunt = dm + rann(1, 1) * sqrt(dS);	//proposal
		if(nite>100)
		{
			if(m_k>=0 && !m_flML) m_cnu2++;
			print("nu2! ");
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
		if(ranu(1,1) < dfrac)
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
				- log(1 - probn((log(4) - dm) / log(dS)));

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
SVGH::SampZ(const dphi, const dsig, const drho, const dmu, const dbeta, const dnu, const vh, const vzo)
{
	decl dmuz, dsigz, db, vzn, vzbo, vzbn, vy1bo, vy1bn, vlo, vln, vac, vhbd;

	dmuz = dnu / (dnu - 2);
	dsigz = sqrt(2*dnu^2 / ((dnu-2)^2*(dnu-4)));
	db = dbeta^2*dsigz^2+dmuz;
	vzn  = 1 ./ rangamma(m_ns, 1, (dnu+1)/2, dnu/2); //candidate

	vhbd = (vh[1:] - dmu) - dphi * (vh[:m_ns-2] - dmu);

	vzbo = vzo - dmuz;
	vzbn = vzn - dmuz;
	vy1bo = (sqrt(db)*m_vy.*exp(-vh/2) - dbeta*vzbo) ./ sqrt(vzo);
	vy1bn = (sqrt(db)*m_vy.*exp(-vh/2) - dbeta*vzbn) ./ sqrt(vzn);
						   
	vlo = -0.5 * vy1bo.^2
		- 0.5 * (((vhbd - dsig*drho *vy1bo[:m_ns-2]).^2) | 0) / (dsig^2 * (1 - drho^2));
	vln = -0.5 * vy1bn.^2
		- 0.5 * (((vhbd - dsig*drho*vy1bn[:m_ns-2]).^2) | 0) / (dsig^2 * (1 - drho^2));
 						   
	vac = ranu(m_ns, 1) .< exp(vln - vlo); 
	vzn = vac .? vzn .: vzo;

	if(!m_flML && m_k>=0) m_cmz += sumc(vac);
						
	return vzn;
}


/*
**  MultiMove implements multi-move sampler
**  by Omori and Watanabe(2008)
*/
SVGH::MultiMove(const dphi, const dsig, const drho, const dgam, const dbeta, const dnu, const vh, const vz, const vzb)
{
	decl dmu_z, nite, car, vk, ir, id, ird, fli, vh_s, vetao, vetaoi;
	decl vyi, vho, vhn, dhrd1, vhh, i, j, t;
	decl vmu, vsig2, dlfhrd1, vder, vder1, vdel, vA, vB;
	decl vzi, vzbi, vDh, vKh, vJh, vbh, vgn;
	decl vgh, vyh, dLh, mBh, mQh;
	decl vD, vZ, mG, vH, vH0, da, da0, dP, dH2, dV, du0;
	decl ve, vDinv, vK, vL, mJ, vu, dr, dU;
	decl deps, dC, dCinv;
	decl vmun, vsign2, dlfhrdn1, vmuo, vsigo2, dlfhrdo1;
	decl vetan, dposn, dposo, dpron, dproo, dfrac;
   
	dmu_z = dnu / (dnu - 2);

	nite = 5;	//# of iteration for h_hat
			
	vh_s = vh;	//sample box

	vk = 0	  						//stochastic knots
	   ~ floor(m_ns * (range(1, m_nKm) + ranu(1, m_nKm))/(m_nKm+2))
	   ~ m_ns-1;

	vetao = ((vh[1:] - dphi*vh[:m_ns-2]) / dsig) | 0;	//eta_old

	for(i=0 ; i<m_nKm+1 ; i++)
	{
		ir  = vk[i];
		id  = vk[i+1] - vk[i]; 
		ird = vk[i+1];
	
		if(id<2)
		{
			if(!m_flML && m_k>=0) m_cah++;
		}
		else
		{
			fli = i < m_nKm;		 
		
			vyi = m_vy[ir:ird];	//r:r+d, (id+1)*1 vector
			vzi  = vz[ir:ird];
			vzbi = vzb[ir:ird];
	
			if(i)
			{
				vetaoi = vetao[ir:ird];
		
				vho			//current(old) point
				 = vhn		//new point, r:r+d
				 = vh_s[ir] | zeros(id, 1);
		
				for(t=0 ; t<id ; t++)
					vho[t+1] = dphi * vho[t] + vetaoi[t]*dsig;	//update r+1:r+d	
				if(fli) dhrd1 = dphi * vho[id] + vetaoi[id]*dsig;
			}
			else	//initial case
			{
				ir = -1;
				id++;
				vyi = 0 | vyi;				//add y_0
				vzi = dmu_z | vzi;
				vzbi = 0 | vzbi;
				vho = 0 | vh_s[:ird];		//add h_0
				vhn = zeros(id+1, 1);
				dhrd1 = vh_s[ird+1];
				vetaoi = vh[0]*sqrt(1-dphi^2)/dsig | vetao[:ird];
			}
	
			/* finding mode	& draw candidate */
			for(j=0 ; j<=nite ; j++)
			{
				vhh = j ? vhn : vho;		//h_hat, r:r+d
				   
				vmu = vyi[:id-1] - dbeta * vzbi[:id-1] .* exp(vhh[:id-1]/2) * dgam
					- drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam .* (vhh[1:] - dphi * vhh[:id-1]) / dsig;
				vsig2 = (1 - drho^2) * vzi .* exp(vhh) * dgam^2;
			
				vder = -0.5 * dbeta * vzbi[:id-1] .* exp(vhh[:id-1]/2) * dgam
					 - drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam
							.* ((vhh[1:] - dphi * vhh[:id-1])/2 - dphi) / dsig;	//dmu_t/dh_t	

				vder1 = - drho * sqrt(vzi[:id-1]) .* exp(vhh[:id-1]/2) * dgam / dsig;
				vder1 = i == 0 ? (0|0|vder1[1:]) : (0|vder1);	//dmu_{t-1}/dh_t

				if(fli)	// if not last obs of last block
				{
					vmu |= vyi[id] - dbeta * vzbi[id] * exp(vhh[id]/2) * dgam
						 - drho * sqrt(vzi[id])	* exp(vhh[id]/2) * dgam	* (dhrd1 - dphi * vhh[id]) / dsig;
					vder |= -0.5 * dbeta * vzbi[id] * exp(vhh[id]/2) * dgam
						  - drho * sqrt(vzi[id]) * exp(vhh[id]/2) * dgam
								 * ((dhrd1 - dphi * vhh[id])/2 - dphi) / dsig;
					dlfhrd1 = -0.5 * (dhrd1 - dphi * vhh[id])^2 / dsig^2;
				}
				else
				{
					vmu |= vyi[id] - dbeta * vzbi[id] * exp(vhh[id]/2) * dgam;
					vsig2[id] = vzi[id] * exp(vhh[id]) * dgam^2;
					vder |= -0.5 * dbeta * vzbi[id] * exp(vhh[id]/2) * dgam;
					dlfhrd1 = 0;
				}
				if(!i) vmu[0] = 0;
				
				//vmu, vsig2, vder, vder1 are t=(r:r+d),(id+1)*1 vector
			
				vdel = -0.5									//delta
					 + 0.5 * vmu[1:].^2 ./ vsig2[1:]
					 - vmu[1:] ./ vsig2[1:] .* vder[1:]
					 - vmu[:id-1] ./ vsig2[:id-1] .* vder1[1:];
				vA = 0.5 + vder[1:].^2 ./ vsig2[1:] + vder1[1:].^2 ./ vsig2[:id-1];
				vB = 0 | (vder[1:id-1] .* vder1[2:] ./ vsig2[1:id-1]);
			
				if(fli)
				{
					vdel[id-1] += dphi * (dhrd1 - dphi * vhh[id]) / dsig^2;
					vA[id-1] += dphi^2 / dsig^2;
				}
					  
				//vdel, vA, vB are t=(r+1:r+d), id*1 vector
	 
				vDh = vJh = vKh = vbh = vZ = zeros(id, 1); //r+1:r+d
				
				vDh[0] = vA[0];
				vKh[0] = sqrt(vDh[0]);
				vbh[0] = vdel[0];
				for(t=1 ; t<id ; t++){
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
										
				for(t=0 ; t<id ; t++)	//Kalman filter
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
							vhn[t+2] = dphi * vhn[t+1] + vu[t];
							 
						//AR step
						dLh = sumc(-0.5 * vhh - 0.5 * vmu.^2 ./ vsig2) + dlfhrd1;
					
						mBh = (0 ~ diag(vB[1:])) | 0;
						mQh = diag(vA) + mBh + mBh';
					
						dpron = dLh + vdel ' (vhn - vhh)[1:]
						 	  - 0.5 * (vhn - vhh)[1:] ' mQh * (vhn - vhh)[1:];
					
						vmun = vyi[:id-1] - dbeta * vzbi[:id-1] .* exp(vhn[:id-1]/2) * dgam
							 - drho * sqrt(vzi[:id-1]) .* exp(vhn[:id-1]/2) * dgam
							 .* (vhn[1:] - dphi * vhn[:id-1]) / dsig;
						vsign2 = (1 - drho^2) * vzi .* exp(vhn) * dgam^2;
	
						if(fli)
						{
							vmun |= vyi[id] - dbeta * vzbi[id] * exp(vhn[id]/2) * dgam
								  - drho * sqrt(vzi[id]) * exp(vhn[id]/2) * dgam
										 * (dhrd1 - dphi * vhn[id]) / dsig;
							dlfhrdn1 = -0.5 * (dhrd1 - dphi * vhn[id])^2 / dsig^2;
						}
						else
						{
							vmun |= vyi[id] - dbeta * vzbi[id] * exp(vhn[id]/2) * dgam;
							vsign2[id] = vzi[id] * exp(vhn[id]) * dgam^2;
							dlfhrdn1 = 0;
						}
						if(!i) vmun[0] = 0;
	
						dposn = sumc(-0.5 * vhn	- 0.5 * vmun.^2 ./ vsign2) + dlfhrdn1;
					
						if(!m_flML && m_k>=0) m_cah++;

						if(car>=100)
							print("h! ");

					}while(ranu(1, 1) >= exp(dposn - dpron) && car<100);
				}	//end draw candidate
			}	//end finding mode & draw candidate	
	
			//MH step
			dproo = dLh + vdel ' (vho - vhh)[1:] - 0.5 * (vho - vhh)[1:] ' mQh * (vho - vhh)[1:];
									  
			vmuo = vyi[:id-1] - dbeta * vzbi[:id-1] .* exp(vho[:id-1]/2) * dgam
				 - drho * sqrt(vzi[:id-1]) .* exp(vho[:id-1]/2) * dgam
				 .* (vho[1:] - dphi * vho[:id-1]) / dsig;
			vsigo2 = (1 - drho^2) * vzi .* exp(vho) * dgam^2;
			if(fli)
			{
				vmuo |= vyi[id] - dbeta * vzbi[id] * exp(vho[id]/2) * dgam
					  - drho * sqrt(vzi[id]) * exp(vho[id]/2) * dgam
							 * (dhrd1 - dphi * vho[id]) / dsig;
				dlfhrdo1 = -0.5 * (dhrd1 - dphi * vho[id])^2 / dsig^2;
			}
			else{
				vmuo |= vyi[id] - dbeta * vzbi[id] * exp(vho[id]/2) * dgam;
				vsigo2[id] = vzi[id] * exp(vho[id]) * dgam^2;
				dlfhrdo1 = 0;
			}
			if(!i) vmuo[0] = 0;
			 
			dposo = sumc(-0.5 * vho	- 0.5 * vmuo.^2 ./ vsigo2) + dlfhrdo1;
	
			dfrac = exp(  dposn + min(dposo|dproo) - dposo - min(dposn|dpron) );
									
			if(ranu(1, 1) < dfrac)
			{
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
SVGH::ML(const nI, const nR, const nsim)
{
	decl vlik, dpos, dposse, dpri, iBm, i;
	decl dvm, vfl, mlsh, mOm, mOmsum, vtaua;
	decl mml;
	decl times = timer();

	m_vposm = loadmat(sprint(m_sdirectory, "/posstat_",m_sfilename,".csv"))[][0];
		//posterior mean of phi, sigma, rho, mu, beta (if skew t), nu (if t)
	iBm = min(500, nsim/2);			//bandwidth

	/*--- likelihood part ---*/
	
	vlik = PFilter(nI, nR);

	
	/*--- posterior part ---*/

	m_flML = 1;
	MCMC(nsim, 0.1*nsim);		//reduced run
	m_flML = 0;

	dvm = meanc(m_mmls);
	vfl = <1; -1; 1; -1; 1>;	// phi, (sigma,rho), mu
	if(m_flt)
	{
		if(m_flb)
			vfl |= <1; -1>;	// beta

		vfl |= <1; -1>;	// nu
	}
	dpos = log(dvm) * vfl;
						 
	mlsh = m_mmls - dvm;
	mOmsum = zeros(5+2*m_flt+2*m_flb, 5+2*m_flt+2*m_flb);
	for(i=0 ; i<iBm ; i++){
		mOm = mlsh[i:][] ' mlsh[:nsim-i-1][] / nsim;
		mOmsum += (1 - i/(iBm+1)) * (mOm + mOm' * (i>0));
	}

	mOmsum /= nsim;
	vtaua = vfl ./ dvm';
	dposse = sqrt(vtaua ' mOmsum * vtaua);


	/*--- prior part ---*/

	dpri = log(0.5*densbeta((m_vposm[0]+1)/2,
							m_dap0, m_dbp0))
		 + log(1/m_vposm[1]^4
		 	   * densgamma(1/m_vposm[1]^2, m_dn0, m_dS0))
		 + log(0.5*densbeta((m_vposm[2]+1)/2,
		 					m_dar0, m_dbr0))
		 + log(densn((m_vposm[3]-m_dmm0)/m_dsm0)
		 	   / m_dsm0);
	if(m_flt)
	{
		dpri += log(densgamma(m_vposm[4+m_flb], m_dnn0, m_dnS0) / (1 - probgamma(4, m_dnn0, m_dnS0)));
		if(m_flb)
			dpri +=	log(densn((m_vposm[4]-m_db0)/sqrt(m_dsb02)) / sqrt(m_dsb02));
	}

	/*--- output ---*/

	print("\n\n========================================"
		  "\n     Marginal likelihood estimation"
		  "\n========================================"
		  "\n                                 svgh.ox"
		  "\nModel: SV");
	if(m_flt)
	{
		if(m_flb) print("SK");
		print("t model");
	}
	else
		print("n model");
	print("\n\n");
//	print("Prior: (phi+1)/2  ~ Beta(", m_dap0, ", ", m_dbp0, ")\n");
//	print("       sigma^(-2) ~ Gamma(", m_dn0, ", ", m_dS0, ")\n");
//	print("       (rho+1)/2  ~ Beta(", m_dar0, ", ", m_dbr0, ")\n");
//	print("       mu         ~ N(", m_dmm0, ", ", m_dsm0^2, ")\n");
//	if(m_flt)
//	{
//		if(m_flb)
//			print("       beta       ~ N(", m_db0, ", ", m_dsb02, ")\n");
//		print("       nu         ~ Gamma(", m_dnn0, ", ", m_dnS0, ")\n");
//		print("%c", {"Mean", "Stdev"}, "%r", {"       phi", "       sigma^(-2)", "       rho", "       nu"}, "%cf", {"%7.2f"},
//		  fBetaMoment(m_dap0,m_dbp0,1)*2-1 ~ sqrt(fBetaMoment(m_dap0,m_dbp0,2))*2 |
//		  fGammaMoment(m_dn0,m_dS0,1) ~ sqrt(fGammaMoment(m_dn0,m_dS0,2)) |
//		  fBetaMoment(m_dar0,m_dbr0,1)*2-1 ~ sqrt(fBetaMoment(m_dar0,m_dbr0,2))*2 |
//		  fGammaMoment(m_dnn0,m_dnS0,1) ~ sqrt(fGammaMoment(m_dnn0,m_dnS0,2)));
//	}
//	else
//	{
//		print("%c", {"Mean", "Stdev"}, "%r", {"       phi", "       sigma^(-2)", "       rho"}, "%cf", {"%7.2f"},
//		  fBetaMoment(m_dap0,m_dbp0,1)*2-1 ~ sqrt(fBetaMoment(m_dap0,m_dbp0,2))*2 |
//		  fGammaMoment(m_dn0,m_dS0,1) ~ sqrt(fGammaMoment(m_dn0,m_dS0,2)) |
//		  fBetaMoment(m_dar0,m_dbr0,1)*2-1 ~ sqrt(fBetaMoment(m_dar0,m_dbr0,2))*2 );
//	}
	print("\n    ESTIMATION RESULT    (S.E.)"
		  "\n--------------------------------");
	print("%r",{" Likelihood", " Prior", " Posterior",
			    " log-ML"},
		  "%cf", {"%10.3f", "%8.3f"},
		    (meanc(vlik) ~ sqrt(varc(vlik)))
		  | (dpri ~ 0)
		  | (dpos ~ dposse)
		  | ((meanc(vlik)+dpri-dpos)
		  	  ~ sqrt(varc(vlik)+dposse^2)));
	print("--------------------------------"
		  "\n[Particle filter]\nParticle", "%7.0f", nI,
		  "\nRepricants", "%5.0f", nR,
		  "\n\n[Reduced MCMC]\nSample", "%9.0f", nsim);
  
	mml = (meanc(vlik) ~ sqrt(varc(vlik)))
		  | (dpri ~ 0)
		  | (dpos ~ dposse)
		  | ((meanc(vlik)+dpri-dpos) ~ sqrt(varc(vlik)+dposse^2));

	savemat(sprint(m_sdirectory, "/ml_",m_sfilename,".csv"), mml);

	println("\n\nTime: ", timespan(times));
}

/*
** PFilter estimates likelihood by particle filter
*/
SVGH::PFilter(const nI, const nR)	// nI: number of particles
{
	decl dphi, dsig, drho, dmu, dbeta, dnu, dmuz, dsigz, db;
	decl vf, vlik, vz, vh, vw, vfhat, vmu, vgj, vuc, vuci;
	decl vzi, vind, cuci, cuc, vhi, vfhati, vgji, i, j, t;

	dphi  = m_vposm[0];
	dsig  = m_vposm[1];
	drho  = m_vposm[2];
	dmu   = m_vposm[3];
	dbeta = m_flb ? m_vposm[4] : 0;
	dnu   = m_flt ? m_vposm[4+m_flb] : M_INF_POS;

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

		vh = dmu + rann(nI, 1) * dsig / sqrt(1 - dphi^2); // generate h_1
	
		vw = sqrt(db)*densn((sqrt(db)*m_vy[0] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2))) ./ (sqrt(vz) .* exp(vh/2));
	
		vf[0] = meanc(vw);
		vfhat = vw / sumc(vw);
		
		/*---------- Filtering START -----------*/
		for(t=0 ; t<m_ns-1 ; t++)
		{
			vmu = dmu + dphi * (vh - dmu) + drho * dsig * (sqrt(db)*m_vy[t] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2));
			vgj = sqrt(db)*densn(sqrt(db)*m_vy[t+1]./(sqrt(dmuz)*exp(vmu/2))) ./ (sqrt(dmuz)*exp(vmu/2)) .* vfhat;
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

			vmu = dmu + dphi * (vhi - dmu) + drho * dsig * (sqrt(db)*m_vy[t] - dbeta*(vzi-dmuz).*exp(vhi/2)) ./ (sqrt(vzi) .* exp(vhi/2));
			vh = vmu + rann(nI, 1) * dsig * sqrt(1-drho^2);
			
			if(m_flt)
				vz = 1 ./ rangamma(nI, 1, dnu/2, dnu/2);
				
			vw = sqrt(db)*densn((sqrt(db)*m_vy[t+1] - dbeta*(vz-dmuz).*exp(vh/2)) ./ (sqrt(vz) .* exp(vh/2)))
				 ./ (sqrt(vz) .* exp(vh/2)) .* vfhati ./ vgji;
		
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
SVGH::fTsvar(const vx, const iBm)
{
	decl ns, vsp;

	ns = rows(vx);
	vsp = periodogram(vx, iBm, 2, 1) / ns;
	
	return M_2PI * vsp[0];
}

SVGH::fGeweke(const vx, const iBm)
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


SVGH::fBetaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(da+db));
	else if(cm==2) return(da*db/((da+db)^2*(da+db+1)));
}

SVGH::fGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/db);
	else if(cm==2) return(da/db^2);
}

SVGH::fInvGammaMoment(const da, const db, const cm)
{
	if(cm==1) return(da/(db-1));
	else if(cm==2) return(db^2/(sqr(da-1)*(da-2)));
}