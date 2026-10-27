#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<RSV_skewFS.ox>
// Free-psi variant, N225 COVID-period (2012-2020) sample, RSV-FS-ST (skew-t, default dist).
// Switch index/distribution via the data string and SetIsT/SetIsSkew, as in the base run file.
// IJF R2 (2026): single-process sequential rolling loop (not the parallel
// 11-way chunk split used for the DJIA main-sample FS run) — matches how
// the AZ/GH COVID-period scripts run, for consistency.

main()
{
	decl dTime = timer();
	ranseed(7);
	
	decl m_cT = 2201;// the size of data; If the dataset changed, change this value.
	decl dmu = -0.5;
	decl dphi = 0.9;
	decl drho = -0.4;
	decl dgamma = 0.5;
	decl dsigeta = 0.3;
	decl dxi = -0.2;
	decl dsigu = 0.05;
	decl dnu = 10;
	decl vh = zeros(1, m_cT);
	decl vlambda = zeros(1, m_cT);
	decl dpsi = 1;  // initial value of slope coefficient psi (estimated when free)

	decl dm0_mu = 0;
	decl dsig0_mu = 10000;
	decl dm0_xi = 0;
	decl dsig0_xi = 10000;

	decl dn0_eta = 0.1;
	decl dS0_eta = 0.1;
	decl dn0_u = 0.1;
	decl dS0_u = 0.1;

	decl da0_phi = 1;
	decl db0_phi = 1;
	decl da0_rho = 1;
	decl db0_rho = 1;

	decl dn0_gamma = 0.00001;
	decl dS0_gamma = 0.00001;
	decl dn0_nu = 0.01;
	decl dS0_nu = 0.2;

	decl mode = "pred";// This controls the mode of this code,
	/*
	"est": estimation;
	"ML": calculate the marginal likelihood;
	"pred": prediction, single sequential process (IJF R2 2026: no parallel split).
	*/

	//These are the settings for prediction
	decl prediction_size = 242;
	decl prediction_number;
	decl isini;

	if(mode == "est"){
		chdir("./output_emp_N225_2020_psi_st");// This is the filename of which the results will be contained.
		decl rsv = new RSV_skewFS("Nikkei-RV-2020", m_cT, 0, 0, 1);// Select the data to use

		//rsv.SetIsT(0); When using the normal dist, Uncomment 
		//rsv.SetIsSkew(0); When using the symmetric dist, Uncomment
		rsv.SetIsPsi(1);  // estimate free psi
			rsv.SetParam(dmu, dphi, drho, dgamma, dsigeta, dxi, dsigu, dnu, vh, vlambda, dpsi);
		rsv.SetPrior(dm0_mu, dsig0_mu, dm0_xi, dsig0_xi,
						 dn0_eta, dS0_eta, dn0_u, dS0_u,
						 da0_phi, db0_phi, da0_rho, db0_rho,
						 dn0_gamma, dS0_gamma, dn0_nu, dS0_nu);
		//rsv.SimObs(); When using the simulation data, Uncomment 
		rsv.MCMC(90000, 60000);// MCMC size, burn-in
	}

	if(mode == "pred"){
		chdir("./output_pred_N225_2020_psi_st"); // This is the filename of which the results will be contained.
		for(prediction_number = 0; prediction_number < prediction_size; prediction_number++){
			isini = (prediction_number == 0);
			decl rsv = new RSV_skewFS("Nikkei-RV-2020", m_cT, prediction_size, prediction_number, isini);
			rsv.SetIsPsi(1);  // estimate free psi
			rsv.SetParam(dmu, dphi, drho, dgamma, dsigeta, dxi, dsigu, dnu, vh, vlambda, dpsi);
			rsv.SetPrior(dm0_mu, dsig0_mu, dm0_xi, dsig0_xi,
							 dn0_eta, dS0_eta, dn0_u, dS0_u,
							 da0_phi, db0_phi, da0_rho, db0_rho,
							 dn0_gamma, dS0_gamma, dn0_nu, dS0_nu);
			if(isini == 1){
				rsv.MCMC(90000, 60000);
			}else{
				rsv.MCMC(10000, 5000);	
			}
			rsv.GetParam(&dmu, &dphi, &drho, &dgamma, &dsigeta, &dxi, &dsigu, &dnu, &vh, &vlambda, &dpsi);
			delete rsv;
		}
	}
	if(mode == "ML"){
		/*
		When calculating the marginal likelihood, start the initial value estimated.
		*/
		dgamma = meanc(loadmat("./90000_result_DJIA_2009/dgamma_Samples_it_90000.dat"));
		dmu = meanc(loadmat("./90000_result_DJIA_2009/dmu_Samples_it_90000.dat"));
		dnu = meanc(loadmat("./90000_result_DJIA_2009/dnu_Samples_it_90000.dat"));
		dphi = meanc(loadmat("./90000_result_DJIA_2009/dphi_Samples_it_90000.dat"));
		drho = meanc(loadmat("./90000_result_DJIA_2009/drho_Samples_it_90000.dat"));
		dsigeta = meanc(loadmat("./90000_result_DJIA_2009/dsigeta_Samples_it_90000.dat"));
		dsigu = meanc(loadmat("./90000_result_DJIA_2009/dsigu_Samples_it_90000.dat"));
		dxi = meanc(loadmat("./90000_result_DJIA_2009/dxi_Samples_it_90000.dat"));
		vh = loadmat("./90000_result_DJIA_2009/h_mean.csv")';
		chdir("./output_emp_DJIA_2009_ML");// This is the filename of which the results will be contained.

		dm0_mu = 0;
		dsig0_mu = 10;
		dm0_xi = 0;
		dsig0_xi = 1;

		da0_phi = 20;
		db0_phi = 1.5;
		da0_rho = 1;
		db0_rho = 2;
		
		dS0_eta = 0.1;;
		dn0_eta = 10;
		dS0_u = 0.4;
		dn0_u = 10;
		
		dS0_gamma = (1 / dgamma);
		dn0_gamma = (dS0_gamma * dgamma);
		dS0_nu = 0.5;
		dn0_nu = 5;
		
		decl rsv = new RSV_skewFS("Nikkei-RV-2020", m_cT, 0, 0, 1);
		decl nsim = 50000;
		decl nburn = 0;
		rsv.SetML(sqrt(nsim), dmu, dphi, drho, dgamma, dsigeta, dxi, dsigu, dnu, vh);
		rsv.SetIsPsi(1);  // estimate free psi
			rsv.SetParam(dmu, dphi, drho, dgamma, dsigeta, dxi, dsigu, dnu, vh, vlambda, dpsi);
		rsv.SetPrior(dm0_mu, dsig0_mu, dm0_xi, dsig0_xi,
						 dn0_eta, dS0_eta, dn0_u, dS0_u,
						 da0_phi, db0_phi, da0_rho, db0_rho,
						 dn0_gamma, dS0_gamma, dn0_nu, dS0_nu);
		decl dprior = rsv.CalculatePrior();
		decl dlikelihood = rsv.CalculateLikelihood(3000);
		decl dpost = rsv.CalculateML(nsim, nburn);
		println("Marginal Likelihood = ",dprior + dlikelihood - dpost);
		savemat("result.csv", dprior~dlikelihood~dpost~(dprior+dlikelihood-dpost));
	}
	decl dCompTime_seconds = (timer() - dTime) / 100;
	decl dCompTime_minutes = idiv(dCompTime_seconds, 60);
	dCompTime_seconds =	imod(dCompTime_seconds, 60);
	decl dCompTime_hours = idiv(dCompTime_minutes, 60);
	dCompTime_minutes =	imod(dCompTime_minutes, 60);
	
	println("=============================================");
	//println("Computation Time: ", (timer() - dTime) / 100);
	println("Computation Time: ", dCompTime_hours, "h/", dCompTime_minutes, "m/", dCompTime_seconds, "s");

}