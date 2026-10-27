#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<RSV_skewFS.ox>
// Free-psi (slope coefficient of log-RV equation) variant: sensitivity/robustness analysis.
// Switch index/distribution via the data string and SetIsT/SetIsSkew, as in the base run file.

StoI(const s){
	decl i;
	decl dig = sizerc(s);
	//println(dig);
	decl n = 0;
	for(i=0;i<dig;i++){
	n += 10^(dig-i-1)  * (int(s[i]) - int("0"));
	}
	return n;
}

main()
{
	decl dTime = timer();
	ranseed(7);
	
	decl m_cT = 2596;// the size of data; If the dataset changed, change this value.
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
	"pred": prediction. parallel computing;
	*/

	//These are the settings for prediction
	decl number_parallel = 11;
	decl prediction_size = 603;
	decl prediction_number, process_size;
	decl ininumber_process;
	decl isini;
	
	if(mode == "est"){
		chdir("./output_emp_DJIA_2009_psi");// This is the filename of which the results will be contained.
		decl rsv = new RSV_skewFS("DJIA-RV-2009", m_cT, 0, 0, 1);// Select the data to use

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
		chdir("./output_pred_DJIA_2009_psi"); // This is the filename of which the results will be contained.
		decl process_number = StoI(replace(oxfilename(2), "RSV_skewFS_run","")) - 1;//starts from zero, ends with (number_parallel)
		/*
		When predicting, create files such that "RSV_skewFS_run[x].ox" ,e.g. "RSV_skewFS_run1.ox", "RSV_skewFS_run2.ox", ...
		Then, run the codes in parallel.
		*/
		process_size = ((process_number == number_parallel) ? imod(prediction_size, number_parallel) : idiv(prediction_size, number_parallel));
		ininumber_process = ((process_number == number_parallel) ? prediction_size - process_size : process_number * process_size); 
		for(prediction_number = ininumber_process; prediction_number < ininumber_process + process_size; prediction_number++){
			isini = 0 + (prediction_number == ininumber_process) * 1;
			decl rsv = new RSV_skewFS("DJIA-RV-2009", m_cT, prediction_size, prediction_number, isini);
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
		
		decl rsv = new RSV_skewFS("DJIA-RV-2009", m_cT, 0, 0, 1);
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