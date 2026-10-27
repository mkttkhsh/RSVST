#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import <maximize>
#include<RSV_skewFS.ox>
// 2026 revision: single full-sample in-sample estimation (main-sample
// period, 2009-06..2019-09, "DJIA-RV-2009"), RSV-FS-ST (skew-t, default
// dist), psi free. Companion to the RSVAZ/RSVGH psi-free fullsample
// scripts -- FS was missing from that set; added here using the class's
// built-in "est" (prediction_size=0) single-run mode, matching the "est"
// block already used by RSV_skewFS_run_psi_n225_2020_st.ox for the
// COVID-period counterpart.

main()
{
	decl dTime = timer();
	ranseed(7);

	decl m_cT = 2596;// DJIA_rcc_rv_2009_2019.csv, main-sample size
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

	chdir("./output_emp_DJIA_fs_sn_psi"); // This is the directory the results will be contained in.
	decl rsv = new RSV_skewFS("DJIA-RV-2009", m_cT, 0, 0, 1);// Select the data to use; prediction_size=0 -> single full-sample "est" run

	rsv.SetIsT(0); // SN (skew-normal) dist
	//rsv.SetIsRSV(0); When using SV (no realized-measure equation), uncomment
	rsv.SetIsPsi(1);  // estimate free psi
	rsv.SetParam(dmu, dphi, drho, dgamma, dsigeta, dxi, dsigu, dnu, vh, vlambda, dpsi);
	rsv.SetPrior(dm0_mu, dsig0_mu, dm0_xi, dsig0_xi,
					 dn0_eta, dS0_eta, dn0_u, dS0_u,
					 da0_phi, db0_phi, da0_rho, db0_rho,
					 dn0_gamma, dS0_gamma, dn0_nu, dS0_nu);
	rsv.MCMC(90000, 60000);// MCMC size, burn-in (FS is much heavier than AZ/GH's 15000/5000)

	decl dCompTime_seconds = (timer() - dTime) / 100;
	decl dCompTime_minutes = idiv(dCompTime_seconds, 60);
	dCompTime_seconds =	imod(dCompTime_seconds, 60);
	decl dCompTime_hours = idiv(dCompTime_minutes, 60);
	dCompTime_minutes =	imod(dCompTime_minutes, 60);

	println("=============================================");
	println("Computation Time: ", dCompTime_hours, "h/", dCompTime_minutes, "m/", dCompTime_seconds, "s");
}
