/*
**
**  Reference:
**    Takahashi, M., Yamauchi, Y., Watanabe, T., & Omori, Y. (2025).
**    Realized Stochastic Volatility Model with Skew-t Distributions for
**    Improved Volatility and Quantile Forecasting. arXiv:2501.13179
**    https://arxiv.org/abs/2501.13179
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <oxfloat.h>
#import <maximize>

main()
{
	decl sdata_dir;
	decl sdata, ndata, idata, sdir, nprd;
	decl flmm, smodels, nmodels;
	decl i, smodel, stemp, t;
	decl mvar, mes, vh1, mvolf, mposstat, mposh, maprob;

	sdata = {"N225"};
	ndata = sizec(sdata);

	flmm = 1;
 // smodels = {"rsv_az_n", "rsv_az_t", "rsv_az_sn", "rsv_az_st"};
    smodels = {"sv_az_n", "sv_az_t", "sv_az_sn", "sv_az_st"};
	nmodels = sizec(smodels);

	for(idata = 0; idata < ndata; idata++){
		println("----------------------------------------------");
		println("Data:  ", sdata[idata]);
		sdir = sprint(sdata[idata], "_prediction_2020/files/");
        nprd = 242;
	
		for(i=0; i<nmodels; i++){
			smodel = smodels[i];
			println("----------------------------------------------");
			println("Model:  ", smodel);
	
			mvar = mes = mvolf = mposstat = mposh = maprob = <>;
			for(t=0; t<nprd; t++)
			{
				stemp = sprint(smodel, "_", t);
				// quantile forecasts
				mvar |= loadmat(sprint(sdir, "var_", stemp, ".csv"));
				mes |= loadmat(sprint(sdir, "es_", stemp, ".csv"));
				// volatility forecasts
				vh1 = loadmat(sprint(sdir, "h1_y_", stemp, ".csv"))[][0];
				mvolf |= loadmat(sprint(sdir, "volf_", stemp, ".csv"))[][0] ~ quantilec(vh1) // mean(h), median(h)
							~ meanc(exp(vh1)) ~ quantilec(exp(vh1)); // mean(exp(h)), median(exp(mh1))
				// estimation results
				mposstat |= loadmat(sprint(sdir, "posstat_", stemp, ".csv"));
//				mposstat |= loadmat(sprint(sdir, "posstatmoments_", stemp, ".csv"));
				mposh ~= loadmat(sprint(sdir, "posh_", stemp, ".csv"))[][0];
				maprob |= loadmat(sprint(sdir, "aprob_", stemp, ".csv"));
			}
	
			println("Size of var: ", sizer(mvar), " ", sizec(mvar));
			println("Size of es: ", sizer(mes), " ", sizec(mes));
			println("Size of volf: ", sizer(mvolf), " ", sizec(mvolf));
			println("Size of posstat: ", sizer(mposstat), " ", sizec(mposstat));
			println("Size of posh: ", sizer(mposh), " ", sizec(mposh));
			println("Size of aprob: ", sizer(maprob), " ", sizec(maprob));
	
			savemat(sprint(sdata[idata], "_prediction_2020/var_", smodel, ".csv"), mvar);
			savemat(sprint(sdata[idata], "_prediction_2020/es_", smodel, ".csv"), mes);
			savemat(sprint(sdata[idata], "_prediction_2020/volf_", smodel, ".csv"), mvolf);
			savemat(sprint(sdata[idata], "_prediction_2020/posstat_", smodel, ".csv"), mposstat);
			savemat(sprint(sdata[idata], "_prediction_2020/posh_", smodel, ".csv"), mposh);
			savemat(sprint(sdata[idata], "_prediction_2020/aprob_", smodel, ".csv"), maprob);
		}
	}
}
