//
//
#include <oxstd.h>
#include <oxfloat.h>
#include <oxprob.h>
#include <oxdraw.h>
#include <ReportMCMC.h>

//#ifndef	ReportMCMC_INCLUDED
#define ReportMCMC_INCLUDED

ReportMCMC::ReportMCMC(mdata)
{
	decl i;
	/* Version */
	m_sversion = "\n [Report MCMC Class] (Last Updated 2006/11/08)"; 

	/*	 Following paragraph may be omitted when we always use SetData */
	
	m_msample = mdata;
	mdata = <>;

	m_cdim = columns(m_msample);
	m_cRep = rows(m_msample);	  
	// Default bandwidth
	m_dBm = 0.1*m_cRep;

	m_asname = {"Param1"};
	if( m_cdim > 1){
		m_asname = {};
    	for(i=0; i<m_cdim; ++i){
			m_asname ~= sprint("Param", i+1);
		}
	}	

	/*                                                                */
	
	// Default format
	m_sfmt = "%#13.5g";
	m_sfigtype = "eps";
	m_soutfilename = "Output";
	m_fcorr = 1;   // print correlation matrix of parameters


	m_cmaxnumfig =  16;  // each figure file has at most 16 figures 
	m_cpage   = floor((m_cdim-1)/m_cmaxnumfig)+1;	

}


	
ReportMCMC::Report()
{

	DrawDensityPlot();
	DrawACFPlot();
	DrawSamplePath();
	CalculateStatistics();
    OutputFile();
}

ReportMCMC::DrawDensityPlot()
{
	decl i, j, k, asfilename;
	// sfilename = sprint(m_soutfilename,"_Density_it_", m_cRep,".",m_sfigtype);
	asfilename = new array[m_cpage];

	if(m_cpage == 1){
		asfilename[0] = sprint(m_soutfilename,"_Density_it_", m_cRep,".",m_sfigtype);
	}
	else{
		for( i = 0; i < m_cpage; i++){
			asfilename[i] = sprint(m_soutfilename,"_Density_",i+1,"_it_", m_cRep,".",m_sfigtype);
		}
	}

	for( i = 0; i < m_cpage; i++){
		if(i < m_cpage-1){
			k = 0;
			for(j = i*m_cmaxnumfig; j < (i+1)*m_cmaxnumfig; j++){
				DrawDensity(k, m_msample[][j]', m_asname[j],1,0,0,0,0,0,1);
//    			DrawTitle(k,"Posterior Density");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}

		else{
			k = 0;
			for(j = i*m_cmaxnumfig; j < m_cdim; j++){
				DrawDensity(k, m_msample[][j]', m_asname[j],1,0,0,0,0,0,1);
    		//	DrawTitle(k,"Posterior Density");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}

	}
}

ReportMCMC::DrawACFPlot(){
 //acf
	decl in_acf  = 1.5*m_dBm;
	decl cn_draw = 1000;  // number of points to draw ACF
	decl vacf;
	decl i, j, k, asfilename;
	asfilename = new array[m_cpage];
	// sfilename = sprint(m_soutfilename,"_ACF_it_",m_cRep,".",m_sfigtype);

	if(m_cpage == 1){
		asfilename[0] = sprint(m_soutfilename,"_ACF_it_", m_cRep,".",m_sfigtype);
	}
	else{
		for( i = 0; i < m_cpage; i++){
			asfilename[i] = sprint(m_soutfilename,"_ACF_",i+1,"_it_", m_cRep,".",m_sfigtype);
		}
	}

	for( i = 0; i < m_cpage; i++){
		if(i < m_cpage-1){
			k = 0;
			for(j = i*m_cmaxnumfig; j < (i+1)*m_cmaxnumfig; j++){
				if(in_acf <= cn_draw){
					DrawAcf(k, m_msample[][j]', m_asname[j], in_acf);
				}
				else{
					vacf = acf(m_msample[][j], in_acf);
					vacf = thinc(vacf[1:]', cn_draw);
					DrawMatrix(k, vacf, m_asname[j], 1, in_acf/cn_draw);
					DrawAdjust(ADJ_INDEX, 1, 0);
				}
				//	DrawTitle(k, "Correlogram");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}
		else{
			k = 0;
			for(j = i*m_cmaxnumfig; j < m_cdim; j++){
				if(in_acf <= cn_draw){
					DrawAcf(k, m_msample[][j]', m_asname[j], in_acf);
				}
				else{
					vacf = acf(m_msample[][j], in_acf);
					vacf = thinc(vacf[1:]', cn_draw);
					DrawMatrix(k, vacf, m_asname[j], 1, in_acf/cn_draw);
					DrawAdjust(ADJ_INDEX, 1, 0);
				}
				//	DrawTitle(k, "Correlogram");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}


	}
	
}

ReportMCMC::DrawSamplePath()
{
	decl sfilename, slabel;
	decl cn_draw, cc;   // number of points to draw sample path
	decl cstep = 1;
	cn_draw = cc = 1000;
	decl i, j, k, asfilename;
	asfilename = new array[m_cpage];
	//sfilename = sprint(m_soutfilename,"_SamplePath_it_", m_cRep,".",m_sfigtype);

	if(m_cpage == 1){
		asfilename[0] = sprint(m_soutfilename,"_SamplePath_it_", m_cRep,".",m_sfigtype);
	}
	else{
		for( i = 0; i < m_cpage; i++){
			asfilename[i] = sprint(m_soutfilename,"_SamplePath_",i+1,"_it_", m_cRep,".",m_sfigtype);
		}
	}


	for( i = 0; i < m_cpage; i++){
		if(i < m_cpage-1){
			k = 0;
			for(j = i*m_cmaxnumfig; j < (i+1)*m_cmaxnumfig; j++){
		    	if(m_cRep <= cn_draw)
					{
     					DrawTMatrix(k, m_msample[][j]', m_asname[j],1,1,1);
					}
				else{
					if(fmod(m_cRep, cn_draw) == 0)
						{
							cstep =  idiv(m_cRep, cn_draw);	
						}
					else{
							cstep =  idiv(m_cRep, cn_draw) + 1;
						}
					cc =  idiv( m_cRep, cstep);

					DrawXMatrix(k, thinc((m_msample[][j]') , cc), m_asname[j],
							range(1, m_cRep, cstep), "", 0,  1);
    			}
				// DrawTitle(k, "Sample Path");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}
		else{
			k = 0;
			for(j = i*m_cmaxnumfig; j < m_cdim; j++){
		    	if(m_cRep <= cn_draw){	
     				DrawTMatrix(k, m_msample[][j]', m_asname[j],1,1,1);
				}
				else{
					if(fmod(m_cRep, cn_draw) == 0)
						{
							cstep =  idiv(m_cRep, cn_draw);	
						}
					else{
							cstep =  idiv(m_cRep, cn_draw) + 1;
						}
					cc =  idiv( m_cRep, cstep);

					DrawXMatrix(k, thinc((m_msample[][j]') , cc), m_asname[j],
							range(1, m_cRep, cstep), "", 0,  1);
    			}
				// DrawTitle(k, "Sample Path");
				k += 1;
			}
			SaveDrawWindow(asfilename[i]);		  
			CloseDrawWindow();
		}
	}

}

ReportMCMC::CalculateStatistics()
{
 //individual
	decl i, sfilename, amstat, asitem, asname;
	decl asname_, vflg;
	
	asname = asname_ = m_asname;
	asitem = new array[2];
	amstat = new array[3];

	
	// Check whether 0 is outside 95% (*) or 99% (**) credible interval
	vflg = fCI(m_msample);
	for(i = 0; i < m_cdim; i++){
//		if(vflg[i] == 0){asname_[i]  = sprint(asname[i], "   ") ;}
//		else if(vflg[i] == 1){asname_[i]  = sprint(asname[i], "  *") ;}
//		else if(vflg[i] == 2){asname_[i] = sprint(asname[i], " **");}
		if(vflg[i] == 0){asname_[i]  = sprint("    ", asname[i]) ;}
		else if(vflg[i] == 1){asname_[i]  = sprint("*   ", asname[i]) ;}
		else if(vflg[i] == 2){asname_[i] = sprint("**  ", asname[i]);}
		}

	amstat[0] = meanc(m_msample)'~( varc(m_msample).^0.5 )'
	    	~quantilec(m_msample, <0.025,0.5,0.975>)' ;

	amstat[1] = fGeweke(m_msample, m_dBm)~fTsvar(m_msample, m_dBm)./varc(m_msample)';

	amstat[2] = upper(correlation(m_msample)) - unit(m_cdim)
			  + lower(variance(m_msample));  

	asitem[0] = {"Mean", "Stdev", "95%L", "Median", "95%U"};
	asitem[1] = {"Geweke", "Inef"};

	println("\n -----------------------------Summary Statistics:",
					"-----------------------------");
	println(m_sversion);
	println( "\n  Date:",timestr(today()));
	println("\n  Number of obs. is ",m_cRep);
	println( "\n [Basic Summary Statistics] ");
	println( "%r", asname_, "%c", asitem[0], amstat[0]);
	println( "\n Note: *(**) -> 95(99)% credible interval does not contain 0.");
	
	println( "\n [Diagnostic Statistics] ");
	println( "%r", asname, "%c", asitem[1],  amstat[1]);	
	if(m_fcorr == 1){
		println( "\n [Covariance (Lower), Variance (Diagonal) and Correlation (Upper)]" );
		println( "%r", asname, "%c", asname, amstat[2]);
	}
	
	decl file;
	sfilename  = sprint(m_soutfilename,"_Result_it_",m_cRep,".txt");
	file=fopen(sfilename, "a");
	fprintln(file, "\n -----------------------------Summary Statistics:",
					"-----------------------------");
	fprintln(file, m_sversion);
	fprintln(file, "\n  Date:",timestr(today()));
	fprintln(file, "\n  Number of obs. is ",m_cRep);
	fprintln(file, "\n [Basic Summary Statistics] ");
	fprintln(file, "%r", asname_, "%c", asitem[0], amstat[0]);
	fprintln(file,  "\n Note: *(**) -> 95(99)% credible interval does not contain 0.");

	fprintln(file, "\n [Diagnostic Statistics] ");
	fprintln(file, "%r", asname, "%c", asitem[1], amstat[1]);
	if(m_fcorr == 1){
		fprintln(file, "\n [Covariance (Lower), Variance (Diagonal) and Correlation (Upper)]" );
		fprintln(file, "%r", asname, "%c", asname, amstat[2]);
	}
	fclose(file);
}
 
ReportMCMC::OutputFile()
{
	format(1024); // max number of columns
	decl sfilename, file;
	sfilename = sprint(m_soutfilename,"_Samples_it_",m_cRep,".dat");
	file = fopen(sfilename,"w");
	fprint(file, m_sfmt, m_msample);
	fclose(file);
}

ReportMCMC::fTsvar(const mX, const dBm)
{
	// mX: T x m matrix (T: number of obs, m: number of series) 
	// Calculate Spectral density function using
	// Bandwidth=dBm at 2 points (0, pi) and Parzen	Window
	// 1st row of msp is spectral density evaluated at 0
	decl msp = periodogram(mX, dBm, 2, 1)/rows(mX);
	return M_2PI*(msp[0][])';
}

ReportMCMC::fTsvar_Batch(const mX, const dBm)
{
	// n=dBm*M  dBm:# of obs within each batch, M:# of batches
	decl i, cnobs, cM, my, vxbar;
	cnobs = rows(mX);
	cM = cnobs/dBm;
	
	vxbar=meanc(mX);
	my=zeros(cM,columns(mX));
	for(i=0; i<cM; ++i){
		my[i][]=meanc(mX[i*dBm:(i+1)*dBm-1][]);
	}
	return dBm/(cM-1)*sumsqrc(my-vxbar)';
}
 

ReportMCMC::fGeweke(const mX, const dBm)
{
 // mX: T x m matrix (T: number of obs, m: number of series) 
	decl cT, cn1, cn2, vvar1, vvar2, mX1, mX2, vx1bar, vx2bar, vz;
	cT = rows(mX);
	cn1 = floor(0.1*cT);
	cn2 = floor(0.5*cT);

	mX1 = mX[:cn1-1][];
	mX2 = mX[cT-cn2:][];
	
	vx1bar = meanc(mX1)';
	vx2bar = meanc(mX2)';

	vvar1=fTsvar(mX1,dBm);
	vvar2=fTsvar(mX2,dBm);

	vz=(vx1bar - vx2bar) ./(vvar1/cn1+vvar2/cn2).^(0.5);

	return 2*tailn(fabs(vz));
} 

ReportMCMC::fCI(const mX)
{
	// Return 0 if 0 is  inide 95% Credible	Interval
	//        1 if 0 is outide 95% Credible	Interval
	//        2 if 0 is outide 99% Credible	Interval
 	decl vflg1, vflg2, mCI;
	mCI = quantilec(mX, <0.025, 0.975>)' ;
	vflg1 = (mCI[][0] .> 0 .|| mCI[][1] .< 0) .?  1 .: 0; 
	mCI = quantilec(mX, <0.005, 0.995>)' ;
	vflg2 = (mCI[][0] .> 0 .|| mCI[][1] .< 0) .?  1 .: 0; 
	return vflg1+vflg2;
}
ReportMCMC::SetBandwidth(const dBm){ m_dBm = dBm;}
ReportMCMC::SetFormat(const sformat){	m_sfmt = sformat;}
ReportMCMC::SetFigureType(const sfigtype){ m_sfigtype = sfigtype;}
ReportMCMC::SetOutfileName(const soutfilename){ m_soutfilename = soutfilename;}
ReportMCMC::SetVarNames(const asname){ m_asname = asname;}
ReportMCMC::SetCorrOutput(const fcorr){ m_fcorr = fcorr; }
ReportMCMC::SetMaxNumFig(const cmaxnumfig){ m_cmaxnumfig = cmaxnumfig;}

//#endif
