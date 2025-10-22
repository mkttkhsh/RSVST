// models_impl.ox : implement RSV & EGARCH-GED (skeleton)
// Replace placeholders with your likelihoods / MCMC routines

#include <oxstd.h>
#import "utils.ox"

FitRSV(Y, X, args)
{
    decl out = new matrix(1,1);
    out[0][0] = 1.0; // placeholder
    println("RSV: [skeleton] fitting done.");
    return out;
}

FitEGARCH_GED(Y, args)
{
    decl out = new matrix(1,1);
    out[0][0] = 1.0; // placeholder
    println("EGARCH-GED: [skeleton] fitting done.");
    return out;
}

ForecastAndEvaluate(fit_rsv, fit_eg, args)
{
    println("Evaluation: [skeleton] DM(1995) etc.");
}

ExportFiguresAndTables(fit_rsv, fit_eg, args)
{
    println("Export: [skeleton] write CSV/PNG/TEX to results/");
}
