// run_estimation.ox : entry point
#include <oxstd.h>
#import  "utils.ox"
#import  "simulate_data.ox"
#import  "models.oxh"
#import  "models_impl.ox"
#import  "postproc.ox"

main()
{
    decl args = ParseArgs();        // --mode, --config, --data_dir, --seed, etc.
    SetSeed(args.seed);

    println("== IJF_Revision_Ox ==");
    println("mode: " + args.mode);

    decl Y, X;

    if (Streq(args.mode, "synthetic")) {
        println(">> Running SYNTHETIC pipeline");
        (Y, X) = SimulateDataset(args);
    } else if (Streq(args.mode, "real")) {
        println(">> Running REAL-DATA pipeline");
        (Y, X) = LoadRealData(args);
    } else {
        oxrunerror("Unknown --mode. Use synthetic or real.");
    }

    // Fit models (skeleton)
    decl fit_rsv = FitRSV(Y, X, args);
    decl fit_eg  = FitEGARCH_GED(Y, args);

    // Forecast & evaluate
    ForecastAndEvaluate(fit_rsv, fit_eg, args);

    // Export tables/figures
    ExportFiguresAndTables(fit_rsv, fit_eg, args);

    println(">> Done. See results/ for outputs.");
}
