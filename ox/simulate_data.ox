// simulate_data.ox : synthetic data generator (skeleton)
#include <oxstd.h>
#import "utils.ox"

SimulateDataset(args)
{
    decl T = args.T ? args.T : 4000;
    decl Y = zeros(T,1);
    decl X = zeros(T,2); // e.g., rv and rk

    // Very simple illustrative DGP (replace with your RSV-like generator)
    decl h = zeros(T,1);
    h[0] = 0.0;
    for (decl t = 1; t < T; ++t) {
        h[t] = 0.95 * h[t-1] + 0.1 * rann();
    }
    for (decl t = 0; t < T; ++t) {
        Y[t][0] = exp(0.5 * h[t]) * rann(); // return
        X[t][0] = exp(h[t]) + 0.01 * rann(); // rv (proxy)
        X[t][1] = exp(h[t]) + 0.01 * rann(); // rk (proxy)
    }
    return (Y, X);
}
