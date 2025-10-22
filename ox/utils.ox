// utils.ox : tiny CLI/config helpers (skeleton)
#include <oxstd.h>

struct Args {
    decl mode;
    decl seed;
    decl T;
    decl outdir;
    decl data_dir;
};

ParseArgs()
{
    decl a = new Args;
    a.mode = "synthetic";
    a.seed = 1234567;
    a.T = 4000;
    a.outdir = "results";
    a.data_dir = "";

    // naive parse, replace with ini reader if needed
    decl ca = GetCmdLine();
    for (decl i = 1; i < sizeof(ca); ++i) {
        if (StartsWith(ca[i], "--mode=")) a.mode = SubStr(ca[i], 7);
        if (StartsWith(ca[i], "--seed=")) a.seed = atoi(SubStr(ca[i], 7));
        if (StartsWith(ca[i], "--data_dir=")) a.data_dir = SubStr(ca[i], 11);
        if (StartsWith(ca[i], "--T=")) a.T = atoi(SubStr(ca[i], 4));
    }
    return a;
}

LoadRealData(args)
{
    // placeholder loader - read your CSVs here (not included in public repo)
    println("LoadRealData: [skeleton] data_dir=" + args.data_dir);
    decl T = 2000;
    decl Y = zeros(T,1);
    decl X = zeros(T,2);
    return (Y, X);
}
