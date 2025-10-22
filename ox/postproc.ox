// postproc.ox : evaluation and export (skeleton)
#include <oxstd.h>

ExportCSV(filename, M)
{
    decl f = fopen(filename, "w");
    if (!f) { println("Cannot open " + filename); return; }
    decl R = rows(M), C = columns(M);
    for (decl i = 0; i < R; ++i) {
        for (decl j = 0; j < C; ++j) {
            fprint(f, "%g", M[i][j]);
            if (j < C-1) fprint(f, ",");
        }
        fprint(f, "\n");
    }
    fclose(f);
}
