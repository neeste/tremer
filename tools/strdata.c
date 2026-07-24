// strdata - rewrite text file with included STRDATA text

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>

#define MAXFNL 1024
#define MAXLEN 8192

int
main(int argc, char **argv)
{
    char ifn[MAXFNL],nfn[MAXFNL], line[MAXLEN];
    int i, j, cmp1, cmp2;
    FILE *ifp, *tfp, *nfp;
    static char tfn[MAXFNL] = "strdata.txt";

    if (argc < 2) {
        printf("usage: strdata filename\n");
	exit(1);
    }
    while (argc > 2) {
        if (argv[1][0] == '-') {
            if (argv[1][1] == 'o') {
                if (argc > 2) {
                    strncpy(tfn, argv[2], MAXFNL);
                    argv++;
                    argc--;
                }
            }
        }
        argv++;
        argc--;
    }
    strncpy(ifn, argv[1], MAXFNL);
    printf("output file: %s\n",tfn);
    ifp = fopen(ifn, "r");
    tfp = fopen(tfn, "w");
    while (fgets(line, MAXLEN, ifp)) {
        cmp1 = strncasecmp(line, "/STRDATA", 8);
        cmp2 = strncasecmp(line, "/STRTREE", 8);
        i = 0;
        while (line[i] > ' ') i++;
        while(line[i] == ' ') i++;
        if ((!cmp1 || !cmp2)  && (line[i] == '*')) {
            i++;
            while(line[i] == ' ') i++;
            j = 0;
            while (line[i] > ' ') nfn[j++] = line[i++];
            nfn[j] = 0;
            strcat(nfn, ".txt");
            printf("include %s\n", nfn);
            nfp = fopen(nfn, "r");
            if (nfp == NULL) {
                char alt[MAXFNL];
                snprintf(alt, MAXFNL, "strdata_out/%s", nfn);
                nfp = fopen(alt, "r");
            }
            if (nfp == NULL) {
                printf("WARNING: can't open %s\n",nfn);
            } else {
                fputs(line, tfp);         // copy input line
                while (fgets(line, MAXLEN, nfp)) {
                    fputs(line, tfp);
                }
                while (fgets(line, MAXLEN, ifp)) {
                    if (line[0] == '/') break;
                    if (!cmp2 && (line[0] == '*')) {
                        //fputs(line, tfp);
                    }
                }
                fclose(nfp);
            }
        }
        fputs(line, tfp);
    }
    fclose(ifp);
    fclose(tfp);

    exit(0);
}
