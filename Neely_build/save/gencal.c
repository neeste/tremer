//  gencal - create SAPP calibration from Newick file

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>

#define MAXFNL 1024
#define MAXLEN 8192

int
group_number(char *nfn)
{
   while (*nfn && !isdigit(*nfn)) nfn++;
   return (atoi(nfn));
}

void
calibrate_nodes(char *nfn)
{
    char line[MAXLEN], label[20], ofn[80];
    int i, j, nc, nn, yr, gn;
    FILE *nfp, *ofp;

    gn = group_number(nfn);
    sprintf(ofn, "calibrate%d.txt", gn);
    ofp = fopen(ofn, "w");
    nfp = fopen(nfn, "r");
    fgets(line, MAXLEN, nfp);
    nc = strlen(line);
    for (i = 0; i < nc; i++) {
        while (line[i] && (line[i] != 'N')) i++;
        if (i == nc) break;
        if (strncmp(line + i, "Node-#", 6) == 0) {
            i += 6;
            nn = atoi(line + i);
            while ((line[i] >= '0') && (line[i] <= '9')) i++;
            if (line[i] != '{') continue;
            i += 2;
            while (line[i] && (line[i] != '"')) i++;
            if (i == nc) break;
            if (line[i - 5] == '.') {
                // copy node label
                j = i - 5;
                while ((line[j] != '-') && (line[j] != '"')) j--;
                strncpy(label, line + j + 1, 20);
                j = 0;
                while(label[j] != '.') j++;
                label[j] = 0;
                // print year
                yr = atoi(line + i - 4);
                if (yr) {
                    fprintf(ofp, "NODE %3d %d * %s\n", nn, yr, label);
                }
            }
        }
    }
    fclose(nfp);
    fclose(ofp);
}

int
main(int argc, char **argv)
{
    char nfn[MAXFNL];

    if (argc < 2) {
        printf("usage: gencal newick_file\n");
	exit(1);
    }
    while (argc > 1) {
        strncpy(nfn, argv[1], MAXFNL);
        calibrate_nodes(nfn);
        argv++;
        argc--;
    }
    exit(0);
}

