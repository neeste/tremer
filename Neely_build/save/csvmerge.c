/* csvmerge.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
csv_read(char *fn, char **kitnam, int **kitstr)
{
    char s[8192], *b, *e;
    int nlen, nstr, nint;
    FILE *fp;

    strcpy(s, fn);
    e = strstr(s, "_YDNA_");
    *e = 0;
    for (b = e; b > s; b--) if (b[-1] == '/') break;
    *kitnam = strdup(b);

    fp = fopen(fn, "r");
    fgets(s, 8192, fp);
    fgets(s, 8192, fp);
    nlen = strlen(s);
    nstr = 1;
    for (b = s; b < (s + nlen); b++) {
        if (b[0] == ',') nstr++;
        if ((b[0] != ' ') && (b[1] == '-')) nstr++;
    }
    *kitstr = calloc(nstr, sizeof(int));
    nint = 0;
    for (b = s; b < (s + nlen); b++) {
        if ((*b >= '0') && (*b <= '9')) {     // STR repetitions
            kitstr[0][nint] = atoi(b);
            nint++;
            while ((*b >= '0') && (*b <= '9')) b++;
        }
        if ((b[0] == ' ') && (b[1] <= '-')) { // missing value
            kitstr[0][nint] = 0;
            nint++;
            b += 2;
        }
    }
    if (nstr != nint) {
        fprintf(stderr, "ERROR: nstr=%d nint=%d\n", nstr, nint);
    }
    return (nstr);
}

int
main(int ac, char **av)
{
    char *kitnam[512];
    int csvcnt, i, j, n, strcnt[512], *kitstr[512];
    FILE *fp;

    csvcnt = 0;
    while (ac > 1) {
        strcnt[csvcnt] = csv_read(av[1], kitnam + csvcnt, kitstr + csvcnt);
        csvcnt++;
        ac--;
        av++;
    }
    printf("Found %d CVS files.\n", csvcnt);
    // write download file
    fp = fopen("csvmerge.txt", "w");
    fprintf(fp, "/STRDATA\n");
    for (i = 0; i < csvcnt; i++) {
        fprintf(fp, "* Kit %s has %d STRs.\n",kitnam[i], strcnt[i]);
        fprintf(fp, "%s",kitnam[i]);
        for (j = 0; j < strcnt[i]; j++) {
            n = kitstr[i][j];
            if (n) {
                fprintf(fp, " %d", kitstr[i][j]);
            } else {
                fprintf(fp, " N");
            }
        }
        fprintf(fp, "\n");
        free(kitnam[i]);
        free(kitstr[i]);
    }
    fclose(fp);
    return (0);
}
