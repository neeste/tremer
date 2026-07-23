// strdif - find differences between two STRDATA files

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>

#define MAXKIT   40
#define MAXHDR   80
#define MAXFNL  256
#define MAXSTR  840
#define MAXLIN 8192

struct STRDAT {
    int nkit, nlab, tnmv;
    short nsvs[MAXKIT], nmvs[MAXKIT];
    char *kits[MAXKIT], *kstrs[MAXKIT];
    char top[MAXHDR], *labs[MAXSTR];
};

int
getval(char *line, char *str)
{
    char *s;
    int nsv = 0;

    nsv = 0;
    s = line;
    // create list of STR values
    while (*s && (nsv < MAXSTR)) {
        while(*s == ' ') s++;
        if (*s == 'N') {
            str[nsv++] = 0;
        } else if ((*s >= '1') && (*s <= '9')) {
            str[nsv++] = (char) atoi(s);
        }
        while(*s && (*s != ' ')) s++;
    }
    return (nsv);
}

int
get_kit(FILE *ifp, char *hdr, char *kit, char *kstr, char *line)
{
    char *s;
    int nsv;

    strncpy(hdr, line, MAXHDR);
    fgets(line, MAXLIN, ifp);
    strncpy(kit, line, 20);
    s = kit;
    while (*s > ' ') s++; 
    *s = 0;
    s = line;
    while (*s > ' ') s++; 
    while (*s == ' ') s++; 
    nsv = getval(s, kstr);
    return (nsv);
}

void
put_kit(FILE *ofp, char *kit, char *kstr, int nsv)
{
    int i;
 
    fprintf(ofp, "* Kit %s has %d STRs.\n", kit, nsv);
    fputs(kit, ofp);
    for (i = 0; i < nsv; i++) {
        fputs(" ", ofp);
        if (kstr[i] == 0) {
	    fputs("N", ofp);
        } else {
	    fprintf(ofp, "%d", kstr[i]);
        }
    }
    fputs("\n", ofp);
}

int
get_strlab(char *strlab[])
{
    char line[MAXLIN], lab[20];
    int i, j, n, nlab;
    FILE *lfp;
    static char lfn[] = "str_label.csv";

    lfp = fopen(lfn, "r");
    fgets(line, MAXLIN, lfp);
    nlab = 0;
    n = strlen(line);
    for (i = 0; i < n; i++) {
        j = 0;
        while ((line[i] > ' ') && (line[i] != ',')) {
            lab[j++] = line[i++];
        }
        lab[j] = 0;
        strlab[nlab] = strdup(lab);
        nlab++;
        //printf("%d %s\n",nlab,lab);
    }
    fclose(lfp);
    return(nlab);
}

int
find_misval(struct STRDAT sd, short *nmvs)
{
    int i, j, nmv, tnmv;

    tnmv = 0;
    for (i = 0; i < sd.nkit; i++) {
        nmv = 0;
        for (j = 0; j < sd.nsvs[i]; j++) {
            if (sd.kstrs[i][j] == 0) nmv++;
        }
        nmvs[i] = nmv;
        tnmv += nmv;
    }
    return (tnmv);
}

struct STRDAT
get_strdat(char *ifn)
{
    char hdr[MAXHDR], line[MAXLIN], kit[20], kstr[MAXSTR];
    int j, nkit, nsv;
    FILE *ifp;
    struct STRDAT sd;

    ifp = fopen(ifn, "r"); // input file
    fgets(sd.top, MAXHDR, ifp);
    nkit = 0;
    while (fgets(line, MAXLIN, ifp)) {
        nsv = get_kit(ifp, hdr, kit, kstr, line);
        sd.kits[nkit] = strdup(kit);
        sd.kstrs[nkit] = (char *)malloc(nsv);
        memcpy(sd.kstrs[nkit], kstr, nsv);
        sd.nsvs[nkit] = (short) nsv;
        //printf("%d %s\n", nkit, sd.kits[nkit]);
        nkit++;
    }
    fclose(ifp);
    for (j = 0; j < nkit; j++) {
        sd.kstrs[j][11] -= sd.kstrs[j][9]; // adjust 
    }
    sd.nkit = nkit;
    sd.nlab = get_strlab(sd.labs);
    sd.tnmv  = find_misval(sd, sd.nmvs);
    return (sd);
}

int
find_kit(char *kit, char *kits[], int nkit)
{
    int j;

    for (j = 0; j < nkit; j++) {
        if (strcmp(kit, kits[j]) == 0) break;
    }
    return(j);
}

int
find_str(char *str, char *strs[], int nstr)
{
    int j;

    for (j = 0; j < nstr; j++) {
        if (strcmp(str, strs[j]) == 0) break;
    }
    return(j);
}

struct STRDAT
fix_misval(struct STRDAT sd, char force)
{
    char *s, *fs, line[MAXLIN], kit[20], str[20];
    int i, j, v0, v1;
    FILE *mfp;
    static char mfn[] = "misval.txt";

    mfp = fopen(mfn, "r");
    while (fgets(line, MAXLIN, mfp)) {
        s = line;
        while (*s == ' ') s++;
        i = 0;
        while (*s > ' ') kit[i++] = *s++;
        kit[i] = 0;
        j = find_kit(kit, sd.kits, sd.nkit);
        //printf("%d %s %s\n", j, kit, sd.kits[j]);
        while (*s++ >= ' ') {
            while (*s == ' ') s++;
            i = 0;
            while ((*s > ' ') && (*s != '=')) str[i++] = *s++;
            str[i] = 0;
            i = find_str(str, sd.labs, sd.nlab);
            if ((j < sd.nkit) && (i < sd.nsvs[j]) && (*s == '=')) {
                v0 = sd.kstrs[j][i]; // old value
                v1 = atoi(++s);      // new value
                if (!v0 || force) { // only change missing value unless forced
                    sd.kstrs[j][i] = v1;
                }
                if (v0 && (v0 != v1)) {
                    fs = force ? "forced" : "not forced";
                    printf("%s=%d->%d (%s)\n", str, v0, v1, fs);
                }
                while (*s > ' ') s++;
            }
        }
    }
    fclose(mfp);
    return(sd);
}

void
put_strdat(char *ofn, struct STRDAT sd)
{
    int i;
    FILE *ofp;

    ofp = fopen(ofn, "w"); // output file
    fputs(sd.top, ofp);
    for (i = 0; i < sd.nkit; i++) {
        put_kit(ofp, sd.kits[i], sd.kstrs[i], (int)sd.nsvs[i]);
    }
    fclose(ofp);
}

void
free_strdat(struct STRDAT sd)
{
    int i;

    for (i = 0; i < sd.nkit; i++) {
        free(sd.kits[i]);
        free(sd.kstrs[i]);
    }
    for (i = 0; i < sd.nlab; i++) {
        free(sd.labs[i]);
    }
}

int
find_diff(struct STRDAT sd1,  struct STRDAT sd2, FILE *fp)
{
    char kit[20], lab[20];
    int i, j, nfnd, v1, v2, nkit1, nkit2, nsvs, ndif;

    nkit1 = sd1.nkit;
    nkit2 = sd2.nkit;
    if (nkit1 != nkit2) {
        fprintf(stderr, "different number of kits: %d!=%d\n", nkit1, nkit2);
        exit(1);
    }
    ndif = 0;
    for (j = 0; j < nkit1; j++) {
        nsvs = sd1.nsvs[j];
        nfnd = 0;
        for (i = 0; i < nsvs; i++) {
            v1 = sd1.kstrs[j][i];
            v2 = sd2.kstrs[j][i];
            if (v1 != v2) {
                if (!nfnd) {
                    strcpy(kit, sd1.kits[j]);
                    fputs(kit, fp);
                    nfnd++;
                }
                strcpy(lab, sd1.labs[i]);
                if (v1) {
                    fprintf(fp, " %s=%d->%d", lab, v1, v2);
                } else {
                    fprintf(fp, " %s=%d", lab, v2);
                }
                ndif++;
            }
        }
        if (nfnd) fputs("\n", fp);
    }
    return (ndif);
}

int
main(int argc, char **argv)
{
    char ifn1[MAXFNL], ifn2[MAXFNL];
    int ndif;
    struct STRDAT sd1, sd2;

    if (argc < 3) {
        printf("usage: strdif filename1 filename2\n");
	exit(1);
    }
    // input STRDATA files
    strncpy(ifn1, argv[1], MAXFNL);
    sd1 = get_strdat(ifn1);
    strncpy(ifn2, argv[2], MAXFNL);
    sd2 = get_strdat(ifn2);
    // find differences
    ndif = find_diff(sd1, sd2, stdout);
    if (!ndif) printf("no diffferences\n");
    // clean up
    free_strdat(sd1);
    free_strdat(sd2);
}
