// snpyear - compile list of SNP years

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "tremer.h"

static char ifn[MAXFNL] = "snpy.txt";
static char ofn[MAXFNL] = "snps.txt";
static int debug = 0;

static int
getval(char *line, char *str)
{
    char *s;
    int nsv = 0;

    nsv = 0;
    s = line;
    // create list of STR values
    while (*s && (nsv < MAXSTR)) {
        while(isspace(*s)) s++;
        if (*s == 'N') {
            str[nsv++] = 0;
        } else if ((*s >= '1') && (*s <= '9')) {
            str[nsv++] = (char) atoi(s);
        }
        while(*s && !isspace(*s)) s++;
    }
    return (nsv);
}

static int
get_kit(FILE *ifp, char *hdr, char *kit, char *kstr, char *line)
{
    char *s;
    int nsv;

    strncpy(hdr, line, MAXHDR);
    fgets(line, MAXLIN, ifp);
    strncpy(kit, line, 80);
    s = kit;
    while (isalnum(*s)) s++; 
    *s = 0;
    s = line;
    while (isalnum(*s)) s++; 
    while (isspace(*s)) s++; 
    nsv = getval(s, kstr);
    return (nsv);
}

static void
add_kit(struct KITDAT *kdp, char *name, int year, char *sv, int nsv)
{
    int nkit, nstr;

    nkit = kdp->nkit;
    nstr = kdp->nstr;
    if (nstr < nsv) nstr = nsv;
    kdp->kits[nkit] = (char *)malloc(40);
    strncpy(kdp->kits[nkit], name, 40);
    kdp->kstrs[nkit] = (char *)malloc(MAXSTR);
    kdp->ostrs[nkit] = (char *)malloc(MAXSTR);
    kdp->nsvs[nkit] = (short) nsv;
    kdp->kity[nkit] = (short) year;
    kdp->nlft[0] = 0;
    memcpy(kdp->kstrs[nkit], sv, nsv);
    memcpy(kdp->ostrs[nkit], sv, nsv);
    kdp->nkit = ++nkit;
    kdp->nstr = nstr;
}

static void
get_strdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, hdr[MAXHDR], kstr[MAXSTR], kit[80];
    int nsv;

    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        nsv = get_kit(ifp, hdr, kit, kstr, line);
        add_kit(kdp, kit, 0, kstr, nsv);
        //printf("%8s %2d %2d %2d\n", kit, kstr[32], kstr[33], kstr[34]);
    }
    if (s == NULL) line[0] = EOF;
}

static int
find_kit(struct KITDAT *kdp, char *kit)
{
    int j, n;

    n = strlen(kit);
    for (j = 0; j < kdp->nkit; j++) {
         if (strncmp(kit, kdp->kits[j], n) == 0) break;
    }
    return(j);
}

static void
get_snpdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kit[80], snpnam[80];
    int i, j, n, ikit, isnp, nkit, nsnp;

    nkit = kdp->nkit;
    nsnp = kdp->nsnp;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        j = 0;
        while (isalnum(line[i])) kit[j++] = line[i++];
        if (j == 0) continue;
        kit[j] = 0;
        ikit = find_kit(kdp, kit);
        if (ikit == nkit) {
            printf("WARNING: kit %s has no STR\n", kit);
            continue;
        }
        while ((i < n) && (line[i] != '(')) i++;
        while ((i < n) && (line[i] != ')')) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) snpnam[j++] = line[i++];
            if (j) {
                snpnam[j] = 0;
                for (isnp = 0; isnp < nsnp; isnp++) {
                    if (strcmp(snpnam, kdp->snps[isnp]) == 0) break;
                }
                if (isnp == nsnp) {
                    kdp->snps[nsnp] = strdup(snpnam);
                    nsnp++;
                }
                if (line[i] == '+') kdp->ksnp[ikit][isnp] = 1;
                if (line[i] == '-') kdp->ksnp[ikit][isnp] = 2;
                //printf("%8s %8s %c\n", kit, snpnam, line[i]);
            }
            while ((i < n) && !isalnum(line[i]) && (line[i] != ')')) i++;
        }
    }
    if (s == NULL) line[0] = EOF;
    kdp->nsnp = nsnp;
}

void
fetch_snptree(struct KITDAT *kdp)
{
    char snpnam[80], snppar[80], line[MAXLIN];
    int i, j, k, y;
    FILE *ifp;

    ifp = fopen(ifn, "r");
    while (fgets(line, MAXLIN, ifp)) {
       i = 0;
       while (isspace(line[i])) i++;
       j = 0;
       while (isalnum(line[i])) snpnam[j++] = line[i++];
       snpnam[j] = 0;
       while (isspace(line[i])) i++;
       y = atoi(line + i);
       while (isdigit(line[i])) i++;
       while (isspace(line[i])) i++;
       j = 0;
       while (isalnum(line[i])) snppar[j++] = line[i++];
       snppar[j] = 0;
       for (k = 0; k < kdp->nsnp; k++) {
           if (debug && 0) {
               printf("%8s %4d : %8s %2d %2d\n", snpnam, y,
                   kdp->snps[k], k, kdp->nsnp);
           }
           if (strcmp(snpnam, kdp->snps[k]) == 0) {
               kdp->snpy[k] = y;
               kdp->snpp[k] = strdup(snppar);
               //printf("%9s  %4d  %s\n", snpnam, y, snppar);
               break;
           }
       }
    }
    fclose(ifp);
}

static struct KITDAT *
get_sapp(char *ifn, int group)
{
    char *s, line[MAXLIN];
    FILE *ifp;
    static struct KITDAT kd;

    kd.nkit = kd.nstr = kd.nsnp = 0;
    ifp = fopen(ifn, "r"); // input file
    if (ifp == NULL) {
        printf("ERROR: can't open %s\n", ifn);
        exit(1);
    }
    while (fgets(line, MAXLIN, ifp)) {
        if (line[0] == '/') break;
    }
    while (line[0] == '/') {
        if (strncmp(line, "/STRDATA", 8) == 0) {
            get_strdata(&kd, ifp, line);
        } else if (strncmp(line, "/SNPDATA", 8) == 0) {
            get_snpdata(&kd, ifp, line);
        } else {
            while ((s = fgets(line, MAXLIN, ifp))) {
                if (line[0] == '/') break;
            }
            if (s == NULL) line[0] = EOF;
        }
    }
    fclose(ifp);
    // get SNP year & parent
    fetch_snptree(&kd);
    // done
    return (&kd);
}

static void
free_kitdat(struct KITDAT *kdp)
{
    int i, n;

    n = kdp->nkit + kdp->ngen + kdp->nsnp + kdp->nsbg;
    for (i = 0; i < n; i++) {
        free_null(kdp->kits[i]);
        free_null(kdp->kstrs[i]);
        free_null(kdp->ostrs[i]);
    }
    for (i = 0; i < kdp->nlab; i++) {
        free_null(kdp->labs[i]);
    }
    for (i = 0; i < kdp->ngen; i++) {
        free_null(kdp->gens[i]);
    }
    for (i = 0; i < kdp->nsnp; i++) {
        free_null(kdp->snps[i]);
    }
    for (i = 0; i < kdp->nsbg; i++) {
        free_null(kdp->sbgs[i]);
    }
}

static int
group_number(char *nfn)
{
   while (*nfn && !isdigit(*nfn)) nfn++;
   return (atoi(nfn));
}

int
add_snps(struct KITDAT *kdp, char *snps[], char *snpp[], 
    int *snpy, int nsnps)
{
    int i, nsnp;

    nsnp = kdp->nsnp;
    for (i = 0; i < nsnp; i++) {
        snps[nsnps + i] = strdup(kdp->snps[i]);
        snpy[nsnps + i] = kdp->snpy[i];
        snpp[nsnps + i] = kdp->snpp[i];
    }
    return (nsnp);
}

int
sort_snps(struct KITDAT *kdp, char *snps[], char *snpp[], 
    int *snpy, int nsnps)
{
    char *s, *p;
    int i, j, y;

    // sort SNPs
    for (i = 0; i < nsnps; i++) {
        for (j = i; j < nsnps; j++) {
            if (strcmp(snps[i], snps[j]) > 0) {
                s = snps[i];
                p = snpp[i];
                y = snpy[i];
                snps[i] = snps[j];
                snpp[i] = snpp[j];
                snpy[i] = snpy[j];
                snps[j] = s;
                snpp[j] = p;
                snpy[j] = y;
            }
        }
    }
    if (nsnps) { // remove duplicate SNPs
        j = 0;
        for (i = 1; i < nsnps; i++) {
            if (strcmp(snps[i], snps[j]) > 0) j++;
            snps[j] = snps[i];
            snpp[j] = snpp[i];
            snpy[j] = snpy[i];
        }
        nsnps = j + 1;
    }
    return (nsnps);
}

void
report_snps(struct KITDAT *kdp, char *snps[], char *snpp[], 
    int *snpy, int nsnps, char *ofn)
{
    int i;
    FILE *ofp;

    ofp = fopen(ofn, "w");
    for (i = 0; i < nsnps; i++) {
        fprintf(ofp, "%-9s %4d  %s\n", snps[i], snpy[i], snpp[i]);
    }
    fclose(ofp);
}

static void
report_snptree(struct KITDAT *kdp, int group)
{
    char ofn[MAXFNL];
    int i, nsnp;
    FILE *ofp;
    static char *name = "snptree";

    if (group) {
        sprintf(ofn, "%s%d.txt", name, group);
    } else {
        sprintf(ofn, "%s.txt", name);
    }
    nsnp = kdp->nsnp;
    ofp = fopen(ofn, "w");
    for (i = 0; i < nsnp; i++) {
        fprintf(ofp, "%-9s %4d  %s\n",
            kdp->snps[i], kdp->snpy[i], kdp->snpp[i]);
    }
    fclose(ofp);
}

int
main(int argc, char **argv)
{
    char sfn[MAXFNL], *snps[MAXSNP], *snpp[MAXSNP];
    int nsnps, snpy[MAXSNP];
    struct KITDAT *kdp;
    static int group = 0;
    static int force = 0;

    if (argc < 2) {
        printf("usage: snpyear filename ...\n");
	exit(1);
    }
    strcpy(ofn, "snpys.txt");
    nsnps = 0;
    while (argc > 1) {
        if (argv[1][0] == '-') {
            if (argv[1][1] == 'f') {
                force++;
            } else if (argv[1][1] == 'i') {
                if (argc > 2) {
                    strncpy(ifn, argv[2], MAXFNL);
                    argv++;
                    argc--;
                }
            } else if (argv[1][1] == 'o') {
                if (argc > 2) {
                    strncpy(ofn, argv[2], MAXFNL);
                    argv++;
                    argc--;
                }
            }
        } else {
            strncpy(sfn, argv[1], MAXFNL);
            group = group_number(sfn);
            // input KIT data
            kdp = get_sapp(sfn, group);
            nsnps += add_snps(kdp, snps, snpp, snpy, nsnps);
            report_snptree(kdp, group);
            // clean up
            free_kitdat(kdp);
        }
        argv++;
        argc--;
    }
    nsnps = sort_snps(kdp, snps, snpp, snpy, nsnps);
    if (force) {
        report_snps(kdp, snps, snpp, snpy, nsnps, ofn);
    } else {
        *ofn = 0;
    }
    printf("snpyear: nsnps=%d ", nsnps);
    printf("ifn=%s ofn=%s\n", ifn, ofn);
} 

