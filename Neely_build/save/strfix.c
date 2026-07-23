// strfix - specify missing values in STRDATA file

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>

#define BIGY    500
#define NOPT    100
#define MAXKIT   40
#define MAXHDR   80
#define MAXFNL  256
#define MAXSTR  840
#define MAXLIN 8192

struct STRDAT {
    int nkit, nstr, nlab, tnmv;
    short nsvs[MAXKIT], nmvs[MAXKIT];
    char *kits[MAXKIT], *kstrs[MAXKIT], *ostrs[MAXKIT];
    char top[MAXHDR], *labs[MAXSTR], *modal;
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

int
get_kit(FILE *ifp, char *hdr, char *kit, char *kstr, char *line)
{
    char *s;
    int nsv;

    strncpy(hdr, line, MAXHDR);
    fgets(line, MAXLIN, ifp);
    strncpy(kit, line, 20);
    s = kit;
    while (isalnum(*s)) s++; 
    *s = 0;
    s = line;
    while (isalnum(*s)) s++; 
    while (isspace(*s)) s++; 
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
    if (lfp == NULL) {
        printf("can't open %s\n", lfn);
        exit(1);
    }
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
find_misval(struct STRDAT sd)
{
    short *nmvs;
    int i, j, nmv, tnmv;

    nmvs = sd.nmvs;
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

int
unique(char *v, char *uv, int n)
{
    int i, j, nuv;

    nuv = 0;
    for (i = 0; i < n; i++) {
        if (!v[i]) continue; // only consider non-zero values
        for (j = 0; j < nuv; j++) {
            if (v[i] == uv[j]) break;
        }
        if (j == nuv) {
            uv[nuv] = v[i];
            nuv++;
        }
    }
    return (nuv);
}

int
modal(char *v, int n)
{
    char uv[20];
    int i, j, nuv, no, mxno, mv;

    nuv = unique(v, uv, n);
    mxno = mv = 0;
    for (j = 0; j < nuv; j++) {
        no = 0;
        for (i = 0; i < n; i++) {
            if (v[i] == uv[j]) no++;
        }
        if (mxno < no) {
            mxno = no;
            mv = uv[j];
        }
    }
    return (mv);
}

void
find_modal(struct STRDAT sd)
{
    char v[MAXKIT];
    int i, j;

    for (i = 0; i < sd.nstr; i++) {
        for (j = 0; j < sd.nkit; j++) {
            v[j] = (i < sd.nsvs[j]) ? sd.kstrs[j][i] : 0;
        }
        sd.modal[i] = modal(v, sd.nkit);
        //printf("modal[%d]=%d\n", i, sd.modal[i]);
    }
}

struct STRDAT
get_strdat(char *ifn)
{
    char hdr[MAXHDR], line[MAXLIN], kit[20], kstr[MAXSTR];
    int nkit, nstr, nsv;
    FILE *ifp;
    struct STRDAT sd;

    ifp = fopen(ifn, "r"); // input file
    if (ifp == NULL) {
        printf("can't open %s\n", ifn);
        exit(1);
    }
    fgets(sd.top, MAXHDR, ifp);
    nkit = nstr = 0;
    while (fgets(line, MAXLIN, ifp)) {
        nsv = get_kit(ifp, hdr, kit, kstr, line);
        if (nstr < nsv) nstr = nsv;
        sd.kits[nkit] = strdup(kit);
        sd.kstrs[nkit] = (char *)malloc(nsv);
        sd.ostrs[nkit] = (char *)malloc(nsv);
        memcpy(sd.kstrs[nkit], kstr, nsv);
        memcpy(sd.ostrs[nkit], kstr, nsv);
        sd.nsvs[nkit] = (short) nsv;
        //printf("%d %s\n", nkit, sd.kits[nkit]);
        nkit++;
    }
    fclose(ifp);
    sd.modal = (char *)malloc(nstr);
    sd.nkit = nkit;
    sd.nstr = nstr;
    sd.nlab = get_strlab(sd.labs);
    sd.tnmv  = find_misval(sd);
    find_modal(sd);
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

void
mod_misval(struct STRDAT sd)
{
    int i, j;

    for (j = 0; j < sd.nkit; j++) {
        for (i = 0; i < sd.nsvs[j]; i++) {
            if (sd.kstrs[j][i] == 0) {
                sd.kstrs[j][i] = sd.modal[i];
            }
        }
    }
}

int
gendst(int j1, int j2, struct STRDAT sd, int *ii, int *jj, int nkit, int nstr)
{
    int i, v1, v2, sum;

    sum = 0;
    for (i = 0; i < nstr; i++) {
        if (ii[i] >= sd.nsvs[jj[j1]]) continue;
        if (ii[i] >= sd.nsvs[jj[j2]]) continue;
        v1 = sd.kstrs[jj[j1]][ii[i]];
        v2 = sd.kstrs[jj[j2]][ii[i]];
        if (v1 && v2) {
            sum += (v1 != v2); // absolute allelle ?
        }
    }
    //printf("gendst: %d %d %d\n", j1, j2, sum);
    return sum;
}

float
avgmin(struct STRDAT sd)
{
    char v[MAXKIT], uv[20];
    float avmngd;
    int ii[MAXSTR], jj[MAXKIT];
    int i, j, k, gd, mngd, nkit, nstr, nuv;

    // select BigY kits
    nkit = 0;
    for (j = 0; j < sd.nkit; j++) {
        if (sd.nsvs[j] > BIGY) {
            jj[nkit] = j;
            nkit++;
        }
    }
    // select mutated STRs
    nstr = 0;
    for (i = 0; i < sd.nstr; i++) {
        for (k = 0; k < nkit; k++) {
            v[k] = (i < sd.nsvs[jj[k]]) ? sd.kstrs[jj[k]][i] : 0;
        }
        nuv = unique(v, uv, nkit);
        if (nuv > 1) {
            ii[nstr] = i;
            nstr++;
        }
    }
    //printf("avgmin: nkit=%d nstr=%d\n", nkit, nstr);
    // compute pairwise genetic distances
    avmngd = 0;
    for (j = 0; j < nkit; j++) {
        mngd = 999;
        for (i = 0; i < nkit; i++) {
            if (i == j) continue;
            gd = gendst(i, j, sd, ii, jj, nkit, nstr);
            if (mngd > gd) {
                mngd = gd;
            }
        }
        avmngd += mngd;
    }
    avmngd /= nkit;
    return (avmngd);
}

void
opt_misval(struct STRDAT sd)
{
    char v[MAXKIT], uv[20];
    float av, pav, avmngd;
    int i, j, k, o, nkit, nuv, nsvs, v1, v2;

    avmngd = avgmin(sd);
    printf("avmngd=%.2f", avmngd);
    nkit = sd.nkit;
    pav = avmngd;
    for (o = 0; o < NOPT; o++) {
        for (j = 0; j < nkit; j++) {
            if (sd.nsvs[j] < BIGY) continue;
            nsvs = sd.nsvs[j];
            for (i = 0; i < nsvs; i++) {
                for (k = 0; k < nkit; k++) {
                    v[k] = 0;
                    if (sd.nsvs[k] < BIGY) continue;
                    if (sd.nsvs[k] <= i) continue;
                    v[k] = sd.kstrs[k][i];
                }
                nuv = unique(v, uv, nkit);
                v1 = sd.ostrs[j][i];
                v2 = sd.kstrs[j][i];
                if ((nuv > 1) && !v1 && v2) {
                    for (k = 0; k < nuv; k++) {
                        sd.kstrs[j][i] = uv[k];
                        av = avgmin(sd);
                        if (avmngd > av) {
                            avmngd = av;
                            v2 = uv[k];
                        }
                    }
                    sd.kstrs[j][i] = v2;
                }
            }
        }
        if (pav == avmngd) break;
        pav = avmngd;
    }
    printf("->%.2f o=%d\n", avmngd, o);
}

int
find_diff(struct STRDAT sd,  char *kstrs[], char *ofn)
{
    char kit[20], lab[20], v[MAXKIT], uv[20];
    int i, j, k, nkit, nfnd, nuv, nsvs, ndif, v1, v2;
    FILE *ofp;

    ofp = fopen(ofn, "w");
    nkit = sd.nkit;
    ndif = 0;
    for (j = 0; j < nkit; j++) {
        nsvs = sd.nsvs[j];
        nfnd = 0;
        for (i = 0; i < nsvs; i++) {
            for (k = 0; k < nkit; k++) {
                v[k] = (i < sd.nsvs[k]) ? sd.kstrs[k][i] : 0;
            }
            nuv = unique(v, uv, nkit);
            v1 = kstrs[j][i];
            v2 = sd.kstrs[j][i];
            if ((nuv > 1) && !v1 && v2) {
                if (!nfnd) {
                    strcpy(kit, sd.kits[j]);
                    fputs(kit, ofp);
                    nfnd++;
                }
                strcpy(lab, sd.labs[i]);
                fprintf(ofp, " %s=%d", lab, v2);
                ndif++;
            }
        }
        if (nfnd) fputs("\n", ofp);
    }
    fclose(ofp);
    return (ndif);
}

void
fix_misval(struct STRDAT sd, int group, int force)
{
    char *s, *fs, line[MAXLIN], kit[20], str[20], mfn[MAXFNL];
    int i, j, nrep, v0, v1;
    FILE *mfp;

    sprintf(mfn, "misval%d.txt", group);
    mfp = fopen(mfn, "r");
    if (mfp == NULL) return;
    printf("replace: %s\n", mfn);
    fs = force ? "=>" : "->"; // indicator of forced replacement
    while (fgets(line, MAXLIN, mfp)) {
        s = line;
        while (isspace(*s)) s++;
        i = 0;
        while (isalnum(*s)) kit[i++] = *s++;
        kit[i] = 0;
        j = find_kit(kit, sd.kits, sd.nkit);
        //printf("%d %s %s\n", j, kit, sd.kits[j]);
        nrep = 0;
        while (*s++ >= ' ') {
            while (isspace(*s)) s++;
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
                    if (!nrep) printf("%s", kit);
                    printf(" %s=%d%s%d", str, v0, fs, v1);
                    nrep++;
                }
                while (isalnum(*s)) s++;
            }
        }
        if (nrep) printf("\n");
    }
    fclose(mfp);
}

void
rpl_misval(struct STRDAT sd, int group, int modrep, int optrep)
{
    char ofn[MAXFNL];
    int ndif;

    if (!modrep && !optrep) return;
    // replacement of missing values
    mod_misval(sd);
    if (optrep) {
        fix_misval(sd, group, 1);
        opt_misval(sd);
    }
    // output recommended replacements
    sprintf(ofn, "%srep%d.txt", optrep ? "opt" : "mod", group);
    ndif = find_diff(sd, sd.ostrs, ofn);
    printf(" output: %s ndif=%d\n", ofn, ndif);
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
        free(sd.ostrs[i]);
    }
    for (i = 0; i < sd.nlab; i++) {
        free(sd.labs[i]);
    }
    free(sd.modal);
}

void
str_adjust(struct STRDAT sd)
{
    int j;

    for (j = 0; j < sd.nkit; j++) {
        sd.kstrs[j][11] -= sd.kstrs[j][9]; // adjust 
    }
}

void
str_all_csv(struct STRDAT sd, int group)
{
    char ofn[MAXFNL];
    int i, j, nlab, val;
    FILE *ofp;

    if (group) {
        sprintf(ofn, "sv_all%d.csv", group);
    } else {
        return;
    }
    nlab = sd.nlab;
    ofp = fopen(ofn, "w");
    fputs("kit,", ofp);
    for (i = 0; i < nlab; i++) {
        fputs(sd.labs[i], ofp);
        if (i < (nlab - 1)) {
            fputs(",", ofp);
        } else {
            fputs("\n", ofp);
        }
    }
    for (j = 0; j < sd.nkit; j++) {
        fprintf(ofp,"%s,", sd.kits[j]);
        for (i = 0; i < nlab; i++) {
            if (i < sd.nsvs[j]) {
                val = sd.kstrs[j][i];
                fprintf(ofp, "%d", val);
            } else {
                fputs("0", ofp);
            }
            if (i < (nlab - 1)) {
                fputs(",", ofp);
            }
        }
        fputs("\n", ofp);
    }
    fclose(ofp);
}

void
select_subset(struct STRDAT sd, char *selkit, int group, int subset)
{
    char line[MAXLIN], select[MAXKIT], tfn[MAXFNL], lab[20], kit[20];
    int i, j, k, n, nkit;
    FILE *tfp;

    if (!group || !subset) return;
    sprintf(tfn, "subset%d.txt", group);
    tfp = fopen(tfn, "r");
    if (tfp == NULL) {
        printf("can't open %s\n", tfn);
        return;
    }
    while (subset--) fgets(line, MAXLIN, tfp);
    fclose(tfp);
    nkit = sd.nkit;
    for (k = 0; k < nkit; k++) select[k] = 0; // clear selections
    n = strlen(line);
    k = 0;                      // start at the begining
    while (isspace(line[k])) k++; // skip space
    while (isalnum(line[k])) {
        lab[k] = line[k];       // copy subset label
        k++;
    }
    lab[k] = 0;
    printf("subset: %s %s\n", tfn, lab);
    while (line[k] && (line[k] != '(')) k++;  // advance to (
    for (i = k + 1; i < n; i++) {
        while (isspace(line[i])) i++;
        if (line[i] == ')') break;
        j = 0;
        while ((line[i] > ' ') && (line[i] != ')')) kit[j++] = line[i++];
        kit[j] = 0;
        for (k = 0; k < nkit; k++) {
            if (strcmp(sd.kits[k], kit) == 0) break;
        }
        if (k == nkit) continue;
        select[k] = 1;
    }
    for (k = 0; k < nkit; k++) {
        if (!select[k]) selkit[k] = 0;
    }
}

void
select_mutated_strs(struct STRDAT sd, char *selkit, char *selstr, int nlab)
{
    int i, j, mnv, mxv, val;

    for (i = 0; i < nlab; i++) {
      mnv = 99;
      mxv = 0;
      for (j = 0; j < sd.nkit; j++) {
          val = sd.kstrs[j][i];
          if (val && selkit[j]) {
              if (mnv > val) {
                  mnv = val;
              }
              if (mxv < val) {
                  mxv = val;
              }
          }
      }
      selstr[i] = mxv ? (mnv != mxv) : 0;
    }
}

float
str_set_csv(struct STRDAT sd, int group, int subset)
{
    char selkit[MAXKIT], selstr[MAXSTR];
    char ofn[MAXFNL], dfn[MAXFNL];
    float avmngd;
    int ii[MAXSTR], jj[MAXKIT], gd[MAXKIT][MAXKIT];
    int i, j, nlab, val, mngd, mnns, nkit, nstr, tnmv;
    FILE *ofp, *dfp;
    static char ss[2] = "a"; // subset label

    if (group) {
        sprintf(ofn, "sv_set%d.csv", group);
        sprintf(dfn, "gd_set%d.csv", group);
    } else {
        strcpy(ofn, "/dev/null");
        strcpy(dfn, "/dev/null");
    }
    nlab = sd.nlab;
    // select kits & STRs
    for (j = 0; j < sd.nkit; j++) {
        selkit[j] = (sd.nsvs[j] > BIGY);
    }
    select_subset(sd, selkit, group, subset);
    select_mutated_strs(sd, selkit, selstr, nlab);
    mnns = MAXSTR;
    nkit = nstr = 0;
    for (j = 0; j < sd.nkit; j++) {
        if (selkit[j]) {
            jj[nkit] = j;
            if (mnns > sd.nsvs[j]) {
                mnns = sd.nsvs[j];
            }
            nkit++;
        }
    }
    for (i = 0; i < mnns; i++) {
        if (selstr[i]) {
            ii[nstr] = i;
            nstr++;
        }
    }
    if (nstr == 0) {
        printf("empty set: nkit=%d mnns=%d\n", nkit, mnns);
        return (0);
    }
    tnmv = 0;
    for (j = 0; j < nkit; j++) {
        for (i = 0; i < nstr; i++) {
            if (ii[i] >= sd.nsvs[jj[j]]) continue;
            if (sd.kstrs[jj[j]][ii[i]] == 0) tnmv++;
        }
    }
    ss[0] = subset ? ('a' + subset - 1) : 0;
    printf(" select: nkit=%d nstr=%d group=%d%s tnmv=%d ",
        nkit, nstr, group, ss, tnmv);
    if (!nkit || !nstr) return(0);
    // output CSV
    ofp = fopen(ofn, "w");
    fputs("kit,", ofp);
    for (i = 0; i < mnns; i++) {
        if (!selstr[i]) continue;
        fputs(sd.labs[i], ofp);
        if (i < (mnns - 1)) {
            fputs(",", ofp);
        }
    }
    fputs("\n", ofp);
    for (j = 0; j < sd.nkit; j++) {
        if (!selkit[j]) continue;
        fprintf(ofp,"%s,", sd.kits[j]);
        for (i = 0; i < mnns; i++) {
            if (!selstr[i]) continue;
            if (i < sd.nsvs[j]) {
                val = sd.kstrs[j][i];
                fprintf(ofp, "%d", val);
            } else {
                fputs("0", ofp);
            }
            if (i < (mnns - 1)) {
                fputs(",", ofp);
            }
        }
        fputs("\n", ofp);
    }
    fclose(ofp);
    // compute pairwise genetic distances
    dfp = fopen(dfn, "w");
    fputs("kit,", dfp);
    for (j = 0; j < nkit; j++) {
        fprintf(dfp, "%s,", sd.kits[jj[j]]);
    }
    fputs("mngd\n", dfp);
    avmngd = 0;
    for (j = 0; j < nkit; j++) {
        fprintf(dfp, "%s,", sd.kits[jj[j]]);
        mngd = 999;
        for (i = 0; i < nkit; i++) {
            gd[i][j] = gendst(i, j, sd, ii, jj, nkit, nstr);
            fprintf(dfp, "%d,",gd[i][j]);
            if ((i != j) && (mngd > gd[i][j])) {
                mngd = gd[i][j];
            }
        }
        fprintf(ofp, "%d\n", mngd);
        avmngd += mngd;
    }
    fclose(dfp);
    avmngd /= nkit;
    printf("avmngd=%.1f\n", avmngd);
    return (avmngd);
}

int
main(int argc, char **argv)
{
    char ss[2], ifn[MAXFNL];
    int force, modrep, optrep, allcsv, group, subset;
    struct STRDAT sd;
    static char ofn[] = "strfix.txt";

    if (argc < 2) {
        printf("usage: strfix [option] filename [group]\n");
        printf("options: \n");
        printf("    -a   write all kits & STRs to CSV\n");
        printf("    -f   force replacement of missing values\n");
        printf("    -m   modal replacement of missing values\n");
        printf("    -o   optimal replacement of missing values\n");
	exit(1);
    }
    force = modrep = optrep = allcsv = subset = 0;
    if (argc > 2) {
        while (argv[1][0] == '-') {
            if (argv[1][1] == 'a') {
                allcsv++;
            } else if (argv[1][1] == 'f') {
                force++;
            } else if (argv[1][1] == 'm') {
                modrep++;
            } else if (argv[1][1] == 'o') {
                optrep++;
            }
            argv++;
            argc--;
        }
    }
    if (argc > 2) {
        group = atoi(argv[2]);
        strncpy(ss, argv[2] + 1, 1);
        if ((ss[0] >= 'a') && (ss[0] <= 'z')) {
            subset = ss[0] - 'a' + 1;
        } 
    } else {
        group = 0;
    }
    // input STRDATA
    strncpy(ifn, argv[1], MAXFNL);
    sd = get_strdat(ifn);
    printf("  input: %s nkit=%d nstr=%d tnmv=%d group=%d%s\n",
        ifn, sd.nkit, sd.nstr, sd.tnmv, group, ss);
    // replace missing values
    rpl_misval(sd, group, modrep, optrep);
    fix_misval(sd, group, force);
    // output STRDATA
    put_strdat(ofn, sd);
    // output CSV
    str_adjust(sd);
    if (allcsv) str_all_csv(sd, group);
    str_set_csv(sd, group, subset);
    // clean up
    free_strdat(sd);
}
