// grpsum - summarize group features

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>

#define BIGY    500
#define NOPT    100
#define MAXKIT  160
#define MAXSNP  160
#define MAXGEN  160
#define MAXHDR  160
#define MAXSIG  160
#define MAXFNL  256
#define MAXSTR  900
#define MAXLIN 8192

struct SIGDAT {
    int n, s[MAXSIG], p[MAXSIG], m[MAXSIG];
};
struct KITDAT {
    int nkit, nstr, nsnp, ngen, nlab, tnmv;
    short nsvs[MAXKIT], nmvs[MAXKIT];
    short snpy[MAXSNP], geny[MAXGEN];
    char *kits[MAXKIT], *kstrs[MAXKIT], *ostrs[MAXKIT];
    char *snps[MAXSNP], *gens[MAXGEN], pars[MAXSNP];
    char ksnp[MAXKIT][MAXSNP], kgen[MAXKIT][MAXGEN], gsnp[MAXGEN][MAXSNP];
    char top[MAXHDR], *labs[MAXSTR];
    int nmem[MAXGEN], nkid[MAXGEN], nlft[MAXGEN];
    char parent[MAXGEN], kids[MAXGEN][20];
    char mems[MAXGEN][MAXGEN], lfts[MAXGEN][MAXGEN];
    struct SIGDAT sig[MAXGEN];
};

static int debug = 0;
static int lefts = 0;
static int group = 0;
static int modal = 0;
static short so[MAXSTR];

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
find_misval(struct KITDAT *kdp)
{
    short *nmvs;
    int i, j, nmv, tnmv;

    nmvs = kdp->nmvs;
    tnmv = 0;
    for (i = 0; i < kdp->nkit; i++) {
        nmv = 0;
        for (j = 0; j < kdp->nsvs[i]; j++) {
            if (kdp->kstrs[i][j] == 0) nmv++;
        }
        nmvs[i] = nmv;
        tnmv += nmv;
    }
    return (tnmv);
}

int
unique_nonzero(char *v, char *uv, int n)
{
    int i, j, nuv;

    nuv = 0;
    for (i = 0; i < n; i++) {
        if (v[i] == 0) continue; // only consider non-zero values
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
modal_val(char *v, int n)
{
    char uv[20];
    int i, j, nuv, no, mxno, mv;

    nuv = unique_nonzero(v, uv, n);
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

int
count_val(char *v, int n, int u)
{
    int i, cnt;

    cnt = 0;
    for (i = 0; i < n; i++) {
        if (v[i] == u) cnt++;
    }
    return (cnt);
}

void
init_modal(struct KITDAT *kdp)
{
    char v[MAXKIT];
    int i, j, nsv, nkit, nstr;

    if (modal) return; // check whether already initialied
    nkit = kdp->nkit;
    nstr = kdp->nstr;
    kdp->kits[nkit] = (char *)malloc(20);
    kdp->kstrs[nkit] = (char *)malloc(nstr);
    kdp->ostrs[nkit] = (char *)malloc(nstr);
    kdp->nsvs[nkit] = (short) nstr;
    for (i = 0; i < nstr; i++) {
        for (j = 0; j < nkit; j++) {
            nsv = kdp->nsvs[j];
            v[j] = (i < nsv) ? kdp->kstrs[j][i] : 0;
        }
        kdp->kstrs[nkit][i] = modal_val(v, nkit);
    }
    kdp->nlft[0] = 0;
    modal++;
}

int
node_name(struct KITDAT *kdp, char *name, int node)
{
    int y;

    // find node name
    if (node < kdp->ngen) {
        sprintf(name, "%s", kdp->gens[node]);
        y = kdp->geny[node];
    } else {
        node -= kdp->ngen;
        sprintf(name, "%s", kdp->snps[node]);
        y = kdp->snpy[node];
    }
    if (y) sprintf(name + strlen(name), ".%d", y); // append year
    return(strlen(name));
}

int
eat_leftovers(struct KITDAT *kdp, int node, int i, int mv)
{
    char *str, *lfts, *mems, u[MAXKIT], mvk[MAXGEN], name[40];
    int nsv, nlft, nkid, kid, nmem, nmvk;
    int cmv, k, m, n, tgt, tran, mc, lc, lft, mem;

    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    lc = 0;
    for (m = 0; m < nlft; m++) {
        lft = lfts[m];
        nsv = kdp->nsvs[lft];
        str = kdp->kstrs[lft];
        u[m] = (i < nsv) ? str[i] : 0;
        if (u[m] == mv) lc++;
    }
    //return (lc); // <- DEBUG
    nkid = kdp->nkid[node];
    node_name(kdp, name, node);
    if (nkid == 0) {
        if (debug && (lc > 1) && (lc < nlft)) {
            printf("*** insert STR branch after %s: %s=%d lc=%d/%d\n",
                name, kdp->labs[i], mv, lc, nlft);
        }
        return (lc);
    }
    // fetch node name
    // find number of kids containing modal value
    nmvk = 0;
    for (k = 0; k < nkid; k++) {
        kid = kdp->kids[node][k];
        nmem = kdp->nmem[kid];
        mems = kdp->mems[kid];
        for (m = 0; m < nmem; m++) {
            mem = mems[m];
            nsv = kdp->nsvs[mem];
            str = kdp->kstrs[mem];
            u[m] = (i < nsv) ? str[i] : 0;
        }
        mc = count_val(u, nmem, mv);
        if (mc) {
            mvk[nmvk] = mc;
            nmvk++;
        }
    }
    cmv = 0;
    if (nmvk == nkid) {
        cmv = (nkid == 1) ? lc : 1;
    } else if (nmvk == 1) {
        cmv = lc;
        tgt = 0;
    } else if (nmvk > 0) {
        tgt = kdp->ngen + kdp->nsnp;
        if (debug) {
            printf("*** insert STR branch after %s: %s=%d\n",
                name, kdp->labs[i], mv);
        }
        return (1);
    }
    if (lc && !cmv) { // move leftovers to kid
        if (debug) {
            printf("*** move leftovers from %s to kid%d: %s=%d\n",
                name, nkid, kdp->labs[i], mv);
        }
        // copy mv from node to kid
        kid = kdp->kids[node][tgt];
        nmem = kdp->nmem[kid];
        mems = kdp->mems[kid];
        nlft = kdp->nlft[kid];
        lfts = kdp->lfts[kid];
        tran = 0;
        for (m = 0; m < nlft; m++) {
            n = kdp->lfts[node][m];
            if (kdp->kstrs[n][i] == mv) {
                mems[nmem + tran] = n;
                lfts[nlft + tran] = n;
                kdp->lfts[node][m] = 0; // clear copied leftover
                tran++;
            }
        }
        kdp->nmem[kid] += tran;
        kdp->nlft[kid] += tran;
        // delete mv from node
        nmem = kdp->nmem[node];
        mems = kdp->mems[node];
        nlft = kdp->nlft[node];
        lfts = kdp->lfts[node];
        tran = 0;
        for (m = 0; m < nmem; m++) {
            n = mems[m];
            if (kdp->kstrs[n][i] == mv) {
                tran++;
            } else {
                mems[m - tran] = mems[m];
            }
        }
        kdp->nmem[node] -= tran;
        tran = 0;
        for (m = 0; m < nlft; m++) {
            n = lfts[m];
            if (kdp->kstrs[n][i] == mv) {
                tran++;
            } else {
                lfts[m - tran] = lfts[m];
            }
        }
        kdp->nlft[node] -= tran;
    }
    return (cmv);
}

void
sort_order(struct KITDAT *kdp, short *so, int rank)
{
    char v[MAXKIT], uv[MAXKIT], cnt[MAXKIT], nmut[MAXSTR];
    int i, j, k, jmx, sum, nstr, nkit, nuv, temp;

    nstr = kdp->nstr;
    nkit = kdp->nkit;
    for (i = 0; i < nstr; i++) {
        so[i] = (short)i;
    }
    if (!rank) return;
    // count STR mutations
    for (i = 0; i < nstr; i++) {
        for (k = 0; k < nkit; k++) {
            v[k] = (i < kdp->nsvs[k]) ? kdp->kstrs[k][i] : 0;
        }
        nuv = unique_nonzero(v, uv, nkit);
        sum = jmx = 0;
        for (j = 0; j < nuv; j++) {
            cnt[j] = count_val(v, nkit, uv[j]);
            sum += (int)cnt[j];
            if (cnt[jmx] < cnt[j]) jmx = j;
        }
        nmut[i] = nuv ? sum - cnt[jmx] : 0; // number of mutations
    }
    // sort by number of mutations
    for (i = 0; i < nstr; i++) {
        for (j = i; j < nstr; j++) {
            if (nmut[so[i]] < nmut[so[j]]) {
                temp = so[i];
                so[i] = so[j];
                so[j] = temp;
            }
        }
    }
    //printf("%s %s\n", kdp->labs[so[0]], kdp->labs[so[1]]); 
}

void
inherit_strs(struct KITDAT *kdp, int node)
{
    char *s, *str, *mems, kit[20], par[20], v[MAXKIT];
    int nsv, nsig, nkid, ngen, nkit, nmem, nstr;
    int i, ii, b, k, m, p, mem, cmv, y, mc, mv, pc, pv;
 
    init_modal(kdp);
    nkit = kdp->nkit;
    ngen = kdp->ngen;
    nstr = kdp->nstr;
    k = node + nkit;
    p = kdp->parent[node];
    s = (p < ngen) ? kdp->gens[p] : kdp->snps[p - ngen];
    y = (p < ngen) ? kdp->geny[p] : kdp->snpy[p - ngen];
    sprintf(par, "%s.%d", s, y);
    s = (node < ngen) ? kdp->gens[node] : kdp->snps[node - ngen];
    y = (node < ngen) ? kdp->geny[node] : kdp->snpy[node - ngen];
    sprintf(kit, "%s.%d", s, y);
    nsv = kdp->nsvs[nkit];
    kdp->nsvs[k] = nsv;
    kdp->kits[k] = (char *)malloc(20);
    kdp->kstrs[k] = (char *)malloc(nsv);
    kdp->ostrs[k] = (char *)malloc(nsv);
    strcpy(kdp->kits[k], kit);
    nmem = kdp->nmem[node];
    mems = kdp->mems[node];
    for (i = 0; i < nstr; i++) {
        ii = so[i];  // sorted order
        for (m = 0; m < nmem; m++) {
            mem = mems[m];
            nsv = kdp->nsvs[mem];
            str = kdp->kstrs[mem];
            v[m] = (ii < nsv) ? str[ii] : 0;
        }
        pv = kdp->kstrs[p + nkit][ii];
        mv = modal_val(v, nmem);
        pc = count_val(v, nmem, pv);
        mc = count_val(v, nmem, mv);
        if ((mv && pv) && (mv != pv) && (mc > pc)) {
            cmv = eat_leftovers(kdp, node, ii, mv);
        } else {
            cmv = 0;
        }
        if (cmv && (mc > 1)) {
            nsig = kdp->sig[node].n;
            kdp->sig[node].s[nsig] = ii;
            kdp->sig[node].p[nsig] = pv;
            kdp->sig[node].m[nsig] = mv;
            kdp->sig[node].n++;
        }
        kdp->kstrs[k][ii] = cmv ? mv : pv;
    }
    if (debug) {
        nsig = kdp->sig[node].n;
        if (nsig) {
            printf("%2d %17s", node, kit);
            for (b = 0; b < nsig; b++) {
                str = kdp->labs[kdp->sig[node].s[b]];
                pv = kdp->sig[node].p[b];
                mv = kdp->sig[node].m[b];
                printf(" %s=%d->%d", str, pv, mv);
            }
            printf("\n");
        }
    }
    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        inherit_strs(kdp, kdp->kids[node][k]);
    }
}

void
get_strdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, hdr[MAXHDR], kstr[MAXSTR], kit[20];
    int nchr, nkit, nstr, nsv;

    nkit = kdp->nkit;
    nstr = kdp->nstr;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        nsv = get_kit(ifp, hdr, kit, kstr, line);
        if (nstr < nsv) nstr = nsv;
        nchr = strlen(kit);
        kdp->kits[nkit] = (char *)malloc(nchr + 1);
        kdp->kstrs[nkit] = (char *)malloc(nsv);
        kdp->ostrs[nkit] = (char *)malloc(nsv);
        strcpy(kdp->kits[nkit], kit);
        memcpy(kdp->kstrs[nkit], kstr, nsv);
        memcpy(kdp->ostrs[nkit], kstr, nsv);
        kdp->nsvs[nkit] = (short) nsv;
        //printf("%d %s\n", nkit, kdp->kits[nkit]);
        nkit++;
    }
    if (s == NULL) line[0] = EOF;
    kdp->nkit = nkit;
    kdp->nstr = nstr;
}

void
get_modal(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kstr[MAXSTR];
    int nkit, nsv;

    nkit = kdp->nkit;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        nsv = getval(line, kstr);
        kdp->kits[nkit] = (char *)malloc(20);
        kdp->kstrs[nkit] = (char *)malloc(nsv);
        kdp->ostrs[nkit] = (char *)malloc(nsv);
        strcpy(kdp->kits[nkit], "modal");
        memcpy(kdp->kstrs[nkit], kstr, nsv);
        memcpy(kdp->ostrs[nkit], kstr, nsv);
        kdp->nsvs[nkit] = (short) nsv;
    }
    if (s == NULL) line[0] = EOF;
    kdp->nlft[0] = 0;
    modal++;
    //printf("get_modal: %d %d\n", nkit, kdp->kstrs[nkit][0]);
}

void
get_snpdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kit[20], snpnam[20];
    int i, j, n, ikit, isnp, nchr, nkit, nsnp;

    nkit = kdp->nkit;
    nsnp = kdp->nsnp;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        j = 0;
        while (isalnum(line[i])) kit[j++] = line[i++];
        kit[j] = 0;
        for (ikit = 0; ikit < nkit; ikit++) {
             if (strcmp(kit, kdp->kits[ikit]) == 0) break;
        }
        while ((i < n) && (line[i] != '(')) i++;
        while ((i < n) && (line[i] != ')')) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) snpnam[j++] = line[i++];
            if (j) {
                snpnam[j] = 0;
                nchr = strlen(snpnam);
                for (isnp = 0; isnp < nsnp; isnp++) {
                    if (strcmp(snpnam, kdp->snps[isnp]) == 0) break;
                }
                if (isnp == nsnp) {
                    kdp->snps[nsnp] = (char *)malloc(nchr + 1);
                    strncpy(kdp->snps[nsnp], snpnam, nchr);
                    kdp->pars[nsnp] = nsnp - 1;
                    nsnp++;
                }
                //printf("%2d %s %s\n", k, kit, snpnam);
                if (line[i] == '+') kdp->ksnp[ikit][isnp] = 1;
                if (line[i] == '-') kdp->ksnp[ikit][isnp] = 2;
            }
            while ((i < n) && !isalnum(line[i]) && (line[i] != ')')) i++;
        }
    }
    if (s == NULL) line[0] = EOF;
    kdp->nsnp = nsnp;
}

void
get_snptree(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, snpnam[20];
    int i, j, n, p, w, isnp, nsnp;

    nsnp = kdp->nsnp;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        p = w = 0;
        while (i < n) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) snpnam[j++] = line[i++];
            if (j) {
                snpnam[j] = 0;
                for (isnp = 0; isnp < nsnp; isnp++) {
                    if (strcmp(snpnam, kdp->snps[isnp]) == 0) break;
                }
                if (w) kdp->pars[isnp] = p;
                //printf("%9s i=%2d p=%2d\n", snpnam, isnp, p);
                p = isnp;
                w++;
            }
        }
    }
    if (s == NULL) line[0] = EOF;
}

void
get_gendata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kit[20], gennam[20];
    int i, j, n, ikit, igen, nchr, nkit, ngen;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        j = 0;
        while (isalnum(line[i])) gennam[j++] = line[i++];
        if (j == 0) continue;
        gennam[j] = 0;
        nchr = strlen(gennam);
        kdp->gens[ngen] = (char *)malloc(nchr + 1);
        strncpy(kdp->gens[ngen], gennam, nchr);
        if (line[i] == '.') {
            kdp->geny[ngen] = atoi(line + i + 1);
        }
        igen = ngen;
        ngen++;
        while ((i < n) && (line[i] != '(')) i++;
        while ((i < n) && (line[i] != ')')) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) kit[j++] = line[i++];
            if (j) {
                kit[j] = 0;
                nchr = strlen(kit);
                for (ikit = 0; ikit < nkit; ikit++) {
                     if (strcmp(kit, kdp->kits[ikit]) == 0) break;
                }
                if (ikit < nkit) kdp->kgen[ikit][igen] = 1;
                //printf("%2d %s %s : %d %d\n", ngen, gennam, kit, igen, ikit);
            }
            while ((i < n) && !isalnum(line[i]) && (line[i] != ')')) i++;
        }
    }
    if (s == NULL) line[0] = EOF;
    kdp->ngen = ngen;
}

void
get_snpy(struct KITDAT *kdp)
{
    char snpnam[20], line[MAXLIN];
    int i, j, k, y;
    static char *ifn = "snpy.txt";
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
       for (k = 0; k < kdp->nsnp; k++) {
           //printf("%s %d : %s %d %d\n", snpnam, y, kdp->snps[k], k, kdp->nsnp);
           if (strcmp(snpnam, kdp->snps[k]) == 0) {
               kdp->snpy[k] = y;
               break;
           }
       }
    }
    fclose(ifp);
}

void
init_gendat(struct KITDAT *kdp)
{
    int i;

    kdp->gens[0] = (char *)malloc(20);
    if (group) {
        sprintf(kdp->gens[0], "Group%d", group);
    } else {
        strcpy(kdp->gens[0], "Group");
    }
    for (i = 0; i < MAXKIT; i++) kdp->kgen[i][0] = 1;
    kdp->geny[0] = 0;
    kdp->ngen = 1;
}

int
memsubset(struct KITDAT *kdp, int s1, int s2)
{
    char *a, *b;
    int i, j, m, n;

    a = kdp->mems[s1];
    b = kdp->mems[s2];
    m = kdp->nmem[s1];
    n = kdp->nmem[s2];
    if (m > n) return (0);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (a[i] == b[j]) break;
        }
        if (j == n) return (0);
    }
    return (1);
}

int
pair_snps(int s1, int s2)
{
    static char c[5][5] = {{0, 3, 4, 3, 4},
                           {3, 1, 2, 1, 2},
                           {4, 2, 2, 4, 4},
                           {3, 1, 4, 4, 4},
                           {4, 2, 4, 4, 4}};
    return (c[s1][s2]);
}

int
combine_snps(struct KITDAT *kdp,  int g, int s)
{
    int k, m, nmem, ksnp, gsnp;

    nmem = kdp->nmem[g];
    k = kdp->mems[g][0];
    ksnp = kdp->ksnp[k][s];
    gsnp = ksnp;
    for (m = 1; m < nmem; m++) {
        k = kdp->mems[g][m];
        ksnp = kdp->ksnp[k][s];
        gsnp = pair_snps(gsnp, ksnp);
    }
    return (gsnp);
}

void
kit2gen(struct KITDAT *kdp)
{
    int g, s, nsnp, ngen;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    for (g = 0; g < ngen; g++) {
        for (s = 0; s < nsnp; s++) {
            kdp->gsnp[g][s] = combine_snps(kdp, g, s);
        }
    }
}

void
gen2gen(struct KITDAT *kdp)
{
    int g, s, p, gsnp, nsnp, ngen;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    for (s = 0; s < nsnp; s++) {
        if (kdp->gsnp[0][s] > 2) kdp->gsnp[0][s] -= 2;
    }
    for (g = 1; g < ngen; g++) {
        p = kdp->parent[g];
        for (s = 0; s < nsnp; s++) {
            gsnp = kdp->gsnp[g][s];
            if (gsnp == 0) kdp->gsnp[g][s] = kdp->gsnp[p][s];
            if (gsnp >= 3) kdp->gsnp[g][s] -= 2;
        }
    }
}

void
infer_snps(struct KITDAT *kdp)
{
    kit2gen(kdp);
    gen2gen(kdp);
}

void
find_members(struct KITDAT *kdp)
{
    int ngen, nsnp, nkit, nmem;
    int i, j, n;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    // find GEN members
    for (i = 0; i < ngen; i++) {
        nmem = 0;
        for (j = 0; j < nkit; j++) {
            if (kdp->kgen[j][i]) {
                kdp->mems[i][nmem] = j;
                nmem++;
            }
        }
        kdp->nmem[i] = nmem;
    }
    // find SNP members
    for (i = 0; i < nsnp; i++) {
        n = ngen + i;
        nmem = 0;
        for (j = 0; j < nkit; j++) {
            if (kdp->ksnp[j][i] == 1) {
                kdp->mems[n][nmem] = j;
                nmem++;
            }
        }
        kdp->nmem[n] = nmem;
    }
}

int
downstream(struct KITDAT *kdp, int k, int n)
{
    int i, m, nsnp, ngen;

    nsnp = kdp->nsnp; 
    ngen = kdp->ngen; 
    for (i = 0; i < nsnp; i++) {
        m = ngen + i;
        if ((kdp->gsnp[k][i] == 1) && (kdp->parent[m] == n)) break;
    }
    if (i < nsnp) n = downstream(kdp, k, m);
    return (n);
}

void
find_parents(struct KITDAT *kdp)
{
    int ngen, nsnp, nkit, npar;
    int i, k, n, r, p, parent;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    // find GEN parents
    for (i = 0; i < ngen; i++) {
        npar = nkit;
        parent = 0;
        for (k = 1; k < ngen; k++) {
            if (k == i) continue;
            if (memsubset(kdp, i, k)) {
                if (npar > kdp->nmem[k]) {
                    npar = kdp->nmem[k];
                    parent = k;
                }
            }
        }
        kdp->parent[i] = parent;
    }
    // find SNP parents
    for (i = 1; i < nsnp; i++) {
        n = ngen + i;
        if (kdp->nmem[n] > 0) kdp->parent[n] = ngen + kdp->pars[i];
    }
    // find SNP per GEN
    infer_snps(kdp);
    // merge GEN & SNP parents
    r = 0; while (kdp->gsnp[0][r] == 1) r++;
    for (i = r; i < nsnp; i++) {
        n = ngen + i;
        for (k = 1; k < ngen; k++) {
            p = kdp->parent[k];
            if ((kdp->gsnp[k][i] == 1) && (kdp->gsnp[p][i] == 2)) {
                //kdp->parent[k] = n;
                kdp->parent[k] = downstream(kdp, k, n);
                if (p) kdp->parent[n] = p;
            }
        }
    }
    for (k = 1; k < ngen; k++) {
        p = kdp->parent[k];
        if (p == 0) kdp->parent[k] = ngen;
    }
}

void
find_kids(struct KITDAT *kdp)
{
    int ngen, nsnp, nkid, nnod;
    int i, k;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    // find NOD kids
    for (i = 0; i < nnod; i++) {
        nkid = 0;
        for (k = 1; k < nnod; k++) {
            if ((kdp->parent[k] == i) && (kdp->nmem[k] > 0)) {
                kdp->kids[i][nkid] = k;
                nkid++;
            }
        }
        kdp->nkid[i] = nkid;
    }
}

int
get_descendants(struct KITDAT *kdp, int n, char *dscs, int incmem)
{
    char *mems, *kids;
    int i, nmem, nkid, ndsc;

    nkid = kdp->nkid[n];
    kids = kdp->kids[n];
    nmem = kdp->nmem[n];
    mems = kdp->mems[n];
    if (incmem) {
        for (i = 0; i < nmem; i++) {
            dscs[i] = mems[i];
        }
        ndsc = nmem;
    } else {
        ndsc = 0;
    }
    for (i = 0; i < nkid; i++) {
        ndsc += get_descendants(kdp, kids[i], dscs + ndsc, 1);
    }
    return (ndsc);
}

void
find_leftovers(struct KITDAT *kdp)
{
    char *mems, select[MAXKIT], dscs[MAXKIT];
    int nkit, ngen, nsnp, nnod, ndsc, nmem, nlft;
    int i, k, m;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    // find NOD kids
    for (i = 0; i < nnod; i++) {
        nmem = kdp->nmem[i];
        mems = kdp->mems[i];
        for (k = 0; k < nkit; k++) select[k] = 0;
        for (m = 0; m < nmem; m++) select[(int)mems[m]] = 1;
        ndsc = get_descendants(kdp, i, dscs, 0);
        for (k = 0; k < ndsc; k++) select[(int)dscs[k]] = 0;
        nlft = 0;
        for (k = 0; k < nkit; k++) {
            if (select[k]) {
                kdp->lfts[i][nlft] = k;
                nlft++;
            }
        }
        kdp->nlft[i] = nlft;
    }
}

void
process_nodes(struct KITDAT *kdp)
{
    char *s;
    int ngen, nsnp, nchr, nnod;
    int i, n, y, sumlft;

    find_members(kdp);
    find_parents(kdp);
    find_kids(kdp);
    find_leftovers(kdp);
    // report
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    if (debug) {
        // find max name length
        nchr = 0;
        for (i = 0; i < ngen; i++) {
            n = strlen(kdp->gens[i]);
            if (nchr < n) nchr = n;
        }
        for (i = 0; i < nsnp; i++) {
            n = strlen(kdp->snps[i]);
            if (nchr < n) nchr = n;
        }
        // debug print
        sumlft = 0;
        for (i = 0; i < nnod; i++) {
            s = (i < ngen) ? kdp->gens[i] : kdp->snps[i - ngen];
            y = (i < ngen) ? kdp->geny[i] : kdp->snpy[i - ngen];
            printf("%2d %*s.%4d  ", i, nchr, s, y);
            printf("p=%2d  nmkl=%2d ", kdp->parent[i], kdp->nmem[i]);
            printf("%2d %2d\n", kdp->nkid[i], kdp->nlft[i]);
            sumlft += kdp->nlft[i];
        }
    }
    // sanity check
    sumlft = 0; for (i = 0; i < nnod; i++) sumlft += kdp->nlft[i]; 
    if (sumlft != kdp->nkit) {
        printf("WARNING: sumlft=%d != nkit=%d\n", sumlft, kdp->nkit);
    } else {
        printf("--------------------\n");
    }
}

void
str_adjust(struct KITDAT *kdp)
{
    int j;

    for (j = 0; j < kdp->nkit; j++) {
        kdp->kstrs[j][11] -= kdp->kstrs[j][9]; // adjust 
    }
}

struct KITDAT *
get_sapp(char *ifn)
{
    char *s, line[MAXLIN];
    FILE *ifp;
    static struct KITDAT kd;

    kd.nkit = kd.nstr = kd.nsnp = 0;
    init_gendat(&kd);
    ifp = fopen(ifn, "r"); // input file
    while (fgets(line, MAXLIN, ifp)) {
        if (line[0] == '/') break;
    }
    while (line[0] == '/') {
        if (strncmp(line, "/STRDATA", 8) == 0) {
            get_strdata(&kd, ifp, line);
        } else if (strncmp(line, "/MODAL", 5) == 0) {
            get_modal(&kd, ifp, line);
        } else if (strncmp(line, "/SNPDATA", 8) == 0) {
            get_snpdata(&kd, ifp, line);
        } else if (strncmp(line, "/SNPTREE", 8) == 0) {
            get_snptree(&kd, ifp, line);
        } else if (strncmp(line, "/GENDATA", 8) == 0) {
            get_gendata(&kd, ifp, line);
        } else {
            while ((s = fgets(line, MAXLIN, ifp))) {
                if (line[0] == '/') break;
            }
            if (s == NULL) line[0] = EOF;
        }
    }
    fclose(ifp);
    get_snpy(&kd);
    process_nodes(&kd);
    kd.nlab = get_strlab(kd.labs);
    kd.tnmv  = find_misval(&kd);
    str_adjust(&kd);
    sort_order(&kd, so, 1);
    inherit_strs(&kd, kd.ngen); // start inheritance with first SNP
    return (&kd);
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
free_kitdat(struct KITDAT *kdp)
{
    int i, n;

    n = kdp->nkit + kdp->ngen;
    for (i = 0; i < n; i++) {
        free(kdp->kits[i]);
        free(kdp->kstrs[i]);
        free(kdp->ostrs[i]);
    }
    for (i = 0; i < kdp->nlab; i++) {
        free(kdp->labs[i]);
    }
    for (i = 0; i < kdp->nsnp; i++) {
        free(kdp->snps[i]);
    }
    for (i = 0; i < kdp->ngen; i++) {
        free(kdp->gens[i]);
    }
}

char *
d(int n)
{
    static char *s[5] = {"?", "+", "-", "?+", "?-"};

   return (s[n]);
}

void
kit_snpgen_csv(struct KITDAT *kdp)
{
    char ofn[MAXFNL];
    int i, j, nsnp, ngen;
    FILE *ofp;
    static char *name = "kit_snpgen";

    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" output: %s\n", ofn);
    nsnp = kdp->nsnp;
    ngen = kdp->ngen;
    ofp = fopen(ofn, "w");
    fputs("kit,", ofp);
    for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", kdp->snps[i]);
    for (i = 1; i < ngen; i++) fprintf(ofp, "%s,", kdp->gens[i]);
    fputs("\n", ofp);
    for (j = 0; j < kdp->nkit; j++) {
        fprintf(ofp,"%s,", kdp->kits[j]);
        for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", d(kdp->ksnp[j][i]));
        for (i = 1; i < ngen; i++) fprintf(ofp, "%s,", d(kdp->kgen[j][i]));
        fputs("\n", ofp);
    }
    fclose(ofp);
}

void
gen_snp_csv(struct KITDAT *kdp)
{
    char ofn[MAXFNL];
    int i, j, nsnp, ngen;
    FILE *ofp;
    static char *name = "gen_snp";

    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" output: %s\n", ofn);
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    ofp = fopen(ofn, "w");
    fputs("name,", ofp);
    for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", kdp->snps[i]);
    fputs("\n", ofp);
    for (j = 0; j < ngen; j++) {
        if (j == 0) {
            fprintf(ofp,"%s,", kdp->gens[0]);
        } else {
            fprintf(ofp,"%s.%d,", kdp->gens[j], kdp->geny[j]);
        }
        for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", d(kdp->gsnp[j][i]));
        fputs("\n", ofp);
    }
    fclose(ofp);
}

void
str_all_csv(struct KITDAT *kdp)
{
    char ofn[MAXFNL];
    int i, j, n, nlab, val;
    FILE *ofp;
    static char *name = "sv_all";

    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" output: %s\n", ofn);
    nlab = kdp->nlab;
    ofp = fopen(ofn, "w");
    fputs("kit,", ofp);
    for (i = 0; i < nlab; i++) {
        fputs(kdp->labs[i], ofp);
        if (i < (nlab - 1)) {
            fputs(",", ofp);
        } else {
            fputs("\n", ofp);
        }
    }
    n = kdp->nkit + kdp->ngen + kdp->nsnp;
    for (j = 0; j < n; j++) {
        fprintf(ofp,"%s,", kdp->kits[j]);
        for (i = 0; i < nlab; i++) {
            if (i < kdp->nsvs[j]) {
                val = kdp->kstrs[j][i];
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

int
nodout(struct KITDAT *kdp, FILE *ofp, int n, int ndnt, int snpl)
{
    char *sn, name[40];
    int b, i, k, mv, pv, nsig, beg, len;
    static int in = 2;
    static int end = 72;

    if (n == 0) {
        strcpy(name, kdp->gens[n]);
        ndnt = strlen(name);
    } else {
        node_name(kdp, name, n);
        ndnt += in;
    }
    fputs(name, ofp);
    len = strlen(name);
    beg = ndnt + len - in;
    if (lefts) {
        for (i = 0; i < kdp->nlft[n]; i++) {
            k = kdp->lfts[n][i];
            sprintf(name, " %s", kdp->kits[k]);
            if ((len + strlen(name)) > end) {
                fprintf(ofp, "\n%*s", beg, "");
                len = beg;
            }
            fputs(name, ofp);
            len += strlen(name);
        }
    }
    nsig = (lefts > 1) ? kdp->sig[n].n : 0;
    if (nsig) {
        for (b = 0; b < nsig; b++) {
            sn = kdp->labs[kdp->sig[n].s[b]];
            pv = kdp->sig[n].p[b];
            mv = kdp->sig[n].m[b];
            sprintf(name, " %s=%d->%d", sn, pv, mv);
            if ((len + strlen(name)) > end) {
                fprintf(ofp, "\n%*s", beg, "");
                len = beg;
            }
            fputs(name, ofp);
            len += strlen(name);
        }
    }
    for (i = 0; i < kdp->nkid[n]; i++) {
        fprintf(ofp, "\n%*s", ndnt, ". ");
        snpl = nodout(kdp, ofp, kdp->kids[n][i], ndnt, snpl);
    }
    return (snpl);
}

void
report_snpgen(struct KITDAT *kdp)
{
    char ofn[MAXFNL];
    FILE *ofp;
    static char *name = "tree";

    if (group) {
        sprintf(ofn, "%s%d.txt", name, group);
    } else {
        sprintf(ofn, "%s.txt", name);
    }
    printf(" output: %s\n", ofn);
    ofp = fopen(ofn, "w");
    nodout(kdp, ofp, 0, 0, 0);
    fputs("\n", ofp);
    fclose(ofp);
}

void
report_sapp(char *ifn, struct KITDAT *kdp)
{
    printf("grpsum: %s nkit=%d nstr=%d ", ifn, kdp->nkit, kdp->nstr);
    printf("nsnp=%d ngen=%d tnmv=%d\n", kdp->nsnp, kdp->ngen, kdp->tnmv);
}

int
main(int argc, char **argv)
{
    char ifn[MAXFNL];
    int allcsv;
    struct KITDAT *kdp;

    if (argc < 2) {
        printf("usage: grpsum [options] filename [group]\n");
        printf("options:\n");
        printf("  -a    make CSV continaing all GENs & STRs\n");
        printf("  -d    enable debug print\n");
        printf("  -l    list leftovers\n");
        printf("  -m    list leftovers + STR mutations\n");
	exit(1);
    }
    allcsv = 0;
    if (argc > 2) {
        while (argv[1][0] == '-') {
            if (argv[1][1] == 'a') {
                allcsv++;
            } else if (argv[1][1] == 'd') {
                debug++;
            } else if (argv[1][1] == 'l') {
                lefts++;
            } else if (argv[1][1] == 'm') {
                lefts = 2;
            }
            argv++;
            argc--;
        }
    }
    strncpy(ifn, argv[1], MAXFNL);
    if (argc > 2) {
        group = atoi(argv[2]);
    } else if (strncmp(ifn, "Neely_s", 7) == 0) {
        group = atoi(ifn + 7);
    } else {
        group = 0;
    }
    // input KIT data
    kdp = get_sapp(ifn);
    report_sapp(ifn, kdp);
    // report SNPs & GENs
    report_snpgen(kdp);
    // output CSV
    if (allcsv) str_all_csv(kdp);
    kit_snpgen_csv(kdp);
    gen_snp_csv(kdp);
    // clean up
    free_kitdat(kdp);
} 

