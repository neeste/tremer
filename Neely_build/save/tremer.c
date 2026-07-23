// tremer - merge GEN, SNP, & STR trees

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "strlab.h"
#include "tremer.h"
#include "utils.h" 

static int debug = 0;
static int count = 0;
static int lefts = 2; // report leftovers + STR signature
static int modal = 0;
static int muthr = 2; // STR mutation threshold

static int
get_kit(FILE *ifp, char *hdr, char *kit, char *kstr, char *line)
{
    char *s;
    int nsv;

    if (*line == '*') {
        strncpy(hdr, line, MAXHDR);
        fgets(line, MAXLIN, ifp);
    } else {
        *hdr=0;
    }
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

static int
get_strlab(char *strlab[])
{
    char line[MAXLIN], lab[80];
    int i, j, n, nlab;
    FILE *lfp;
    static char lfn[] = "str_label.csv";

    lfp = fopen(lfn, "r");
    if (!lfp) {
        for (i = 0; i < n_lab; i++) {
            strlab[i] = strdup(str_lab[i]);
        }
	return (n_lab);
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
    if (0) {
        FILE *tfp;
        static char tfn[] = "str_label.txt";
        tfp = fopen(tfn, "w");
        fprintf(tfp, "static int nlab = %d;\n", nlab);
        fprintf(tfp, "static char *str_lab[%d] = {\n    ", nlab);
        for (i = 0; i < nlab; i++) {
            fprintf(tfp, "\"%s\", ", strlab[i]);
            if (((i+1)%9)==0) fprintf(tfp, "\n    ");
        }
        fprintf(tfp, "    \n};\n");
        fclose(tfp);
    }
    return(nlab);
}

static int
find_misval(struct KITDAT *kdp)
{
    short *nmvs;
    int i, j, nmv, tnmv;

    nmvs = kdp->nmvs;
    tnmv = 0;
    for (i = 0; i < kdp->nkit; i++) {
        nmv = 0;
        for (j = 0; j < kdp->nsvs[i]; j++) {
            if (get_strk(i,j) == 0) nmv++;
        }
        nmvs[i] = nmv;
        tnmv += nmv;
    }
    return (tnmv);
}

static int
modal_val(char *v, int n)
{
    char uv[MAXSTR];
    int i, j, nuv, no, mxno, mv;

    nuv = unique(v, uv, n, 1);
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

static int
node_name(struct KITDAT *kdp, char *name, int node)
{
    int y, nkit;

    nkit = kdp->nkit;
    sprintf(name, "%s", kdp->kits[nkit + node]);
    y = kdp->kity[nkit + node];
    if (y) sprintf(name + strlen(name), ".%d", y); // append year
    return(strlen(name));
}

static int
unique_strmut(struct KITDAT *kdp, short *usi, char *mems, int nmem)
{
    char sv[MAXSV], uv[MAXSV];
    int i, j, k, n, nuv, nus, nstr;

    nstr = kdp->nstr;
    nus = 0;
    for (i = 0; i < nstr; i++) {
        for (j = 0; j < nmem; j++) {
            k = mems[j];
            n = kdp->nsvs[k];
            sv[k] = (i < n) ? kdp->kstrs[k][i] : 0;
        }
        nuv = unique(sv, uv, nmem, 1);
        if (nuv) usi[nus++] = (short)i;
    }
    return (nus); // number of unique STR mutations in set of kits
}

static int
sv_count(struct KITDAT *kdp, char *mems, int nmem, int s, int m)
{
    int i, j, v, c;

    c = 0;
    for (i = 0; i < nmem; i++) {
        j = mems[i];
        v = get_stro(j,s);
        if (v == m) c++;
    }
    return (c);
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
str_adjust(struct KITDAT *kdp)
{
    int j, n;

    n = kdp->nkit + kdp->ngen + kdp->nsnp + kdp->nsbg;
    for (j = 0; j < n; j++) {
        if (kdp->kstrs[j][11] > 22)
            kdp->kstrs[j][11] -= kdp->kstrs[j][9]; // adjust 389ii
        if (kdp->ostrs[j][11] > 22)
            kdp->ostrs[j][11] -= kdp->ostrs[j][9]; // adjust 389ii
    }
}

static void
init_modal(struct KITDAT *kdp)
{
    char v[MAXKIT], sv[MAXSTR];
    int i, j, nsv, nkit, nstr;

    if (modal) return; // check whether already initialied
    nkit = kdp->nkit;
    nstr = kdp->nstr;
    for (i = 0; i < nstr; i++) {
        for (j = 0; j < nkit; j++) {
            nsv = kdp->nsvs[j];
            v[j] = (i < nsv) ? get_strk(j,i) : 0;
        }
        sv[i] = modal_val(v, nkit);
    }
    add_kit(kdp, "modal", 0, sv, nstr);
    kdp->nkit--; // ancestral STR not a kit
    modal++;
}

static int
modal_set(struct KITDAT *kdp, char *sv, char *mems, int nmem)
{
    char v[MAXKIT];
    int i, j, k, nsv, nstr;

    nstr = kdp->nstr;
    for (i = 0; i < nstr; i++) {
        for (j = 0; j < nmem; j++) {
            k = mems[j];
            nsv = kdp->nsvs[k];
            v[j] = (i < nsv) ? get_stro(k,i) : 0;
        }
        sv[i] = modal_val(v, nmem);
    }
    return (nstr);
}

static void
patch_strs(struct KITDAT *kdp, int node)
{
    char mv, sv, v[MAXKIT];
    int i, m, k, nsv, kidk, mm, nmem, nkid, nstr;

    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        patch_strs(kdp, kidk);
    }
    nmem = kdp->nmem[node];
    if (nmem < 4) return;
    nstr = kdp->nstr;
    for (i = 0; i < nstr; i++) {
        for (m = 0; m < nmem; m++) {
            mm = kdp->mems[node][m];
            nsv = kdp->nsvs[mm];
            v[m] = (i < nsv) ? get_strk(mm,i) : 0;
        }
        mv = modal_val(v, nmem);       // node-modal value
        for (m = 0; m < nmem; m++) {
            mm = kdp->mems[node][m];
            sv = get_strk(mm,i);
            if (!sv) set_strk(mm,i,mv); // set missing value to modal value
        }
    }
}

static void
name_branch(struct KITDAT *kdp, int node)
{
    char *sbgs;
    int isbg, nkit, ngen, nsnp;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    // check if node is unnamed
    isbg = ngen + nsnp;
    if (node >= isbg) {
        sbgs = kdp->kits[nkit + node];
        if (sbgs[0] == '.') sprintf(sbgs, "STR%02d", (node - isbg) + 1);
        strcpy(kdp->sbgs[node - isbg], sbgs);
    }
}

static int
add_branch(struct KITDAT *kdp, int node)
{
    char *sbgs, *kstr;
    int bran, nsv, nstr, nsbg, nkit, nkid, ngen, nsnp;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkid = kdp->nkid[node];
    // create new branch
    nsbg = kdp->nsbg;
    kdp->nsbg = nsbg + 1;
    sbgs = (char *)malloc(40);
    strcpy(sbgs, "."); // default branch name
    kdp->sbgs[nsbg] = sbgs;
    // add new kit
    kstr = kdp->kstrs[node + nkit];
    nsv = kdp->nsvs[node + nkit];
    kdp->nkit = nkit + ngen + nsnp + nsbg;
    add_kit(kdp, sbgs, 0, kstr, nsv);
    kdp->nkit = nkit;
    // add new node
    bran = ngen + nsnp + nsbg;
    kdp->parent[bran] = node;
    kdp->kids[node][nkid++] = bran;
    kdp->nkid[node] = nkid;
    // copy STRs from node to bran and clear sig
    nstr = kdp->nstr;
    memcpy(kdp->kstrs[bran + nkit], kdp->kstrs[node + nkit], nstr);
    kdp->sig[bran].n = 0;
    return (bran);
}

static int
eat_leftovers(struct KITDAT *kdp, int node, int i, int mv)
{
    char u[MAXKIT], mvl[MAXKIT], mvk[MAXGEN];
    char nodnam[40];
    char *str, *lfts, *mems; 
    int nsv, nlft, nkid, nmem;
    int cmv, k, m, tgt, mc, lc, lft, kid, mem, kc;

    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    lc = 0; // kit count
    for (m = 0; m < nlft; m++) {
        lft = lfts[m];
        nsv = kdp->nsvs[lft];
        str = kdp->kstrs[lft];
        u[m] = (i < nsv) ? str[i] : 0;
        if (u[m] == mv) {
            mvl[lc] = m;
            lc++;
        }
    }
    nkid = kdp->nkid[node];
    // find number of kids containing modal value
    kc = 0; // kid count
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
            mvk[kc] = kid;
            kc++;
        }
    }
    node_name(kdp, nodnam, node);
    if (nkid == 0) return (lc);
    cmv = tgt = 0;
    if (kc == nkid) {
        cmv = (nkid == 1) ? lc : 1;
    } else if (kc == 1) {
        cmv = lc;
        tgt = 0;
    }
    return (cmv);
}

static void
str_sort(struct KITDAT *kdp, short *so, int rank)
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
	    v[k] = (i < kdp->nsvs[k]) ? get_strk(k,i) : 0;
        }
        nuv = unique(v, uv, nkit, 1);
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

static int
sv_histo(struct KITDAT *kdp, char *mems, int nmem, short si,
    char *sv, short *nv)
{
    char uv[MAXSV], svtmp, *s;
    int i, j, n, nuv;
    short nvtmp;

    // make list of all STR values across members
    for (i = 0; i < nmem; i++) {
        j = mems[i];
        n = kdp->nsvs[j];
        sv[i] = (i < n) ? get_stro(j,si) : 0;
    }
    // make list of unique values across list of all STR values
    nuv = unique(sv, uv, nmem, 1);
    for (j = 0; j < nuv; j++) {
        nv[j] = 0;
        for (i = 0; i < nmem; i++) {
            if (sv[i] && (sv[i] == uv[j])) nv[j]++;
        }
    }
    // sort list of unique value counts
    for (j = 0; j < nuv; j++) sv[j] = uv[j]; // copy unique values
    for (j = 0; j < nuv; j++) {              // largest-count-first order
        for (i = j + 1; i < nuv; i++) {
            if (nv[i] > nv[j]) {     // exchange ??
                nvtmp = nv[i];
                nv[i] = nv[j];
                nv[j] = nvtmp;
                svtmp = sv[i];
                sv[i] = sv[j];
                sv[j] = svtmp;
            }
        }
    }
    if (debug && 0) {
        if ((nuv > 1) && (nv[1] > 1)) {
            s = str_lab[si];
            printf("sv_histo: %6s nmem=%2d nuv=%d  ", s, nmem, nuv);
            printf("sv=%2d %2d  nv=%2d %2d\n", sv[0], sv[1], nv[0], nv[1]);
        }
    }
    return (nuv);  // number of histo values
}

static int
comodl(struct KITDAT *kdp, char *mems, int nmem, short si, char *cmsv)
{
    char sv[MAXSTR];
    int nuv, ncv;
    short nv[MAXSTR];

    nuv = sv_histo(kdp, mems, nmem, si, sv, nv);
    if (nuv > 1) {
        if (cmsv) *cmsv = sv[1]; // co-modal STR value
        ncv = nv[1];             // number of co-modal values
    } else {
        ncv = 0;                 // no co-modal values
    }
    return (ncv);
}

static void
inherit(struct KITDAT *kdp, int node)
{
    char *s, *str, *mems, kit[80], par[80], v[MAXKIT];
    int nsv, nkid, ngen, nkit, nmem, nstr, prnt;
    int i, ii, k, m, mem, cmv, cc, mc, mv, pc, pv;

    if (node) {
        nstr = kdp->nstr;
        ngen = kdp->ngen;
        nkit = kdp->nkit;
        k = node + nkit;
        prnt = kdp->parent[node];
        s = (prnt < ngen) ? kdp->gens[prnt] : kdp->snps[prnt - ngen];
        strcpy(par, s);
        s = (node < ngen) ? kdp->gens[node] : kdp->snps[node - ngen];
        strcpy(kit, s);
        nsv = kdp->nsvs[nkit];
        kdp->nsvs[k] = nsv;
        strcpy(kdp->kits[k], kit);
        nmem = kdp->nmem[node];
        mems = kdp->mems[node];
        for (i = 0; i < nstr; i++) {
            ii = kdp->sso[i];  // sorted order
            for (m = 0; m < nmem; m++) {
                mem = mems[m];
                nsv = kdp->nsvs[mem];
                str = kdp->kstrs[mem];
                v[m] = (ii < nsv) ? str[ii] : 0;
            }
            pv = get_strn(prnt, ii);
            mv = modal_val(v, nmem);
            pc = count_val(v, nmem, pv);
            mc = count_val(v, nmem, mv);
            cmv = 0;
            if ((mv && pv) && (mv != pv) && (mc > pc)) {
                cmv = eat_leftovers(kdp, node, ii, mv);
                if (cmv < 0) {
                    inherit(kdp, node);
                    return;
                }
            }
            if (cmv) {
                cc = comodl(kdp, mems, nmem, ii, 0);
                if (mc > 1) { // <-- DEBUG
                    set_strn(prnt, ii, pv);
                    set_strn(node, ii, mv);
                }
                if (debug && 0) {
                    char nnm[40];
                    node_name(kdp, nnm, node);
                    printf("inherit: %18s =%2d->%2d mc=%2d cc=%2d nmem=%2d\n",
                        nnm, pv, mv, mc, cc, nmem);
                }
            } else {
                set_strk(k, ii, pv);
            }
        }
    }
    // inherit kids
    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        inherit(kdp, kdp->kids[node][k]);
    }
}

static void
reset_sig(struct KITDAT *kdp, int node)
{
    char *mems, *lfts;
    int c, d, i, n, s, p, m, prnt, nmem, nkid, nlft, nus;
    short si[MAXSTR];
    static char *fnm = "reset_sig";

    if (node == 0) return;
    nmem = kdp->nmem[node];
    mems = kdp->mems[node];
    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    nkid = kdp->nkid[node];
    prnt = kdp->parent[node];
    nus = unique_strmut(kdp, si, mems, nmem);
    n = 0;
    for (i = 0; i < nus; i++) {
        s = si[i];
        p = get_strn(prnt,s);
        m = get_strn(node,s);
        c = sv_count(kdp, mems, nmem, s, m);
        d = sv_count(kdp, lfts, nlft, s, m);
        if (m && p && (m != p) && ((c > 1) || (d > 1))) {
            kdp->sig[node].s[n] = s;
            kdp->sig[node].p[n] = p;
            kdp->sig[node].m[n] = m;
            kdp->sig[node].c[n] = c;
            kdp->sig[node].d[n] = d;
            n++;
        }
    }
    kdp->sig[node].n = n;
    if (debug && n && 0) {
        char *l, nnm[40];
        node_name(kdp, nnm, node);
        n = kdp->sig[node].n;
        printf("%s: %18s prnt=%2d ", fnm, nnm, prnt);
        printf("nmem=%2d nus=%3d nsig=%d ", nmem, nus, n);
        if (n > 0) {
            s = kdp->sig[node].s[0];
            p = kdp->sig[node].p[0];
            m = kdp->sig[node].m[0];
            c = kdp->sig[node].c[0];
            d = kdp->sig[node].d[0];
            l = (nkid && (c == d)) ? "*" : "";
            printf("%6s=%2d->%2d|%d%s ", kdp->labs[s], p, m, c, l);
        }
        printf("\n");
    }
}

static void
update_sig(struct KITDAT *kdp, int node)
{
    int k, nkid;

    reset_sig(kdp, node);
    nkid = kdp->nkid[node];
    if (debug && 0) {
        char nnm[40];
        node_name(kdp, nnm, node);
        printf("update_sig: %18s nkid=%2d\n", nnm, nkid);
    }
    for (k = 0; k < nkid; k++) {
        update_sig(kdp, kdp->kids[node][k]);
    }
}

static void
get_strdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, hdr[MAXHDR], kstr[MAXSTR], kit[80];
    int nsv;

    while ((s = fgets(line, MAXLIN, ifp))) {
        if (s == NULL) { line[0] = EOF; return; }
        if (line[0] == '/' || line[0] < ' ') break;
        nsv = get_kit(ifp, hdr, kit, kstr, line);
        if (*kit) {
            add_kit(kdp, kit, 0, kstr, nsv);
            //printf("%8s %2d %2d %2d\n", kit, kstr[1], kstr[2], kstr[3]);
        }
    }
    
}

static void
get_modal(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kstr[MAXSTR];
    int nsv;

    while ((s = fgets(line, MAXLIN, ifp))) {
        if (s == NULL) { line[0] = EOF; return; }
        if (line[0] == '/' || line[0] < ' ') break;
        nsv = getval(line, kstr);
        if (nsv) {
            add_kit(kdp, "ancestral", 0, kstr, nsv);
            kdp->nkit--; // ancestral STR not a kit
            modal++;
            //printf("get_modal: %d %d %d\n", kstr[0], kstr[1], kstr[2]);
        }
    }
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
        if (s == NULL) { line[0] = EOF; return; }
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
    kdp->nsnp = nsnp;
}

static void
get_snptree(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, snpnam[80], snppar[80], snpyear[80];
    int i, j, n, p, y, isnp, nsnp;

    nsnp = kdp->nsnp;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (s == NULL) { line[0] = EOF; return; }
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        while (i < n) {
            while (isspace(line[i])) i++;
            if (line[i] == 0 || line[0] == '*') break;
            j = 0;
            while (isalnum(line[i])) snpnam[j++] = line[i++];
            if (j) {
                snpnam[j] = 0;
                for (isnp = 0; isnp < nsnp; isnp++) {
                    if (strcmp(snpnam, kdp->snps[isnp]) == 0) break;
                }
                while (isspace(line[i])) i++;
                j = 0;
                while (isalnum(line[i])) snpyear[j++] = line[i++];
                snpyear[j] = 0;
                y = atoi(snpyear);
                while (isspace(line[i])) i++;
                j = 0;
                while (isalnum(line[i])) snppar[j++] = line[i++];
                snppar[j] = 0;
                for (p = 0; p < nsnp; p++) {
                    if (strcmp(snppar, kdp->snps[p]) == 0) break;
                }
                if (p == nsnp) {
                    p = 0;
                    snppar[0] = 0;
                }
                while (line[i]) i++; // skip to EOL
//printf("%9s i=%2d y=%4d p=%2d %s\n", snpnam, isnp, y, p, snppar);
                if (y) kdp->snpy[isnp] = y;
                if (p) kdp->pars[isnp] = p;
            }
        }
    }
}

static int
isalnumeq(char c)
{
    return (isalnum(c) || (c == '='));
}

static void
get_gendata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kit[80], gennam[80];
    int i, j, n, c, d, ikit, igen, nkit, ngen, nmem;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (s == NULL) { line[0] = EOF; return; }
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        j = 0;
        while (isalnumeq(line[i])) gennam[j++] = line[i++];
        if (j == 0) continue;
        gennam[j] = 0;
        d = 1;
        c = gennam[j - 1];
        if (strchr("+-?", c)) {
            if (c == '?') d = 0;
            if (c == '+') d = 1;
            if (c == '-') d = 2;
            gennam[j - 1] = 0;
        }
        kdp->gens[ngen] = strdup(gennam);
        if (line[i] == '.') {
            kdp->geny[ngen] = atoi(line + i + 1);
        }
        igen = ngen;
        ngen++;
        nmem = 0;
        while ((i < n) && (line[i] != '(')) i++;
        while ((i < n) && (line[i] != ')')) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) kit[j++] = line[i++];
            if (j) {
                kit[j] = 0;
                ikit = find_kit(kdp, kit);
                if (ikit < nkit) {
                    kdp->kgen[ikit][igen] = d;
                    nmem++;
                } else {
                    printf("WARNING: kit %s has no STR\n", kit);
                }
                if (debug && 0) {
                    printf("%2d %17s %7s : %2d %2d %2d\n",
                        ngen, gennam, kit, igen, ikit, nmem);
                }
            }
            while ((i < n) && !isalnum(line[i]) && (line[i] != ')')) i++;
        }
    }
    kdp->ngen = ngen;
}

static void
get_grpdata(struct KITDAT *kdp, FILE *ifp, char *line)
{
    char *s, kit[80], grpnam[80];
    int i, j, n, ikit, igrp, nkit, ngrp, nmem;

    nkit = kdp->nkit;
    ngrp = kdp->ngrp;
    while ((s = fgets(line, MAXLIN, ifp))) {
        if (s == NULL) { line[0] = EOF; return; }
        if (line[0] < ' ') break;
        if (line[0] == '/') break;
        i = 0;
        n = strlen(line);
        j = 0;
        while (isalnum(line[i])) grpnam[j++] = line[i++];
        if (j == 0) continue;
        grpnam[j] = 0;
        kdp->grps[ngrp] = strdup(grpnam);
        igrp = ngrp;
        ngrp++;
        nmem = 0;
        while ((i < n) && (line[i] != '(')) i++;
        while ((i < n) && (line[i] != ')')) {
            while (isspace(line[i])) i++;
            j = 0;
            while (isalnum(line[i])) kit[j++] = line[i++];
            if (j) {
                kit[j] = 0;
                ikit = find_kit(kdp, kit);
                if (ikit < nkit) {
                    kdp->kgrp[ikit] = igrp;
                    nmem++;
                } else {
                    printf("WARNING: kit %s has no STR\n", kit);
                }
            }
            while ((i < n) && !isalnum(line[i]) && (line[i] != ')')) i++;
        }
        if (debug) {
            printf("get_grpdata(%d): %s nmem=%d\n", igrp, grpnam, nmem);
        }
    }
    kdp->ngrp = ngrp;
    if (debug) {
        printf("get_grpdata: ngrp=%d\n", ngrp);
    }
}

static int
intersect(char *mi, char *mn, char *mp, int nn, int np, int in)
{
    char mt[800];
    int i, j, ni;

    ni = 0;
    for (i = 0; i < nn; i++) {
        mt[ni] = mn[i];
        for (j = 0; j < np; j++) {
            if (mn[i] == mp[j]) break;
        }
        if (in && (j < np)) ni++;   // find inclusions
        if (!in && (j == np)) ni++; // find exclusions
    }
    memcpy(mi, mt, ni);
    return (ni);
}

static void
excl_gen(struct KITDAT *kdp)
{
    char *m1, *m2, gkit[MAXGEN][MAXKIT], mi[MAXKIT], me[MAXKIT];
    int i, j, k, n, ni, ne, n1, n2, nkit, ngen, nkpg[MAXGEN];

    // group kits in each GEN node
    nkit = kdp->nkit;
    ngen = kdp->ngen;
    for (i = 0; i < ngen; i++) {
        j = 0;
        for (k = 0; k < nkit; k++) {
            if (kdp->kgen[k][i] == 1) gkit[i][j++] = k;
        }
        nkpg[i] = j;
    }
    // find GEN descendants
    for (i = 1; i < ngen; i++) {
        for (j = 1; j < ngen; j++) {
            if (i == j) continue;
            m1 = gkit[i];
            m2 = gkit[j];
            n1 = nkpg[i];
            n2 = nkpg[j];
            if (!n1 || !n2) continue;
            ni = intersect(mi, m2, m1, n2, n1, 1); // inclusions
            ne = intersect(me, m2, m1, n2, n1, 0); // exclusions
            if (ni && ne) {
                for (n = 0; n < ne; n++) {
                    k = me[n];
                    kdp->kgen[k][i] = 2;
                }
            }
        }
    }
}

void
fetch_snptree(struct KITDAT *kdp)
{
    char snpnam[80], snppar[80], line[MAXLIN];
    int i, j, k, y;
    FILE *ifp;
    static char *ifn = "snptree.txt";
 
    if (debug) printf("->fetch_snptree\n");
    ifp = fopen(ifn, "r");
    if (ifp == NULL) {
        //printf("WARNING: can't open %s\n", ifn);
        return;
    }
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

static void
init_gendat(struct KITDAT *kdp, int group)
{
    char gennam[80];
    int i;

    if (group) {
        sprintf(gennam, "Group%d", group);
    } else {
        strcpy(gennam, "Group");
    }
    kdp->gens[0] = strdup(gennam);
    for (i = 0; i < MAXKIT; i++) kdp->kgen[i][0] = 1;
    kdp->geny[0] = 0;
    kdp->ngen = 1;
    kdp->ngrp = 1;
}

static int
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

static int
pair_snps(int s1, int s2)
{
    char d;
    static char c[5][5] = {{0, 3, 4, 3, 4},
                           {3, 1, 2, 1, 2},
                           {4, 2, 2, 4, 4},
                           {3, 1, 4, 4, 4},
                           {4, 2, 4, 4, 4}};
    d = c[s1][s2];
    return (d);
}

static int
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

static void
kit2gen(struct KITDAT *kdp)
{
    int g, s, nsnp, ngen;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    for (g = 0; g < (ngen + nsnp); g++) {
        for (s = 0; s < nsnp; s++) {
            kdp->gsnp[g][s] = combine_snps(kdp, g, s);
        }
    }
}

static void
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

static void
infer_snps(struct KITDAT *kdp)
{
    kit2gen(kdp);
    gen2gen(kdp);
}

static int
snp_present(struct KITDAT *kdp, char *mems, int nmem, int isnp)
{
    int j, kit, npr;

    npr = 0;
    for (j = 0; j < nmem; j++) {
        kit = mems[j];
        if (kdp->ksnp[kit][isnp] == 1) npr++;
    }
    return (npr);
}

static int
lookup_node(struct KITDAT *kdp, char *sfn)
{
    char name[40];
    int i, ngen, nsnp, nsbg, nnod;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    for (i = 0; i < nnod; i++) {
        node_name(kdp, name, i);
        if (strcmp(sfn, name) == 0) return (i);
    }
    for (i = 0; i < ngen; i++) {
        if (strcmp(sfn, kdp->gens[i]) == 0) return (i);
    }
    for (i = 0; i < nsnp; i++) {
        if (strcmp(sfn, kdp->snps[i]) == 0) return (i + ngen);
    }
    for (i = 0; i < nsbg; i++) {
        if (strcmp(sfn, kdp->sbgs[i]) == 0) return (i + ngen + nsnp);
    }
    return (0);
}

static int
lookup_kits(struct KITDAT *kdp, char *sfn)
{
    char *lfts;
    int i, j, k, ngen, nsnp, nsbg, nnod, nlft;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    for (i = 0; i < nnod; i++) {
        nlft = kdp->nlft[i];
        lfts = kdp->lfts[i];
        for (j = 0; j < nlft; j++) {
            k = lfts[j];
            if (strcmp(sfn, kdp->kits[k]) == 0) return (i);
        }
    }
    return (0);
}

static void
find_members(struct KITDAT *kdp)
{
    int ngen, nsnp, nkit, nmem;
    int i, j, m, n;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    // find GEN members
    for (i = 0; i < ngen; i++) {
        nmem = 0;
        for (j = 0; j < nkit; j++) {
            if ((i == 0) || (kdp->kgen[j][i] == 1)) {
                kdp->mems[i][nmem] = j;
                nmem++;
            }
        }
        kdp->nmem[i] = nmem;
    }
    // find certain SNP members
    for (i = 0; i < nsnp; i++) {
        n = ngen + i;
        m = 0;
        for (j = 0; j < nkit; j++) {
            if (kdp->ksnp[j][i] == 1) { // certain members
                kdp->mems[n][m] = j;
                m++;
            }
        }
        kdp->nmem[n] = m;
    }
}

static void
snp_logic(struct KITDAT *kdp)
{
    int nsnp, nkit;
    int i, j, k, p;

    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    for (k = 0; k < nkit; k++) {
        // include ancestors
        for (i = 0; i < nsnp; i++) {
            if (kdp->ksnp[k][i] != 1) continue;
            p = kdp->pars[i];
            if (p == 0) continue;
            kdp->ksnp[k][p] = 1;
        }
        // exclude nephews
        for (i = 0; i < nsnp; i++) {
            if (kdp->ksnp[k][i] != 0) continue;
            p = kdp->pars[i];
            if (p == 0) continue;
            if (debug) printf("nephew: p=%d ksnp=%d\n", p, kdp->ksnp[k][p]);
            if (kdp->ksnp[k][p] != 1) continue;
            kdp->ksnp[k][i] = 2;
        }
        // exclude decendants
        for (i = 0; i < nsnp; i++) {
            if (kdp->ksnp[k][i] != 2) continue;
            for (j = 0; j < nsnp; j++) {
                if (j == i) continue;
                p = kdp->pars[j];
                if (p != i) continue;
                kdp->ksnp[k][j] = 2;
            }
        }
    }
}

static void
snp_fill(struct KITDAT *kdp)
{
    int ngen, nsnp, nkit;
    int i, j, m, n;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    // find uncertain SNP members
    for (i = 0; i < nsnp; i++) {
        n = ngen + i;
        m = kdp->nmem[n];
        if (m == 0) continue;
        for (j = 0; j < nkit; j++) {
            if (kdp->ksnp[j][i] == 0) { // uncertain members
                kdp->mems[n][m] = j;
                m++;
            }
        }
        kdp->nmem[n] = m;
    }
}

static void
debug60(struct KITDAT *kdp, char *fnm)
{
    int m, n, p, ngen, nsnp, nnm1;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnm1 = ngen + nsnp - 1;
    if (nnm1 >= 60) {
        printf("%-13s:  ",fnm);
        n = 20;
        m = kdp->nmem[n];
        p = kdp->parent[n];
        printf("nmem[%2d]=%2d prnt=%2d  ",n,m,p);
        n = nnm1;
        m = kdp->nmem[n];
        p = kdp->parent[n];
        printf("nmem[%2d]=%2d prnt=%2d\n",n,m,p);
    }
}

static int
ismember(struct KITDAT *kdp, int kit, int node)
{
    char *mems;
    int i, nmem;

    nmem = kdp->nmem[node];
    mems = kdp->mems[node];
    for (i = 0; i < nmem; i++) if (mems[i] == kit) break;
    return (i < nmem);
}

static int
snp_in_gen(struct KITDAT *kdp, int n, int k)
{
    int i, ikit, isnp, nmem, ngen;

    ngen = kdp->ngen;
    nmem = kdp->nmem[n];
    isnp = n - ngen;
    // find kits with SNP 
    for (i = 0; i < nmem; i++) {
        ikit = kdp->mems[n][i];
        if (kdp->ksnp[ikit][isnp] == 1) {
            if (!ismember(kdp, ikit, k)) return (0);
        }
        if (kdp->ksnp[ikit][isnp] == 2) {
            if (ismember(kdp, ikit, k)) return (0);
        }
    }
    return (1);
}

static int
snp_parent(struct KITDAT *kdp, int n, int k)
{
    int ma, mn, mk, mp, pn, pk, mran;

    mran = kdp->mran;
    pn = kdp->parent[n];
    pk = kdp->parent[k];
    mp = kdp->nmem[pn];
    mn = kdp->nmem[n];
    mk = kdp->nmem[k];
    ma = kdp->nmem[mran];
    if ((ma <= 60) || (ma != mp) || (mp <= mn)) return (0); // <-- DEBUG
    if ((mn <= mk) || (mp <= mk)) return (0);
    if (!snp_in_gen(kdp, n, k)) return (0);
    if (debug) {
        printf("   snp_parent:  ");
        printf("   k=%2d       pk=%2d  ", k, pk);
        printf("   n=%2d       pn=%2d  ", n, pn);
        printf("mp=%2d\n", mp);
    }
    return (1);
}

static void
find_parents(struct KITDAT *kdp)
{
    char *memsk, *memsn, mi[800];
    int ngen, nsnp, nkit, npar, nmemk, nmemn;
    int i, j, k, n, ni, o, p, r, z, adopt, prnt;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    // initialize ancestral SNP parents
    p = 0;
    z = kdp->nmem[0];
    for (i = 0; i < nsnp; i++) {
        k = ngen + i;
        if (kdp->nmem[k] == z) {
            kdp->parent[k] = p;
            p = k;
        }
    }
    // initialize ancestral GEN parents
    for (k = 1; k < ngen; k++) {
        if (kdp->nmem[k] == z) {
            kdp->parent[k] = p;
            p = k;
        }
    }
    kdp->mran = p; // most recent ancestral node
    // find GEN parents
    for (i = 1; i < ngen; i++) {
        if (kdp->nmem[i] == z) continue;
        npar = nkit;
        prnt = 0;
        for (k = 1; k < ngen; k++) {
            if (k == i) continue;
            adopt = 0;
            if (memsubset(kdp, i, k)) {
                if (memsubset(kdp, k, i)) { // flag preSNPs ?
                    if ((i < ngen) && (k > ngen)) adopt++;
                } else if (npar > kdp->nmem[k]) {
                    adopt++;
                }
            }
            if (adopt) {
                npar = kdp->nmem[k];
                prnt = k;
            }
        }
        kdp->parent[i] = prnt;
        while (prnt) {
             if (prnt == i) break;
             prnt = kdp->parent[prnt];
        }
        if (i && prnt) {
            char n1[40], n2[40];
            node_name(kdp, n1, prnt);
            node_name(kdp, n2, i);
            printf("ERROR: parent loop, parent(%s)=%s\n", n1, n2);
            exit(1);
        }
    }
    // find SNP parents
    for (i = 0; i < nsnp; i++) {
        n = ngen + i;
        nmemn = kdp->nmem[n];
        memsn = kdp->mems[n];
        if (nmemn == z) continue;
        if (nmemn > 0) kdp->parent[n] = ngen + kdp->pars[i];
        for (k = 1; k < ngen; k++) {
            nmemk = kdp->nmem[k];
            memsk = kdp->mems[k];
            ni = intersect(mi, memsn, memsk, nmemn, nmemk, 1);
            if (ni < nmemk) continue;
            r = snp_present(kdp, memsk, nmemk, i);
            if (r) {
                p = kdp->parent[k];
                if (!p) kdp->parent[k] = n;
            }
        }
    }
    // find SNP per GEN
    infer_snps(kdp);
    // merge GEN & SNP parents
    if (debug) debug60(kdp,"find_parents1");
    for (i = 0; i < nsnp; i++) {
        n = ngen + i;
        if (kdp->nmem[n] == z) continue;
        p = kdp->parent[n];
        if (kdp->nmem[p] == z) kdp->parent[n] = kdp->mran;
        for (k = 1; k < ngen; k++) {
            if (kdp->nmem[k] == z) continue;
            p = kdp->parent[k];
            if (kdp->nmem[p] == z) p = kdp->parent[k] = kdp->mran;
            if (snp_parent(kdp, n, k)) kdp->parent[n] = k; //<-- DEBUG
            if (debug && 0) {
                if ((n == (ngen + nsnp - 1)) && (p == 20)) {
                    int gsk, gsp;
                    gsk = kdp->gsnp[k][i];
                    gsp = kdp->gsnp[p][i];
                    printf("n=%2d gsnp[%2d]=%2d gsnp[%2d]=%2d\n",
                        n, k, gsk, p, gsp);
                }
            }
            if ((kdp->gsnp[k][i] == 1) && (kdp->gsnp[p][i] == 2)) {
                kdp->parent[k] = n;
                if (p) kdp->parent[n] = p;
            }
        }
    }
    if (debug) debug60(kdp,"find_parents2");
    // fix SNP siblings
    for (k = 0; k < ngen; k++) {
        do {
            j = 0; // number of fixes
            p = kdp->parent[k];
            for (i = 0; i < nsnp; i++) {
                r = snp_present(kdp, kdp->mems[k], kdp->nmem[k], i);
                if (r) {
                    o = kdp->parent[i + ngen];
                    if (o == p) { // fix parent
                        kdp->parent[k] = i + ngen;
                        j++;
                    }
                }
            }
        } while (j);
    }
    // initialize unset GEN parents to most recent ancestral node
    for (k = 1; k < ngen; k++) {
        if (kdp->parent[k] == 0) kdp->parent[k] = kdp->mran;
    }
    kdp->parent[0] = 0;
}

static void
report_desc(struct KITDAT *kdp, int n, int i)
{
    if (debug && 1) {
        int ngen, p, m;
        char *nnm, *pnm;
        p = kdp->parent[n];
        m = kdp->nmem[n];
        ngen = kdp->ngen;
        nnm = (n < ngen) ? kdp->gens[n] : kdp->snps[n - ngen];
        pnm = (p < ngen) ? kdp->gens[p] : kdp->snps[p - ngen];
        
        // FIX: Corrected format string arguments and specifiers.
        // Format string: ("report_desc(%2d): %18s n=%2d m=%2d p=%2d %s\n")
        // Arguments:      (i,   nnm,   n,   m,   p,   pnm)
        printf("report_desc(%2d): %18s n=%2d m=%2d p=%2d %s\n",
            i, nnm, n, m, p, pnm);
    }
}

static void
find_intersect(struct KITDAT *kdp)
{
    char *mn, *mp;
    int ngen, nsnp;
    int i, n, ni, n0, nn, np, p;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    // check for kit intersections
    n0 = kdp->nmem[0];
    for (i = 1; i < (ngen + nsnp); i++) {
        n = i;
        p = kdp->parent[n];
        nn = kdp->nmem[n];
        np = kdp->nmem[p];
        if (nn == n0) continue;
        if (!nn || !np || (nn == np)) continue;
        mn = kdp->mems[n];
        mp = kdp->mems[p];
        ni = intersect(mn, mn, mp, nn, np, 1); // inclusions
        if (debug && 0) {
            report_desc(kdp, i, i);
            printf("%2d %2d %2d %2d %2d: ", i, p, nn, np, ni);
            for (n = 0; n < nn; n++) printf(" %2d", mn[n]);
            printf("\n");
        }
        kdp->nmem[i] = ni;
    }
}

static int
count_duplicates(struct KITDAT *kdp, int node)
{
    char *mems;
    int j, k, ndup, nmem;

    // check for kit intersections
    nmem = kdp->nmem[node];
    mems = kdp->mems[node];
    ndup = 0;
    for (j = 0; j < nmem; j++) {
        for (k = 0; k < j; k++) {
            if (mems[j] == mems[k]) ndup++;
        }
    }
    return (ndup);
}

static int
count_matches(struct KITDAT *kdp, char *mat, int node1, int node2)
{
    char *lfts1, *lfts2;
    int j, k, nmat, nlft1, nlft2;

    // check for kit intersections
    nlft1 = kdp->nlft[node1];
    nlft2 = kdp->nlft[node2];
    lfts1 = kdp->lfts[node1];
    lfts2 = kdp->lfts[node2];
    if (node1 == node2) return (nlft1); // self matches all
    nmat = 0;
    for (j = 0; j < nlft1; j++) {
        for (k = 0; k < nlft2; k++) {
            if (lfts1[j] == lfts2[k]) {
                mat[nmat] = lfts1[j];
                nmat++;
            }
        }
    }
    return (nmat);
}

static void
prune_kits(struct KITDAT *kdp, char *umat, int num)
{
    char *mems, *lfts, xman[MAXGEN];
    int nxma, nxmb, ngen, nsnp, nnod, nmem, nlft, nxmm, nxml;
    int i, j, k, n, kit;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    nxmb = 0;
    for (k = 0; k < num; k++) {
        kit = umat[k];
        // find nodes that contain kit
        nxma = 0;
        for (i = 0; i < nnod; i++) {
            nlft = kdp->nlft[i];
            lfts = kdp->lfts[i];
            nxml = 0;
            for (j = 0; j < nlft; j++) {
                if (lfts[j] == kit) break;
            }
            if (j < nlft) {
                xman[nxma] = i;
                nxma++;
            }
        }
        for (i = 0; i < nxma; i++) {
            if (i == 0) { // retain kit only in first node
                 kdp->xmbn[nxmb] = xman[0];
                 kdp->xmbk[nxmb] = umat[k];
                 nxmb++;
                 kdp->nxmb = nxmb;
                 continue;
            }
            n = xman[i];
            nmem = kdp->nmem[n];
            mems = kdp->mems[n];
            nlft = kdp->nlft[n];
            lfts = kdp->lfts[n];
            nxmm = 0;
            for (j = 0; j < nmem; j++) {
                mems[nxmm] = mems[j];
                if (mems[j] != kit) nxmm++;
            }
            kdp->nmem[n] = nxmm;
            nxml = 0;
            for (j = 0; j < nlft; j++) {
                lfts[nxml] = lfts[j];
                if (lfts[j] != kit) nxml++;
            }
            kdp->nlft[n] = nxml;
        }
    }
}

static void
find_duplicates(struct KITDAT *kdp)
{
    char *plr, *nnm, mat[2048], umat[MAXKIT];
    int i, j, nnod, ndup, ngen, nsnp, nmat, tmat, num;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    // check for kit matches
    nnod = ngen + nsnp;
    tmat = 0;
    for (i = 0; i < nnod; i++) {
        ndup = count_duplicates(kdp, i);
        if (debug && ndup) {
            nnm = (i < ngen) ? kdp->gens[i] : kdp->snps[i - ngen];
            plr = (ndup == 1) ? "" : "s";
            printf("WARNING: %s has %d duplicate%s\n", nnm, ndup, plr);
        }
        for (j = i + 1; j < nnod; j++) {
            nmat = count_matches(kdp, mat + tmat, i, j);
            tmat += nmat;
            if (0 && debug && nmat) {
                nnm = (i < ngen) ? kdp->gens[i] : kdp->snps[i - ngen];
                printf("WARNING: %s ", nnm);
                nnm = (j < ngen) ? kdp->gens[j] : kdp->snps[j - ngen];
                plr = (nmat == 1) ? "" : "es";
                printf("and %s have %d cross match%s\n", nnm, nmat, plr);
            }
        }
    }
    num = unique(mat, umat, tmat, 0);
    prune_kits(kdp, umat, num);
    if (debug && 0) printf("total matches: tmat=%3d num=%2d umat=%2d\n",
        tmat, num, umat[0]);
}

static void
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
            if (kdp->nmem[k] > 0) {
                if (kdp->parent[k] == i) {
                    kdp->kids[i][nkid] = k;
                    nkid++;
                }
            }
        }
        kdp->nkid[i] = nkid;
    }
}

static int
get_descendants(struct KITDAT *kdp, int n, char *dscs, int incmem)
{
    char *mems, *kids, name[40];
    int i, k, nmem, nkid, ndsc;

    if (debug && 0) {
        int ngen, p, m;
        char *nnm, *pnm;
        p = kdp->parent[n];
        m = kdp->nmem[n];
        ngen = kdp->ngen;
        nnm = (n < ngen) ? kdp->gens[n] : kdp->snps[n - ngen];
        pnm = (p < ngen) ? kdp->gens[p] : kdp->snps[p - ngen];
        printf("get_descendants: %18s n=%2d m=%2d p=%2d %s\n",
            nnm, n, m, p, pnm);
    }
    nkid = kdp->nkid[n];
    kids = kdp->kids[n];
    nmem = kdp->nmem[n];
    mems = kdp->mems[n];
    if (incmem) {
        memcpy(dscs, mems, nmem);
        ndsc = nmem;
    } else {
        ndsc = 0;
    }
    for (i = 0; i < nkid; i++) {
        k = get_descendants(kdp, kids[i], dscs + ndsc, 1);
        if (k < 1) {
            node_name(kdp, name, n);
            printf(" ERROR: no descendants: %18s\n", name);
            exit(1);
        }
        ndsc += k;
        if (ndsc > MAXDSC) {
            printf(" ERROR: too many descendants: ndsc=%d\n", ndsc);
            break;
        }
        ndsc = unique(dscs, dscs, ndsc, 0);
    }
    return (ndsc);
}

static void
find_leftovers(struct KITDAT *kdp)
{
    char *mems0, *mems, *lfts, dscs[MAXSTR];
    int i, prnt, ngen, nsnp, nnod, ndsc, nmem0, nmem, nlft;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    nmem0 = kdp->nmem[0];
    mems0 = kdp->mems[0];
    // fix loops
    for (i = 0; i < nnod; i++) {
        prnt = kdp->parent[i];
        if (prnt && (kdp->nmem[prnt] == kdp->nmem[0])) {
            //kdp->parent[i] = 0;
        }
    }
    // find NOD kids
    for (i = 0; i < nnod; i++) {
        nmem = kdp->nmem[i];
        mems = kdp->mems[i];
        lfts = kdp->lfts[i];
        ndsc = get_descendants(kdp, i, dscs, 0);
        ndsc = intersect(dscs, mems0, dscs, nmem0, ndsc, 1);
        nlft = intersect(lfts, mems, dscs, nmem, ndsc, 0);
        kdp->nlft[i] = nlft;
    }
}

static void
rebuild_node(struct KITDAT *kdp, int node)
{
    char *mems, *lfts, *kmems;
    int i, j, kid, nkid, nlft, nmem, knmem;

    nmem = 0;
    mems = kdp->mems[node];
    nkid = kdp->nkid[node];
    for (i = 0; i < nkid; i++) {
        kid = kdp->kids[node][i];
        rebuild_node(kdp, kid);
        knmem = kdp->nmem[kid];
        kmems = kdp->mems[kid];
        for (j = 0; j < knmem; j++) {
            mems[nmem++] = kmems[j];
        }
    }
    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    for (i = 0; i < nlft; i++) {
        mems[nmem++] = lfts[i];
    }
    kdp->nmem[node] = nmem;
}

static void
adjust_snp_age(struct KITDAT *kdp)
{
    char *kids, *snps;
    int i, k, prnt, ysnp, ypar, ykid, kidk, node, nkid, ngen, nsnp;
    static char *fnm = "adjust_snp_age";

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    for (i = 0; i < nsnp; i++) {
        ysnp = kdp->snpy[i];
        snps = kdp->snps[i];
        node = ngen + i;
        if (kdp->nmem[node] == 0) continue;
        prnt = kdp->parent[node];
        ypar = (prnt < ngen) ? kdp->geny[prnt] : kdp->snpy[prnt - ngen]; 
        if (ysnp < ypar) {
            kdp->snpy[i] = ypar;
            if (debug) {
                printf("%s: %4d %4d (par) %s\n", fnm, ysnp, ypar, snps);
            }
        }
        nkid = kdp->nkid[node];
        kids = kdp->kids[node];
        for (k = 0; k < nkid; k++) {
            kidk = kids[k];
            ykid = (kidk < ngen) ? kdp->geny[kidk] : kdp->snpy[kidk - ngen]; 
            if (ykid == 0) continue;
            if (ysnp > ykid) {
                kdp->snpy[i] = ykid;
                if (debug && 0) {
                    printf("%s: %4d %4d (kid) %s\n", fnm, ysnp, ykid, snps);
                }
            }
        }
    }
}

static float
left_age(struct KITDAT *kdp, int node)
{
    float yl;
    int nlft;

    nlft = kdp->nlft[node];
    yl = nlft ? 1960 : 0;
    return (yl);
}

static float
kids_age(struct KITDAT *kdp, int node)
{
    float sy, yl, yk;
    int k, kidk, nk, nkit, nkid;

    nkit = kdp->nkit;
    nkid = kdp->nkid[node];
    if (!nkid) return (1960);
    sy = 0;
    nk = 0;
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        yk = kdp->kity[nkit + kidk];
        if (yk == 0) {
            yl = left_age(kdp, kidk);
            yk = kids_age(kdp, kidk);
            if (yl && yk) {
                yk = (yl + yk) / 2;
            } else if (yl) {
                yk = yl;
            }
        }
        sy += yk;
        nk++;
    }
    yk = nk ? sy / nk : 0;
    return (yk);
}

static float
kids_sig(struct KITDAT *kdp, int node)
{
    float nsig;
    int k, kidk, nkid;

    nkid = kdp->nkid[node];
    if (nkid < 1) return (0);
    nsig = 0;
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        nsig += kdp->sig[kidk].n;
    }
    return (nsig / nkid); 
}

static void
date_str(struct KITDAT *kdp)
{
    float yk, ksig, psig;
    int i, isbg, yp, prnt, sbgi, ngen, nsnp, nkit, nsbg;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    nsbg = kdp->nsbg;
    isbg = ngen + nsnp;
    for (i = 0; i < nsbg; i++) {
        sbgi = isbg + i;
        if (isalpha(kdp->kits[nkit + sbgi][0])) {
            prnt = kdp->parent[sbgi];
            yp = kdp->kity[nkit + prnt];
            yk = kids_age(kdp, sbgi);
            psig = 1 + kdp->sig[prnt].n;
            ksig = 1 + kids_sig(kdp, sbgi);
            kdp->kity[nkit + sbgi] = yp + (yk - yp) * psig / (psig + ksig);
        }
        if (kdp->nkid[sbgi] == 0) {
            kdp->nmem[sbgi] = kdp->nlft[sbgi]; // patch nmem
        }
    }
}

static void
fix_gensnp(struct KITDAT *kdp)
{
    int i, j, k, ngen, nsnp;

    // eliminate uncertains
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    for (i = 0; i < ngen; i++) {
        for (j = 0; j < nsnp; j++) {
            if (kdp->gsnp[i][j] == 0) {
                for (k = 0; k < nsnp; k++) {
                    if (kdp->gsnp[i][k] == 1) {
                        kdp->gsnp[i][j] = 2;
                    }
                }
            }
        }
    }
}

static void
merge_gensnp(struct KITDAT *kdp)
{
    char *s, *d, year[8];
    int i, n, y, sumlft, ngen, nsnp, mxnl, nnod;

    if (debug) printf("->merge_gensnp\n");
    excl_gen(kdp);
    snp_logic(kdp); // <-- DEBUG
    init_modal(kdp);
    find_members(kdp);
    snp_fill(kdp); // <-- DEBUG
    find_parents(kdp);
    find_intersect(kdp);
    find_kids(kdp);
    find_leftovers(kdp);
    find_duplicates(kdp); 
    rebuild_node(kdp, 0); 
    adjust_snp_age(kdp);
    fix_gensnp(kdp);
    // report
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nnod = ngen + nsnp;
    if (debug) {
        // find max name length
        mxnl = 0;
        for (i = 0; i < ngen; i++) {
            n = strlen(kdp->gens[i]);
            if (mxnl < n) mxnl = n;
        }
        for (i = 0; i < nsnp; i++) {
            n = strlen(kdp->snps[i]);
            if (mxnl < n) mxnl = n;
        }
        // debug print
        sumlft = 0;
        for (i = 0; i < nnod; i++) {
            s = (i < ngen) ? kdp->gens[i] : kdp->snps[i - ngen];
            y = (i < ngen) ? kdp->geny[i] : kdp->snpy[i - ngen];
            sprintf(year, ".%d", y);
            d = y ? year : "     ";
            printf("%2d %*s%s  ", i, mxnl, s, d);
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

static void
prepare_strs(struct KITDAT *kdp)
{
    char sv[MAXSTR];
    int i, j, nkit, ngen, nsnp, nstr;

    if (debug) printf("->prepare_strs: nmem=%d\n", kdp->nmem[0]);
    // initialize STR names
    kdp->nlab = get_strlab(kdp->labs);
    // patch missing STR values
    patch_strs(kdp, 0);
    // count missing STR values
    kdp->tnmv  = find_misval(kdp);
    // initialize STR sort order (by prevalence)
    str_sort(kdp, kdp->sso, 1);
    // initialize node STRs
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit++; // skip ancestral node
    for (i = 1; i < ngen; i++) { 
        nstr = modal_set(kdp, sv, kdp->mems[i], kdp->nmem[i]);
        add_kit(kdp, kdp->gens[i], kdp->geny[i], sv, nstr);
    }
    for (j = 0; j < nsnp; j++) { 
        i = j + ngen;
        nstr = modal_set(kdp, sv, kdp->mems[i], kdp->nmem[i]);
        add_kit(kdp, kdp->snps[j], kdp->snpy[j], sv, nstr);
    }
    kdp->nkit = nkit; // restore kit count
}

static int
isleftover(struct KITDAT *kdp, int kit, int node)
{
    char *lfts;
    int i, nlft;

    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    for (i = 0; i < nlft; i++) if (lfts[i] == kit) break;
    return (i < nlft);
}

static int
isdescendant(struct KITDAT *kdp, int n1, int n2)
{
    int p1;

    // did n1 descend from n2 ?
    p1 = kdp->parent[n1];
    if (n2 ==  0) return (1);
    if (p1 ==  0) return (0);
    if (n1 == n2) return (1);
    if (p1 == n2) return (1);
    return (isdescendant(kdp, p1, n2));
}

static int
get_strmut(struct KITDAT *kdp, int node, short si, char *sp, char *sm)
{
    int k, kidk, pv, mv, prnt, nkid;

    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        if (get_strmut(kdp, kidk, si, sp, sm)) return (1);
    }
    prnt = kdp->parent[node];
    pv = get_strn(prnt, si);
    mv = get_strn(node, si);
    if (mv != pv) {
        if (sp) sp[0] = pv;
        if (sm) sm[0] = mv;
        return (1);
    }
    return (0);
}

static void
check_kids(struct KITDAT *kdp)
{
   int k, n, kidk, nkid, nmkd, ntot, nmem, ngen, nsnp, nsbg, nnod;

   ngen = kdp->ngen;
   nsnp = kdp->nsnp;
   nsbg = kdp->nsbg;
   nnod = ngen + nsnp + nsbg;
   for (n = 0; n < nnod; n++) {
       nmem = kdp->nmem[n];
       nkid = kdp->nkid[n];
       ntot = kdp->nlft[n];
       for (k = 0; k < nkid; k++) {
           kidk = kdp->kids[n][k];
           nmkd = kdp->nmem[kidk];
           ntot += nmkd;
           if (ntot > nmem) {
               printf("check_kids: "); 
               printf("nt[%d/%d]=%d > nm[%d]=%d; ", k, nkid, ntot, n, nmem);
               printf("nm[%d]=%d prnt=%d\n", kidk, nmkd, kdp->parent[kidk]);
              break;
           }
       }
   }
}

static void
copy_strs(struct KITDAT *kdp, int node)
{
    int i, k, nk, pk, prnt, nlft, nkid, ngen, nstr, nkit;
    static char *fnm = "copy_strs";

    nkid = kdp->nkid[node];
    // copy STRs from parent
    if (node > 0) {
        ngen = kdp->ngen;
        nlft = kdp->nlft[node];
        if ((nlft == 0)  && (nkid == 1) && (node >= ngen)) {
            nstr = kdp->nstr;
            nkit = kdp->nkit;
            prnt = kdp->parent[node];
            nk = nkit + node;
            pk = nkit + prnt;
            memcpy(kdp->kstrs[nk], kdp->kstrs[pk], nstr);
            if (debug) {
                char nnm[40];
                node_name(kdp, nnm, node);
                printf("%s: %13s\n", fnm, nnm);
            }
        }
    }
    // copy STRs from parents to kids
    for (i = 0; i < nkid; i++) {
        k = kdp->kids[node][i];
        copy_strs(kdp, k);
    }
}

static void
infer_strs(struct KITDAT *kdp, int node)
{
    char *mems, sv[MAXSTR];
    int i, k, ngen, nkid, nmem, nstr;

    ngen = kdp->ngen;
    if (node >= ngen) return;
    // infer STRs of kids
    nkid = kdp->nkid[node];
    for (i = 0; i < nkid; i++) {
        k = kdp->kids[node][i];
        infer_strs(kdp, k);
    }
    if ((node == 0) || (nkid > 0)) return;
    // infer STRs from leftovers
    nmem = kdp->nmem[node];
    mems = kdp->mems[node];
    if (nmem == 0) return;
    nstr = modal_set(kdp, sv, mems, nmem);
    for (i = 0; i < nstr; i++) {
        if (sv[i]) set_strn(node,i,sv[i]);
    }
}

static void
clear_strrev(struct KITDAT *kdp)
{
    int nkit, ngen, nsnp, nsbg, nnod, nstr;
    int i, n, k, mv, pv, gv, node, prnt;

    nkit = kdp->nkit;
    nstr = kdp->nstr;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    for (k = 0; k < nkit; k++) {
        for (n = 0; n < nnod; n++) {
            if (isleftover(kdp, k, n)) break;
        }
        if (n == nnod) {
            printf("WARNING: kit %d has no terminal node\n", k);
            continue;
        }
        if ((n < ngen) && (kdp->parent[n] > ngen)) continue;
        for (i = 0; i < nstr; i++) {
            mv = get_stro(k, i);
            pv = get_strn(n,i);
            if (!mv) continue;
            node = n;
            prnt = kdp->parent[n];
            while (prnt) {
                gv = get_strn(prnt,i);
                if ((mv != pv) && (mv == gv)) {
                    set_strn(node, i, mv);
                }
                pv = get_strn(prnt,i);
                node = prnt;
                prnt = kdp->parent[node];
            }
        }
    }
}

static void
elevate_sbgs(struct KITDAT *kdp)
{
    int ngen, nsnp, nsbg, nmem, nkit;
    int i, k, n, isbg, sbgi, prnt;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    isbg = ngen + nsnp;
    for (i = 0; i < nsbg; i++) {
        sbgi = isbg + i;
        nmem = kdp->nmem[sbgi];
        prnt = kdp->parent[sbgi];
        if (nmem == 0) continue;
        if (kdp->nkid[sbgi]) continue;
        if (kdp->nlft[prnt]) continue;
        memcpy(kdp->lfts[prnt], kdp->mems[sbgi], nmem);
        kdp->nlft[prnt] = nmem;
        if (debug) {
            char nnm[40], pnm[40];
            name_branch(kdp, sbgi);
            node_name(kdp, nnm, sbgi);
            node_name(kdp, pnm, prnt);
            printf("%s > %s\n", pnm, nnm);
        }
        kdp->nlft[sbgi] = 0;
        kdp->nmem[sbgi] = 0;
        kdp->nkid[prnt]--; // ???
    }
    nkit = kdp->nkit;
    k = nkit + sbgi;
    n = 0;
    for (i = 0; i < nsbg; i++) {
        sbgi = isbg + i;
        nmem = kdp->nmem[sbgi];
        prnt = kdp->parent[sbgi];
        if (nmem == 0) continue;
        strcpy(kdp->sbgs[n], kdp->sbgs[sbgi]);
        memcpy(kdp->mems[n], kdp->mems[sbgi], nmem);
        memcpy(kdp->lfts[n], kdp->mems[sbgi], nmem);
        kdp->nlft[n] = nmem;
        kdp->nmem[n] = nmem;
        kdp->parent[n] = prnt;
        memcpy(kdp->kstrs[nkit + n], kdp->kstrs[k], kdp->nsvs[k]);
        kdp->nsvs[nkit + n] = kdp->nsvs[k];
        n++;
        kdp->nkid[prnt]--;
    }
    if (debug && n) {
        printf("elevate_sbgs: isbg=%d nsbg=%d->%d\n",
            isbg, kdp->nsbg, n);
    }
    kdp->nsbg = n;
}

static int
get_usig(struct KITDAT *kdp, int node, struct SIGDAT *usig, int nusg)
{
    int i, j, k, s, p, m, kidk, ngkd, nsig;
    int nkid = kdp->nkid[node]; 

    if ((nusg == MAXSIG) || (nkid == 0)) return(nusg);
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        ngkd = kdp->kids[kidk][0];
        nsig = kdp->sig[kidk].n;
        if ((nsig == 0) && (ngkd > 0)) {
            kidk = kdp->kids[kidk][0];
            nsig = kdp->sig[kidk].n;
        }
        for (i = 0; i < nsig; i++) {
            s = kdp->sig[kidk].s[i];
            p = kdp->sig[kidk].p[i];
            m = kdp->sig[kidk].m[i];
            for (j = 0; j < nusg; j++) {
               if ((usig->s[j]== s) &&
                   (usig->p[j]== p) &&
                   (usig->m[j]== m)) break;
            }
            if (j == nusg) {
                usig->s[j] = s;
                usig->p[j] = p;
                usig->m[j] = m;
                usig->c[j] = 0;
                nusg++;
            }
            usig->c[j]++;
        }
     }
    return (nusg);
}

static void
tran_kids(struct KITDAT *kdp, int node, int grpn, struct SIGN sig)
{
    char sp[1], sm[1];
    int k, m, n, nm, kidk, nkid;
    short ki[MAXGEN];

    nkid = kdp->nkid[node];
    // select kids for transfer to grpn
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        ki[k] = get_strmut(kdp, kidk, sig.s, sp, sm);
    }
    ki[nkid - 1] = 0;
    // transfer kids from node to grpn
    n = m = 0;
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        if (ki[k]) {
            kdp->kids[grpn][n] = kidk;
            kdp->parent[kidk] = grpn;
            n++;
            nm = kdp->nmem[kidk];
            memcpy(kdp->mems[grpn] + m, kdp->mems[kidk], nm);
            m += nm;
            get_strmut(kdp, kidk, sig.s, sp, sm);
        } else {
            kdp->kids[node][k - n] = kidk;
            kdp->parent[kidk] = node;
        }
    }
    kdp->nmem[grpn] = m; 
    kdp->nkid[grpn] = n;
    kdp->nkid[node] -= n;
    check_kids(kdp);
    // set grpn STR values
    set_strn(grpn, sig.s, sig.m);
    nkid = kdp->nkid[grpn];
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[grpn][k];
        set_strn(kidk, sig.s, sig.m);
    }
}

static int
group_kids(struct KITDAT *kdp)
{
    int ngen, nsnp, nsbg, nnod, nusg, nkid;
    int s, p, m, c, i, n, grpn;
    struct SIGN sig;
    struct SIGDAT usig;
    static char *set = "";

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    grpn = 0;
    for (n = 0; n < nnod; n++) {
        update_sig(kdp, n);
        nusg = get_usig(kdp, n, &usig, 0);
        for (i = 0; i < nusg; i++) {
            s = usig.s[i];
            p = usig.p[i];
            m = usig.m[i];
            c = usig.c[i];
            nkid = kdp->nkid[n];
            if (c > 1) {
                if ((3 * c) > (2 * nkid)) {
                    if (grpn) {
                        set_strn(grpn,s,m);
                        update_sig(kdp, grpn);
                        set = "<--";
                    } else {
                        set_strn(n,s,m);
                        update_sig(kdp, n);
                        set = "<-";
                    }
                } else {
                    sig.s = s;
                    sig.p = p;
                    sig.m = m;
                    if (grpn) {
                        set_strn(n,s,m);
                        set_strn(grpn,s,m);
                        update_sig(kdp, grpn);
                    } else {
                        grpn = add_branch(kdp, n);
                        name_branch(kdp, grpn);
                        tran_kids(kdp, n, grpn, sig);
                    }
                    set = "++";
                }
                if (debug) {
                    char nnm[20], gnm[20], *str;
	            node_name(kdp, nnm, n);
                    node_name(kdp, gnm, grpn);
                    str = kdp->labs[s];
                    printf("group_kids: ");
                    printf("%-18s>%s ", nnm, gnm);
                    printf("%s=%d->%d|%d %s ", str, p, m, c, set);
                    printf("nkid=%d nusg=%d ", nkid, nusg);
	            printf("\n");
                }
            }
        }
    }
    update_sig(kdp, kdp->mran);
    return (grpn);
}

static void
tran_left(struct KITDAT *kdp, int node, int grpn, char *lfts, char *kids, 
    int mlft, int mkid)
{
    char *memsg, *memsk;
    int k, m, n, kidk, kidm, lftk, nkid, nlft, nmemg, nmemk;

    if (!mlft) return;
    // add kits to STR node
    // delete kits from leftovers
    nlft = kdp->nlft[node];
    n = m = 0;
    for (k = 0; k < nlft; k++) {
        lftk = kdp->lfts[node][k];
        if (memchr(lfts, lftk, mlft)) {
            kdp->lfts[grpn][m] = lftk;
            kdp->mems[grpn][m] = lftk;
            m++;
        } else {
            kdp->lfts[node][n] = lftk;
            n++;
        }
    }
    kdp->nlft[grpn] = m;
    kdp->nmem[grpn] = m;
    kdp->nlft[node] = n;
    if (mkid) {
        // add kids to STR node
        // delete kids from parent node
        nkid = kdp->nkid[node];
        n = m = 0;
        for (k = 0; k < nkid; k++) {
            kidk = kdp->kids[node][k];
            if (memchr(kids, kidk, mkid)) {
                kidm = kidk;
                kdp->kids[grpn][m] = kidm;
                kdp->parent[kidm] = grpn;
                nmemg = kdp->nmem[grpn];
                nmemk = kdp->nmem[kidk];
                memsg = kdp->mems[grpn] + nmemg;
                memsk = kdp->mems[kidk];
                memcpy(memsg, memsk, nmemk);
                kdp->nmem[grpn] = nmemg + nmemk;
                m++;
            } else {
                kdp->kids[node][n] = kidk;
                n++;
            }
        }
        kdp->nkid[grpn] = m;
        kdp->nkid[node] = n;
    }
}

/**
 * Searches for leftovers and kids in a node that share a co-mutation.
 * This is the core logic extracted from brnch_lftovr.
 *
 * @param node The current node index.
 * @param i    The STR marker index to check.
 * @param v    The node's current STR value (potential ancestral value).
 * @param lfts_out Buffer to store leftover kit indices that match the mutation.
 * @param mlft_out Pointer to store the count of matching leftovers.
 * @param kids_out Buffer to store child node indices that match the mutation.
 * @param mkid_out Pointer to store the count of matching kids.
 */
void
find_comut_kits(struct KITDAT *kdp, int node, int i, int v,
    char *lfts_out, int *mlft_out, char *kids_out, int *mkid_out)
{
    int j, k, lftk, kidk, w, x;
    int nlft = kdp->nlft[node];
    int nkid = kdp->nkid[node];
    int mlft = 0;
    int mkid = 0;
    
    // Temporary buffers for gathering kits that share the same mutation
    char temp_lfts[MAXKIT]; 
    int temp_mlft = 0;
    int target_w = 0;

    // Iterate through all leftovers to find a base kit that supports a mutation (w != v)
    for (j = 0; j < nlft; j++) {
        lftk = kdp->lfts[node][j];
        w = get_stro(lftk, i); // Kit's STR value

        if (w && (v != w)) {
            // Found a kit with a mutation 'w'. Now look for support (co-mutations).
            target_w = w;
            temp_mlft = 0;
            temp_lfts[temp_mlft++] = lftk; // Add the base kit
            
            // 1. Check other leftover kits (k > j) for the same mutation 'w'
            for (k = j + 1; k < nlft; k++) {
                // Cast kdp->lfts[node][k] to int to avoid char subscript warning
                if (get_stro((int)kdp->lfts[node][k], i) == target_w) {
                    temp_lfts[temp_mlft++] = kdp->lfts[node][k];
                }
            }

            // 2. Check child nodes for the same modal mutation 'w'
            for (k = 0; k < nkid; k++) {
                kidk = kdp->kids[node][k];
                x = get_strn(kidk, i); // Child node's STR modal value

                if (x == target_w) {
                    kids_out[mkid++] = kidk;
                }
            }

            // If we found any support (more than 1 kit, or any supporting kids)
            if (temp_mlft > 1 || mkid > 0) {
                // Transfer temporary leftovers (including the base kit lftk) to the final output buffer
                for(x = 0; x < temp_mlft; x++) {
                    lfts_out[mlft++] = temp_lfts[x];
                }
                break; // Exit the main leftover loop
            }
            // If only the base kit was found and no kids supported it, continue search.
            mkid = 0; // Reset kid count as this mutation didn't yield a consensus branch
        }
    }
    
    // Clean duplicates (just in case the complex looping introduced any)
    mlft = unique(lfts_out, lfts_out, mlft, 0);

    *mlft_out = mlft;
    *mkid_out = mkid;
}

static int
brnch_lftovr(struct KITDAT *kdp)
{
    char lfts[MAXKIT], kids[MAXGEN];
    int ngen, nsnp, nsbg, nnod, nlft, nkid, nstr;
    int i, grpn, mlft, mkid, n, node, v, w;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    nstr = kdp->nstr;
    grpn = 0;

    for (n = 0; n < nnod; n++) {
        node = n;
        nlft = kdp->nlft[node];
        nkid = kdp->nkid[node];
        
        i = nstr; // Default to an invalid index
        w = 0;
        
        // Iterate through all STR markers to find a co-mutation pattern
        for (int s = 0; s < nstr; s++) {
            v = get_strn(node, s); // Node's STR value (ancestral/modal)
            
            // Skip if node value is missing
            if (v == 0) continue;

            // Use the helper function to find co-mutations for this marker 's'
            find_comut_kits(kdp, node, s, v, lfts, &mlft, kids, &mkid);

            if (mlft || mkid) {
                // Found a co-mutation.
                // We need the mutated value 'w' from one of the kits/kids.
                if (mlft) {
                    w = get_stro((int)lfts[0], s); // Cast to int to avoid char subscript warning
                } else { // Must have mkid > 0
                    w = get_strn(kids[0], s);
                }
                
                i = s; // Store the successful STR index
                break; // Stop searching for markers
            }
        }
        
        if (i >= nstr) continue; // No co-mutations found for this node.
        
        if (mkid && (mkid == nkid)) { // If all kids carry the mutation, elevate the STR mutation to the node
            set_strn(node, i, w);
        } else if (mkid || mlft) {    // Create new STR branch (SBG)
            
            grpn = add_branch(kdp, node);
            name_branch(kdp, grpn);
            tran_left(kdp, n, grpn, lfts, kids, mlft, mkid);
            set_strn(node, i, v);
            set_strn(grpn, i, w);
            update_sig(kdp, node);
	    if (debug) {
                char nnm[20], gnm[20], *str;
	        node_name(kdp, nnm, n);
                node_name(kdp, gnm, grpn);
                str = kdp->labs[i];
	        printf("brnch_lftovr: %18s > %s ", nnm, gnm);
	        if (mlft) printf("%s=%d->%d ", str, v, w);
                printf("mlft=%d/%d ", mlft, nlft);
                printf("mkid=%d/%d ", mkid, nkid);
	        printf("nsig=%d\n", kdp->sig[grpn].n);
	    }
	}
    }
    update_sig(kdp, kdp->mran);
    return (grpn);
}

static int
gendis(struct KITDAT *kdp, int kit1, int kit2)
{
    int i, gd, sv1, sv2, nsvs;

    sv1 = kdp->nsvs[kit1];
    sv2 = kdp->nsvs[kit2];
    nsvs = (sv1 < sv2) ? sv1 : sv2;
    gd = 0;
    for (i = 0; i < nsvs; i++) {
        sv1 = kdp->kstrs[kit1][i];
        sv2 = kdp->kstrs[kit2][i];
        if (sv1 && sv2 && (sv1 != sv2)) gd++;
    }
    return (gd);
}

static void
raise_orphans(struct KITDAT *kdp, int node)
{
    char *kids, *lfts, kmov[MAXKIT];
    int ngen, nsnp, nkid, nlft, nkit, nmov;
    int i, j, k, m, n, lftk, kiti, gdn, gdp, mlft, prnt;
    static char *fnm = "raise_orphans";

    nkid = kdp->nkid[node];
    kids = kdp->kids[node];
    for (k = 0; k < nkid; k++) {
        raise_orphans(kdp, kids[k]);
    }
    prnt = kdp->parent[node];
    nlft = kdp->nlft[node];
    lfts = kdp->lfts[node];
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    nmov = 0;
    for (i = 0; i < nlft; i++) {
        kiti = lfts[i];
        for (j = 1; j < ngen; j++) {
            if (kdp->kgen[kiti][j]) break;
        }
        for (k = 0; k < nsnp; k++) {
            if (kdp->ksnp[kiti][k]) break;
        }
        gdn = gendis(kdp, kiti, nkit + node);
        gdp = gendis(kdp, kiti, nkit + prnt);
        if ((j == ngen) && (k == nsnp) && (gdn >= gdp)) {
            kmov[nmov++] = kiti;
            if (debug) {
                //char *kits;
                //kits = kdp->kits[kiti];
                //printf("%s: %s gdn=%d gdp=%d\n", fnm, kits, gdn, gdp);
            }
        }
    }
    if (nmov && isdescendant(kdp, node, kdp->mran)) {
        // move kits from node to prnt;
        prnt = kdp->parent[node];
        mlft = kdp->nlft[prnt];
        nlft = kdp->nlft[node];
        lfts = kdp->lfts[node];
        m = n = 0;
        for (k = 0; k < nlft; k++) {
            lftk = lfts[k];
            if (memchr(kmov, lftk, nmov)) {
                kdp->lfts[prnt][m + mlft] = lftk;
                m++;
            } else {
                kdp->lfts[node][n] = lftk;
                kdp->mems[node][n] = lftk;
                n++;
            }
        }
        kdp->nlft[prnt] += m;
        kdp->nlft[node] = n;
        kdp->nmem[node] = n;
        if (debug) {
            char nnm[20];
            node_name(kdp, nnm, node);
            printf("%s: %s nmov=%d\n", fnm, nnm, nmov);
        }
    }
}

static void
set_mran(struct KITDAT *kdp)
{
    int av, i, node, prnt, nstr;

    node = kdp->mran;
    nstr = kdp->nstr;
    while (node) {
        prnt = kdp->parent[node];
        for (i = 0; i < nstr; i++) {
            av = get_strn(node, i);
            set_strn(prnt, i, av);
        }
        node = prnt;
    }
    update_sig(kdp, kdp->mran);
}

static void
set_term(struct KITDAT *kdp, int node)
{
    int k, kitk, kidk, nlft, nkid;

    nlft = kdp->nlft[node];
    nkid = kdp->nkid[node];
    for (k = 0; k < nlft; k++) {
        kitk = kdp->lfts[node][k];
        kdp->term[kitk] = node;
    }
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        set_term(kdp, kidk);
    }
}

int
acceptable_kitnod(struct KITDAT *kdp, int kit, int node)
{
    char *kgen, *ksnp, *gsnp;
    int i, j, g, si, sj, ngen, nsnp;

    // check whether node is unused
    if (kdp->nmem[node] == 0) return (0);
    // check GEN constraints
    ngen = kdp->ngen;
    kgen = kdp->kgen[kit];
    for (g = 1; g < ngen; g++) {
        if(kgen[g] == 1) { // include
            if (!ismember(kdp, kit, g)) return (0);
            if (!isdescendant(kdp, node, g)) return (0);
        }
        if(kgen[g] == 2) { // exclude
            if (isdescendant(kdp, node, g)) return (0);
        }
    }
    // check SNP constraints
    nsnp = kdp->nsnp;
    ksnp = kdp->ksnp[kit];
    for (i = 0; i < nsnp; i++) {
        si = ngen + i; // SNP node i
        if(ksnp[i] == 1) { // include
            if (!isdescendant(kdp, node, si)) return (0);
            for (j = 0; j < nsnp; j++) {
                sj = ngen + j; // SNP node j
                if (j == i) continue;
                if(ksnp[j] == 1) { // include
                    if (isdescendant(kdp, sj, si)) return (0);
                }
            }
        }
        if(ksnp[i] == 2) { // exclude
            if (isdescendant(kdp, node, si)) return (0);
        }
    }
    // check GEN>SNP constraints
    nsnp = kdp->nsnp;
    for (i = 1; i < ngen; i++) {
        gsnp = kdp->gsnp[i];
        for (j = 0; j < nsnp; j++) {
            sj = ngen + j; // SNP node j
            if((gsnp[j] == 1) || (gsnp[j] == 3)) { // include
                if (!ismember(kdp, kit, sj)) return (0);
                if (!isdescendant(kdp, node, sj)) return (0);
            }
            if((gsnp[j] == 2) || (gsnp[j] == 4)) { // exclude
                if (isdescendant(kdp, node, sj)) return (0);
            }
        }
    }
    return (1);
}

static int
min_gd(struct KITDAT *kdp, int kit, int node, int *mngd)
{
    int gdmn, gdnd, k, kidk, nkit, nkid;

    nkit = kdp->nkit;
    nkid = kdp->nkid[node];
    gdnd = gendis(kdp, kit, nkit + node);
    gdmn = gendis(kdp, kit, nkit + *mngd);
    if (acceptable_kitnod(kdp, kit, node) && (gdmn > gdnd)) {
        gdmn = gdnd;
        *mngd = node;
    }
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        gdmn = min_gd(kdp, kit, kidk, mngd);
    }
    return (gdmn);
}

static void
move_kit(struct KITDAT *kdp, int kit, int mngd, int term)
{
    int k, kitk, n, nlft, nmem, node;

    if (debug) {
        printf("move_kit(%7s): ", kdp->kits[kit]);
        printf("%2d <- %2d\n", mngd, term);
    }
    // delete kit from term node
    node = term;
    nlft = kdp->nlft[node];
    n = 0;
    for (k = 0; k < nlft; k++) {
        kitk = kdp->lfts[node][k];
        kdp->lfts[node][n] = kitk;
        if (kitk != kit) n++; 
    }
    kdp->nlft[node] = n;
    while (node) {
        nmem = kdp->nmem[node];
        n = 0;
        for (k = 0; k < nmem; k++) {
            kitk = kdp->mems[node][k];
            kdp->mems[node][n] = kitk;
            if (kitk != kit) n++; 
        }
        kdp->nmem[node] = n;
        node = kdp->parent[node];
    }
    // add kit to mngd node
    node = mngd;
    nlft = kdp->nlft[node];
    kdp->lfts[node][nlft] = kit;
    kdp->nlft[node]++;
    while (node) {
        nmem = kdp->nmem[node];
        kdp->mems[node][nmem] = kit;
        kdp->nmem[node]++;
        node = kdp->parent[node];
    }
}

static int
check_kits(struct KITDAT *kdp)
{
    int gdmn, gdnd, k, mngd, nmov, nkit, term;

    //return (0); // disable moves until acceptability is fixed <--

    // initialize all mngd & term
    nkit = kdp->nkit;
    for (k = 0; k < nkit; k++) {
        kdp->term[k] = 0;
    }
    // set term while descending all nodes
    set_term(kdp, 0);
    // report displaced nodes
    nmov = 0;
    for (k = 0; k < nkit; k++) {
        term = kdp->term[k];
        mngd = term;
        gdmn = min_gd(kdp, k, 0, &mngd);
        gdnd = gendis(kdp, k, nkit + term);
        if (acceptable_kitnod(kdp, k, mngd) && (gdmn < gdnd)) {
            move_kit(kdp, k, mngd, term);
            nmov++;
        }
    }
    if (nmov || debug) {
        printf("check_kits: %d moved\n", nmov);
    }
    if (nmov) {
        update_sig(kdp, kdp->mran);
    }
    return (nmov);
}

static void
branch_strs(struct KITDAT *kdp)
{
    if (debug) printf("-> branch_strs: nmem=%d\n", kdp->nmem[0]);
    elevate_sbgs(kdp);
    copy_strs(kdp, 1);
    infer_strs(kdp, 1);
    raise_orphans(kdp, 0);
    clear_strrev(kdp);
    set_mran(kdp);
    group_kids(kdp);
    brnch_lftovr(kdp);
    check_kits(kdp);
    update_sig(kdp, kdp->mran);
    date_str(kdp);
    update_sig(kdp, 0);
}

static void
inherit_strs(struct KITDAT *kdp)
{
	    int i, nnod;

    if (debug) printf("->inherit_strs: nmem=%d\n", kdp->nmem[0]);
    // clear list of STR signature mutations
    nnod = kdp->ngen + kdp->nsnp;
    for (i = 0; i < nnod; i++) {
        kdp->sig[i].n = 0;
    }
    // being inheritance at node=0
    inherit(kdp, 0);
}

static struct KITDAT *
get_sapp(char *ifn, int group)
{
    char *s, line[MAXLIN];
    FILE *ifp;
    static struct KITDAT kd;

    // initialize GEN SNP STR data
    kd.nkit = kd.nstr = kd.ngen =kd.nsnp =  kd.nsbg = 0;
    init_gendat(&kd, group);
    ifp = fopen(ifn, "r"); // input file
    if (ifp == NULL) {
        printf("ERROR: can't open %s\n", ifn);
        exit(1);
    }
    while (fgets(line, MAXLIN, ifp)) {
        if (line[0] == EOF) break;
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
            } else if (strncmp(line, "/GROUPS", 7) == 0) {
                get_grpdata(&kd, ifp, line);
            } else {
                while ((s = fgets(line, MAXLIN, ifp))) {
                    if (line[0] == '/') break;
                }
                if (s == NULL) { line[0] = EOF; break; }
            }
        }
    }
    fclose(ifp);
    // get SNP year & parent
    fetch_snptree(&kd);
    // process GEN & SNP data
    merge_gensnp(&kd);
    // process STR data
    prepare_strs(&kd);
    inherit_strs(&kd);
    branch_strs(&kd);
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

static char *
d(int n)
{
    static char *s[5] = {"?", "+", "-", "?+", "?-"};

   return (s[n]);
}

static void
kit_gensnp_csv(struct KITDAT *kdp, int group)
{
    char ofn[MAXFNL];
    int i, j, nsnp, ngen;
    FILE *ofp;
    static char *name = "kit_gensnp";

    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" %s", ofn);
    nsnp = kdp->nsnp;
    ngen = kdp->ngen;
    ofp = fopen(ofn, "w");
    fputs("kit,", ofp);
    for (i = 1; i < ngen; i++) fprintf(ofp, "%s,", kdp->gens[i]);
    for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", kdp->snps[i]);
    fputs("\n", ofp);
    for (j = 0; j < kdp->nkit; j++) {
        fprintf(ofp,"%s,", kdp->kits[j]);
        for (i = 1; i < ngen; i++) fprintf(ofp, "%s,", d(kdp->kgen[j][i]));
        for (i = 0; i < nsnp; i++) fprintf(ofp, "%s,", d(kdp->ksnp[j][i]));
        fputs("\n", ofp);
    }
    fclose(ofp);
}

static void
gen_snp_csv(struct KITDAT *kdp, int group)
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
    printf(" %s", ofn);
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

static void
name_branches(struct KITDAT *kdp)
{
    int aka, i, ngen, nsnp, nsbg, nkit, isbg, sbgi;

    // name STR branches
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    isbg = ngen + nsnp;
    nkit = kdp->nkit;
    aka = 1960;
    for (i = 0; i < nsbg; i++) {
        sbgi = isbg + i;
        name_branch(kdp, sbgi);
        if (kdp->kity[nkit + sbgi] == 0) {
            kdp->kity[nkit + sbgi] = aka - 20;
        }
    }
}

static void
str_all_csv(struct KITDAT *kdp, int group)
{
    char ofn[MAXFNL], nodnam[40];
    int i, j, n, nstr, nkit, year, val;
    FILE *ofp;
    static char *name = "sv_all";

    // output all STRs
    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" %s", ofn);
    nstr = kdp->nstr;
    ofp = fopen(ofn, "w");
    fputs("kit", ofp);
    for (i = 0; i < nstr; i++) {
        fputs(",", ofp);
        fputs(kdp->labs[i], ofp);
    }
    fputs("\n", ofp);
    nkit  = kdp->nkit;
    n = kdp->nkit + kdp->ngen + kdp->nsnp + kdp->nsbg;
    for (j = 0; j < n; j++) {
        year = kdp->kity[j];
        if ((j < nkit)  || (year == 0)) {
            sprintf(nodnam,"%s,", kdp->kits[j]);
        } else {
            sprintf(nodnam,"%s.%d,", kdp->kits[j], year);
        }
        fprintf(ofp,"%s", nodnam);
        for (i = 0; i < nstr; i++) {
            val = get_stro(j,i);
            fprintf(ofp, "%d", val);
            if (i < (nstr - 1)) {
                fputs(",", ofp);
            }
        }
        fputs("\n", ofp);
    }
    fclose(ofp);
}

static int
select_mems(struct KITDAT *kdp, char *mems, int node)
{
    int nkit, nmem, prnt;

    nkit = kdp->nkit;
    nmem = kdp->nmem[node];
    prnt = kdp->parent[node];
    memcpy(mems, kdp->mems[node], nmem);
    mems[nmem++] = nkit + node; // append node
    mems[nmem++] = nkit + prnt; // append prnt
    return (nmem);
}

static int
select_lfts(struct KITDAT *kdp, char *lfts, int node)
{
    int nkit, nlft, prnt;

    nkit = kdp->nkit;
    nlft = kdp->nlft[node];
    prnt = kdp->parent[node];
    memcpy(lfts, kdp->lfts[node], nlft);
    lfts[nlft++] = nkit + node; // append node
    lfts[nlft++] = nkit + prnt; // append prnt
    return (nlft);
}

static int
select_strs(struct KITDAT *kdp, int *strs, char *mems, int nmem)
{
    char sv[MAXSTR];
    int i, j, n, pn, pv, nstr, nuv;
    short nv[MAXSTR];

    nstr = kdp->nstr;
    pn = mems[nmem - 1] & 0xFF;
    n = 0;
    for (i = 0; i < nstr; i++) {
        pv = get_stro(pn, i);
        nuv = sv_histo(kdp, mems, nmem - 2, i, sv, nv);
        for (j = 0; j < nuv; j++) {
            if (sv[j] && (sv[j] != pv) && (nv[j] >= muthr)) {
                strs[n++] = i;
                break;
            }
        }
    }
    return (n);
}

static int
numsig(struct KITDAT *kdp, int node)
{
    int k, kidk, nsig, nkid;

    nsig = kdp->sig[node].n;
    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        nsig += numsig(kdp, kidk);
    }
    return (nsig);
}

static void
str_set_csv(struct KITDAT *kdp, char *sfn, int group)
{
    char ofn[MAXFNL], nodnam[40], sels[MAXKIT];
    int node, nmem, nkit, nstr, nsig, nsel, nnod;
    int i, j, k, n, year, val, ii[MAXSTR];
    FILE *ofp;
    static char *name = "sv_set";

    // output set STRs
    if (group) {
        sprintf(ofn, "%s%d.csv", name, group);
    } else {
        sprintf(ofn, "%s.csv", name);
    }
    printf(" %s\n", ofn);
    node = lookup_node(kdp, sfn);
    nmem = kdp->nmem[node];
    nsel = select_mems(kdp, sels, node);
    if (node == 0) {
        node = lookup_kits(kdp, sfn);
        if (node) {
            nsel = select_lfts(kdp, sels, node);
        }
    }
    // report subset summary
    nstr = select_strs(kdp, ii, sels, nsel);
    nsig = numsig(kdp, node);
    printf("subset: ");
    if (node) printf("%s ", sfn); 
    printf("node=%d nmem=%d nsel=%d ", node, nmem, nsel);
    printf("nstr=%d nsig=%d\n", nstr, nsig);
    // append descendants of node to set
    nnod = kdp->ngen + kdp->nsnp + kdp->nsbg;
    nkit = kdp->nkit;
    for (n = 0; n < nnod; n++) {
        if ((n != node) && isdescendant(kdp, n, node)) {
            sels[nsel++] = n + nkit;
        }
    }
    // output CSV file
    ofp = fopen(ofn, "w");
    fputs("kit", ofp);
    for (i = 0; i < nstr; i++) {
        fputs(",", ofp);
        fputs(kdp->labs[ii[i]], ofp);
    }
    fputs("\n", ofp);
    for (k = 0; k < nsel; k++) {
        j = sels[k] & 0xFF;
        year = kdp->kity[j];
        if ((j < nkit)  || (year == 0)) {
            sprintf(nodnam,"%s,", kdp->kits[j]);
        } else {
            sprintf(nodnam,"%s.%d,", kdp->kits[j], year);
        }
        fprintf(ofp,"%s", nodnam);
        for (i = 0; i < nstr; i++) {
            val = get_stro(j,ii[i]);
            fprintf(ofp, "%d", val);
            if (i < (nstr - 1)) {
                fputs(",", ofp);
            }
        }
        fputs("\n", ofp);
    }
    fclose(ofp);
}

static int
nodout(struct KITDAT *kdp, FILE *ofp, int n, int ndnt, int ps, int snpl)
{
    char *sn, name[256];
    int b, i, k, cv, mv, nc, nm, pv, presnp;
    int nkid, nlft, nsig, beg, len;
    static int in = 2;
    static int end = 78;

    if (n == 0) {
        sprintf(name, "%s - %d kits", kdp->gens[0], kdp->nkit);
        ndnt = in;
    } else {
        node_name(kdp, name, n);
        if (count) {
            nm = kdp->nmem[n];
            nc = strlen(name);
            sprintf(name + nc, "|%d", nm);
        }
    }
    fputs(name, ofp);
    len = ndnt + strlen(name);
    ndnt += in;
    beg = len - ps;
    if (lefts) {
        for (i = 0; i < kdp->nlft[n]; i++) {
            k = kdp->lfts[n][i];
            sprintf(name, " %s", kdp->kits[k]);
            if ((len + strlen(name) + ps) > end) {
                fprintf(ofp, "\n%*s", beg, "");
                len = beg;
                ps = 0;
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
            cv = kdp->sig[n].c[b];
            if (count) {
                sprintf(name, " %s=%d->%d|%d", sn, pv, mv, cv);
            } else {
                sprintf(name, " %s=%d->%d", sn, pv, mv);
            }
            if ((len + strlen(name) + ps) > end) {
                fprintf(ofp, "\n%*s", beg, "");
                len = beg;
                ps = 0;
            }
            fputs(name, ofp);
            len += strlen(name);
        }
    }
    nkid = kdp->nkid[n];
    nlft = kdp->nlft[n];
    for (i = 0; i < nkid; i++) {
        k = kdp->kids[n][i];
        if (kdp->nmem[k] == 0) continue;
        if (kdp->nmem[k] > kdp->nmem[n]) continue; // ???
        presnp = issnp(n) && isgen(k) && (nkid == 1) && !nlft
            && ((node_year(k) - node_year(n)) < PRESNPY);
        if (presnp) {
            fprintf(ofp, " ");
            ndnt--;
            ps = len - ndnt;
        } else {
            fprintf(ofp, "\n%*s", ndnt, ". ");
            ps = 0;
        }
	snpl = nodout(kdp, ofp, k, ndnt, ps, snpl);
    }
    return (snpl);
}

static void
report_tree(struct KITDAT *kdp, int group)
{
    char ofn[MAXFNL];
    FILE *ofp;
    static char *name = "tree";

    if (group) {
        sprintf(ofn, "%s%d.txt", name, group);
    } else {
        sprintf(ofn, "%s.txt", name);
    }
    printf("output %s", ofn);
    ofp = fopen(ofn, "w");
    nodout(kdp, ofp, 0, 0, 0, 0);
    fputs("\n", ofp);
    fclose(ofp);
}

static void
report_input(char *ifn, struct KITDAT *kdp)
{
    int nmem, nsig;

    nmem = kdp->nmem[0];
    nsig = numsig(kdp, 0);
    printf("tremer: %s\n", ifn);
    printf("      : nkit=%d nstr=%d ", kdp->nkit, kdp->nstr);
    printf("nsnp=%d ngen=%d nsbg=%d ", kdp->nsnp, kdp->ngen, kdp->nsbg);
    printf("mran=%d nmem=%d nsig=%d\n", kdp->mran, nmem, nsig);
}

static void
report_murate(struct KITDAT *kdp)
{
    float yrs, mut;
    int i, aby, fsa, nmut, ngen, nsnp, nkit, nsig;

    // report STR mutation rate
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nkit = kdp->nkit;
    nsig = nmut = 0;
    for (i = 0; i < (ngen + nsnp); i++) {
        nsig += kdp->sig[i].n * kdp->nmem[i];
        nmut += kdp->sig[i].n;
    }
    aby = 1960;         // average member birth year
    fsa = kdp->snpy[0]; // first SNP age
    yrs = aby - fsa;    // range of years
    mut = nkit * yrs / nsig;   
    if (nmut) {
        printf("STR mutation: nmut=%d rate= %d/%.0f/%d = 1/%.0f years\n",
            nmut, nsig, yrs, nkit, mut);
    } else {
        //printf("STR mutation: none\n");
    }
}

static void
report_ancestral(struct KITDAT *kdp, char *sfn, int group)
{
    char nnm[40], ofn[40];
    int i, k, kidk, node, nkid, nkit, nstr;
    FILE *ofp;
    static char *nm = "ancestral";

    // output set STRs
    if (group) {
        sprintf(ofn, "%s%d.txt", nm, group);
    } else {
        sprintf(ofn, "%s.txt", nm);
    }
    node = lookup_node(kdp, sfn);
    //node = kdp->mran;
    nkit = kdp->nkit;
    nstr = kdp->nstr;
    nkid = kdp->nkid[0];
    ofp = fopen(ofn, "w");
    if (ofp == NULL) return;
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        node_name(kdp, nnm, kidk);
        fputs(nnm, ofp);
        for (i = 0; i < nstr; i++) {
            fprintf(ofp, " %d", kdp->kstrs[nkit + kidk][i]);
        }
        fputs("\n", ofp);
    }
    fclose(ofp);
}

static int
group_number(char *nfn)
{
   while (*nfn && !isdigit(*nfn)) nfn++;
   return (atoi(nfn));
}

int
main(int argc, char **argv)
{
    char ifn[MAXFNL], sfn[MAXFNL];
    struct KITDAT *kdp;
    static int group = 0;

    if (argc < 2) {
        printf("usage: tremer [options] filename\n");
        printf("options:\n");
        printf("  -c    report number of mutation occurences\n");
        printf("  -d    enable debug print\n");
        printf("  -s N  report STRs for subset N (mut_thr=1) \n");
        printf("  -S N  report STRs for subset N (mut_thr=2) \n");
	exit(1);
    }
    while (argc > 1) {
        if (argv[1][0] == '-') {
            if (argv[1][1] == 'c') {
                count++;
            } else if (argv[1][1] == 'd') {
                debug++;
            } else if (argv[1][1] == 's') {
                muthr = 1; // select STRs with one mutation
                if (argc > 2) {
                    strncpy(sfn, argv[2], MAXFNL);
                    argv++;
                    argc--;
                }
            } else if (argv[1][1] == 'S') {
                muthr = 2; // select STRs with two mutations
                if (argc > 2) {
                    strncpy(sfn, argv[2], MAXFNL);
                    argv++;
                    argc--;
                }
            }
        } else {
            strncpy(ifn, argv[1], MAXFNL);
            group = group_number(ifn);
            // input KIT data
            kdp = get_sapp(ifn, group);
            report_input(ifn, kdp);
            // report STR
            str_adjust(kdp);
            report_murate(kdp);
            // report SNP GEN
            report_tree(kdp, group);
            report_ancestral(kdp, sfn, group);
            kit_gensnp_csv(kdp, group);
            gen_snp_csv(kdp, group);
            name_branches(kdp);
            // report STR
            str_all_csv(kdp, group);
            str_set_csv(kdp, sfn, group);
            plot_tree(kdp, ifn, group);
            // clean up
            free_kitdat(kdp);
        }
        argv++;
        argc--;
    }
}
