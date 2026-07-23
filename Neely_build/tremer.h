// tremer.h

#define BIGY    500
#define NOPT    100
#define PRESNPY  25
#define MAXSV    80
#define MAXSBG   40
#define MAXGRP   10
#define MAXKIT  160
#define MAXSNP  160
#define MAXGEN  160
#define MAXHDR  160
#define MAXSIG  160
#define MAXFNL  256
#define MAXSTR  900*2
#define MAXDSC  900*2
#define MAXLIN 8192

// Removed free_null definition: now defined in utils.h

// --- START: Implementation of Suggestion 8 (Node Indexing) ---
// Base indices for different node types within the combined node array (kstrs, kits, etc.).
// Node indices (n) are 0 to (ngen + nsnp + nsbg - 1).
// Kit indices (k) are 0 to (nkit - 1).

#define GEN_BASE_INDEX     0
#define SNP_BASE_INDEX     kdp->ngen
#define SBG_BASE_INDEX     (kdp->ngen + kdp->nsnp)
#define NODE_STR_OFFSET    kdp->nkit // The first non-kit index in the kstrs/kits arrays is kdp->nkit

// Macros for node STR value access using the explicit offset
#define get_strn(n,i)   (kdp->kstrs[NODE_STR_OFFSET + n][i])
#define set_strn(n,i,v) (kdp->kstrs[NODE_STR_OFFSET + n][i]=v)
// --- END: Implementation of Suggestion 8 ---

#define get_stro(k,i)  ((k<kdp->nkit)?kdp->ostrs[k][i]:kdp->kstrs[k][i])
#define get_strk(k,i)   (kdp->kstrs[k][i])
#define set_strk(k,i,v) (kdp->kstrs[k][i]=v)
#define isgen(n)        (n<kdp->ngen)
#define issnp(n)        ((n>=kdp->ngen)&&(n<(kdp->ngen+kdp->nsnp)))
#define node_year(n)    (kdp->kity[(kdp->nkit)+n])

struct SIGN {
    int s, p, m;
};
struct SIGDAT {
    int n, s[MAXSIG], p[MAXSIG], m[MAXSIG], c[MAXSIG], d[MAXSIG];
};
struct KITDAT {
    int nkit, nstr, ngen, nsnp, nsbg, nlab, ngrp, nxmb, tnmv, mran;
    short nsvs[MAXKIT], nmvs[MAXKIT], term[MAXKIT], sso[MAXSTR];
    short kity[MAXKIT], geny[MAXGEN], snpy[MAXSNP];
    char *kits[MAXKIT], *kstrs[MAXKIT], *ostrs[MAXKIT];
    char *gens[MAXGEN], *snps[MAXSNP], *sbgs[MAXSBG], *labs[MAXSTR];
    char ksnp[MAXKIT][MAXSNP], kgen[MAXKIT][MAXGEN], gsnp[MAXGEN][MAXSNP];
    char top[MAXHDR], pars[MAXSNP], *grps[MAXGRP], *snpp[MAXSNP];
    int nmem[MAXGEN], nkid[MAXGEN], nlft[MAXGEN];
    char parent[MAXGEN], kids[MAXGEN][MAXGEN];
    char xmbn[MAXGEN], xmbk[MAXGEN];
    char mems[MAXGEN][MAXKIT], lfts[MAXGEN][MAXKIT], kgrp[MAXKIT];
    struct SIGDAT sig[MAXGEN];
};

void plot_tree(struct KITDAT *, char *, int);
void find_comut_kits(struct KITDAT *kdp, int node, int i, int v,
    char *lfts_out, int *mlft_out, char *kids_out, int *mkid_out);
