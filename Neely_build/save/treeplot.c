// treeplot.c - export tree plot as SVG

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "tremer.h"

#define XPAGE 11.0
#define YPAGE  8.0

static int wg = 200;
static int hg = 200;
static int wn = 180;
static int hn =  80;
static int wk = 180;
static int hk = 160;
static int wv =   0;
static int hv =   0;
static int mg =  20;
static int ch =  18;

void
svg_open(FILE *ofp, char *title, char *creator, char *date, int wv, int hv)
{
    static char *s = "square";
     
    fprintf(ofp, " <svg width=\"%d\" height=\"%d\"", wv, hv);
    fprintf(ofp, " viewBox=\"%d %d %d %d\"\n", 0, 0, wv, hv);
    fprintf(ofp, " xmlns=\"http://www.w3.org/2000/svg\">\n");
    fprintf(ofp, " <g style=\"stroke-linecap:%s;", s);
    fprintf(ofp, " stroke-linejoin:%s\">\n", s);
    fprintf(ofp, "<title> %s </title>\n", title);
    fprintf(ofp, "<desc>\n");
    fprintf(ofp, "Creator: %s\n", creator);
    fprintf(ofp, "Date: %s\n", date);
    fprintf(ofp, "</desc>\n");
    fprintf(ofp, "<rect fill=\"#fff\" stroke=\"#000\" ");
    fprintf(ofp, "x=\"%d\" y=\"%d\" ", 0, 0);
    fprintf(ofp, "width=\"%d\" height=\"%d\"/>\n", wv, hv);
}

void
svg_path(FILE *ofp, int *x, int *y, int n, char *c, int w)
{
    char *s;
    int i;

    fprintf(ofp, "<polyline points=\"");
    for (i = 0; i < n; i++) {
        s = (i < (n - 1)) ? " " : "";
        fprintf(ofp, "%d,%d%s", x[i], y[i], s);
    }
    fprintf(ofp, "\" style=\"fill:none;stroke:%s;", c);
    fprintf(ofp, "stroke-width:%d\" />\n", w);
}

void
svg_text(FILE *ofp, int xpos, int ypos, char *s, int a, int h, int c)
{
    static char *align[] = {"middle", "end", "begin"};
    static char *size[] = {"small", "medium", "large"};
    static char *color[] = {"black", "red"};

    fprintf(ofp, "<text font-size=\"%s\" ", size[h]);
    fprintf(ofp, "fill=\"%s\" ", color[c]);
    fprintf(ofp, "x=\"%d\" y=\"%d\" ", xpos, ypos);
    fprintf(ofp, "text-anchor=\"%s\">", align[a]);
    fprintf(ofp, "%s</text>\n", s);
}

void
svg_node(FILE *ofp, int xpos, int ypos, char s[8][20], int n, int c)
{
    int i, x, y, x1, y1;
    static char *fc[] = {"lightblue", "lightcyan"};

    x = mg + xpos;
    y = mg + ypos;
    fprintf(ofp, "<rect x=\"%d\" y=\"%d\" ", x, y);
    fprintf(ofp, "width=\"%d\" height=\"%d\" ", wn, hn);
    fprintf(ofp, "fill=\"%s\" stroke-width=\"1\" ", fc[c]);
    fprintf(ofp, "stroke=\"lightgray\" />\n");
    x1 = x + wn / 2;
    y1 = y + (hn - (n - 1) * ch) / 2 + 5;
    for (i = 0; i < n; i++) {
        svg_text(ofp, x1, y1 + ch * i, s[i], 0, 1, 0);
    }
}

void
svg_left(FILE *ofp, int xpos, int ypos, char s[8][20], int *nv, int n, char *gn)
{
    char nvs[20];
    int i, x, y, x1, y1, x2;

    x = mg + xpos;
    y = mg + ypos;
    fprintf(ofp, "<rect x=\"%d\" y=\"%d\" ", x, y);
    fprintf(ofp, "width=\"%d\" height=\"%d\" ", wk, hk);
    fprintf(ofp, "fill=\"lightyellow\" stroke-width=\"1\" ");
    fprintf(ofp, "stroke=\"lightpink\" />\n");
    x1 = x + wk / 2;
    y1 = y + (hk - (n - 1) * ch) / 2 + 5;
    x2 = x + wk - ch;
    for (i = 0; i < n; i++) {
        sprintf(nvs, "%3d", nv[i]);
        svg_text(ofp, x1, y1 + ch * i, s[i], 1, 1, 0);
        svg_text(ofp, x2, y1 + ch * i, nvs, 1, 1, 0);
    }
    x = mg + xpos + ch / 2;
    y = mg + ypos + ch + hk;
    svg_text(ofp, x, y, gn, 2, 0, 1);
}

void
svg_close(FILE *ofp)
{
    static char *s[] = {"</g> <!-- EOF -->", "</svg>"};

    fprintf(ofp, "%s\n", s[0]);
    fprintf(ofp, "%s\n", s[1]);
}

void
position_node(struct KITDAT *kdp, char *fmt, int *xp, int *yp, int node)
{
    int k, kidk, prnt, xmax, ymax, nlft, nkid, ngen, nsnp, nsbg, nnod;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    prnt = kdp->parent[node];
    xmax = 0;
    ymax = 0;
    if (node == 0) {
        for (k = 0; k < nnod; k++) fmt[k] = 0;
        xp[0] = 0;
        yp[0] = 0;
        fmt[node] = 1;
    } else if (kdp->nmem[node] == 0) {
        xp[0] = 0;
        yp[0] = 0;
        fmt[node] = 0;
    } else {
        for (k = 0; k < nnod; k++) {
            if (fmt[k] == 0) continue;
            if (xmax < xp[k]) xmax = xp[k];
            if (ymax < yp[k]) ymax = yp[k];
        }
        fmt[node] = 1;
        yp[node] = hg + yp[prnt];
        xp[node] = (yp[node] > ymax) ? 0 : wg + xmax;
        nlft = 0;
        for (k = 1; k < nnod; k++) {
            if (fmt[k] == 0) continue;
            if (xp[k] != xp[node]) continue;
            if (yp[k] != (yp[node] - hg)) continue;
            nlft = kdp->nlft[k];
        }
        if (nlft) {
            xp[node] += wg;
            if (xmax < xp[node]) xmax = xp[node];
        }
    }
    nkid = kdp->nkid[node];
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        position_node(kdp, fmt, xp, yp, kidk);
    }
}

int
insert_kids(struct KITDAT *kdp, char *fmt, int *xp, int *yp, int node, int x0)
{
    int k, kidk, nkid, y0;

    xp[node] = x0;
    if (kdp->nlft[node]) x0 += wg;
    nkid = kdp->nkid[node];
    y0 = yp[node] + hg;
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        xp[kidk] = x0 + k * wg; 
        yp[kidk] = (fmt[kidk] == 2) ? y0 - hg : y0; 
        x0 = insert_kids(kdp, fmt, xp, yp, kidk, x0);
        if (fmt[kidk] == 2) yp[node] = yp[kidk];
    }
    return (x0);
}

void
adjust_ypos(struct KITDAT *kdp, char *fmt, int *xp, int *yp)
{
    int k, p, presnp, nkid, nlft, ngen, nsnp, nsbg, nnod;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    for (k = 1; k < nnod; k++) {
        p = kdp->parent[k];
        nkid = kdp->nkid[p];
        nlft = kdp->nlft[p];
        presnp = issnp(p) && isgen(k) && (nkid == 1) && !nlft
            && ((node_year(k) - node_year(p)) < PRESNPY);
        if (presnp) {
            //printf("%s > %s\n", kdp->snps[p - ngen], kdp->gens[k]);
            xp[k] = xp[p];
            yp[k] = yp[p] - hg;
            fmt[p] = 3;
            fmt[k] = 2;
        }
    }
    insert_kids(kdp, fmt, xp, yp, 0, 0);
}


int
xshift(struct KITDAT *kdp, char *fmt, int *xp, int *yp, int node)
{
    int k, kidk, xmin, xmax, xdif, xshf, nlft, nkid, ngen, nsnp;
    static int prvn = 0;

    if (node == prvn) {
        //printf("*** WARNING: xshift loop detected\n");
        return (0);
    }
    prvn = node;
    if (kdp->nmem[node] == kdp->nmem[0]) node = kdp->mran;
    nkid = kdp->nkid[node];
    nlft = kdp->nlft[node];
    if ((nlft == 0) && (nkid == 1)) {
        kidk = kdp->kids[node][0];
        if (kidk) {
            return (xshift(kdp, fmt, xp, yp, kidk));
        }
    }
    if (node == kdp->mran && nkid && !nlft) {
        kidk = kdp->kids[node][0];
        xmin = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
        kidk = kdp->kids[node][nkid - 1];
        xmax = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
        return ((xmin + xmax) / 2);
    }
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    if (node < (ngen + nsnp)) {
        if ((nlft != 0) && (nkid == 1)) {
            kidk = kdp->kids[node][0];
            xmax = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
            return (xmax / 2);
        }
        if ((nlft == 0) && (nkid > 1)) {
            kidk = kdp->kids[node][0];
            xmin = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
            kidk = kdp->kids[node][nkid - 1];
            xmax = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
            return ((xmax + xmin) / 2);
        }
    } else {
        if ((nlft == 0) && (nkid == 1)) {
            kidk = kdp->kids[node][0];
            return (xshift(kdp, fmt, xp, yp, kidk));
        }
        if ((nlft == 0) && (nkid > 1)) {
            kidk = kdp->kids[node][0];
            xmin = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
            kidk = kdp->kids[node][nkid - 1];
            xmax = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
            return ((xmax + xmin) / 2);
        }
    }
    xmax = xmin = 0;
    if ((nlft == 0) && (nkid == 2)) {
        kidk = kdp->kids[node][0];
        xmin = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
    }
    for (k = 0; k < nkid; k++) {
        kidk = kdp->kids[node][k];
        if (fmt[kidk] == 0) continue;
        if (kdp->nmem[kidk] == 0) continue;
        xdif = xp[kidk] + xshift(kdp, fmt, xp, yp, kidk) - xp[node];
        if (xmax < xdif) xmax = xdif;
    }
    xshf = (xmax + xmin) / 2;
    return (xshf);
}

void
plot_path(struct KITDAT *kdp, FILE *ofp, char *fmt, int *xp, int *yp)
{
    int x[MAXGEN], y[MAXGEN];
    int nkid, nlft, ngen, nsnp, nsbg, nnod;
    int i, m, n, hh, xs, xx, yy, w, isbg, p;
    static char *c = "lightblue";

    w = 4;
    hh = hg * 3 / 4;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    isbg = ngen + nsnp;
    nnod = ngen + nsnp + nsbg;
    for (i = 1; i < nnod; i++) {
        if (!fmt[i]) continue;
        if (!kdp->nmem[i]) continue;
        if (!kdp->nlft[i]) continue;
        // connect leftovers to node
        nkid = kdp->nkid[i];
        nlft = kdp->nlft[i];
        xx = xp[i] + wg / 2;
        yy = yp[i] + hg / 2;
        xs = xshift(kdp, fmt, xp, yp, i);
        n = 0;
        //if (nkid || (i < isbg)) {
        if ((i < isbg) || nkid || nlft) {
            x[n] = xx;
            y[n] = yy + hg;
            n++;
            x[n] = x[n - 1];
            y[n] = y[n - 1] - hh;
            n++;
            x[n] = x[n - 1] + xs;
            y[n] = y[n - 1];
            n++;
            x[n] = x[n - 1];
            y[n] = yy;
            n++;
        } else {
            x[n] = xx;
            y[n] = yy;
            n++;
        }
        // connect parents
        p = kdp->parent[i];
        while (p) {
            if (fmt[p]) {
                xx = xp[p]+ wg / 2;
                yy = yp[p]+ hg / 2;
                xs = xshift(kdp, fmt, xp, yp, p);
                //if (p == 6) xs = 175;
                x[n] = x[n - 1];
                y[n] = y[n - 1] - hh;
                n++;
                x[n] = xx + xs;
                y[n] = y[n - 1];
                n++;
                x[n] = x[n - 1];
                y[n] = yy;
                n++;
            }
            p = kdp->parent[p];
        }
        svg_path(ofp, x, y, n, c, w);
    }
    // connect ancestral nodes
    m = kdp->mran;
    x[0] = xp[m]+ wg / 2 + xshift(kdp, fmt, xp, yp, m);
    y[0] = yp[m]+ hg / 2;
    x[1] = xp[0]+ wg / 2 + xshift(kdp, fmt, xp, yp, 0);
    y[1] = yp[0]+ hg / 2;
    n = 2;
    svg_path(ofp, x, y, n, c, w);
}

void
adjust_xpos(struct KITDAT *kdp, char *fmt, int *xp, int *yp)
{
    int k, n, nh, xmax, ymax, ngen, nsnp, nsbg, nnod;

    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    nh = 0;
    for (k = 0; k < nnod; k++) { // check for completion
        if (!fmt[k]) nh++;
    }
    for (xmax = 0; xmax < (wg * (nnod - nh)); xmax += wg) {
        n = 0;
        for (k = 0; k < nnod; k++) { // check for completion
            if (xp[k] < xmax) n++;
        }
        if (n == (nnod - nh)) break;
        n = 0;
        for (k = 0; k < nnod; k++) { // check for gaps
            if (fmt[k] && (xp[k] >= xmax) && (xp[k] < (xmax + wg))) n++;
        }
        if (n) continue;
        for (k = 0; k < nnod; k++) { // fill gap
            if (xp[k] >= (xmax + wg)) xp[k] -= wg;
        }
        //xmax -= wg;
    }
    xmax = xp[0];
    ymax = yp[0];
    for (k = 1; k < nnod; k++) {
        if (fmt[k] == 0) continue;
        if (xmax < xp[k]) {
            xmax = xp[k];
        }
        if (ymax < yp[k]) {
            ymax = yp[k];
        }
    }
    wv = xmax + mg * 2 + wg;
    hv = ymax + mg * 2 + hg * 2;
}

void
plot_node(struct KITDAT *kdp, FILE *ofp, int *xp, int *yp, char *lab,
    int xshf, int node)
{
    char *s, *kits, *gn, labs[8][20];
    int nkid, nlab, nlft, nkit, ngen, nsnp, nsvs[8];
    int c, g, k, kitk, y, n, p, isbg, presnp, yplft;

    nkit = kdp->nkit;
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    isbg = ngen + nsnp;
    if (node  == 0) { // skip label
    } else if (node < isbg) {
        s = (node < ngen) ? kdp->gens[node] : kdp->snps[node - ngen];
        y = (node < ngen) ? kdp->geny[node] : kdp->snpy[node - ngen];
        if (y) {
            sprintf(lab, "%s.%d", s, y);
        } else {
            sprintf(lab, "%s", s);
        }
    } else {
        sprintf(lab, ".%d", node - isbg);
    }
    n = node;
    p = kdp->parent[n];
    nkid = kdp->nkid[p];
    nlft = kdp->nlft[p];
    presnp = issnp(p) && isgen(n) && (nkid == 1) && !nlft
            && ((node_year(n) - node_year(p)) < PRESNPY);
    if (node >= isbg) {
        nlab = 1;
        sprintf(labs[0], "STR%02d", node - isbg + 1);
        s = kdp->kits[node + nkit];
        y = kdp->kity[node + nkit];
        if (y > 0) {
            sprintf(labs[0], "%s.%d", s, y);
        } else {
            strcpy(labs[0], s);
        }
    } else if (presnp) {
        nlab = 2;
        s = kdp->snps[p - ngen];
        y = kdp->snpy[p - ngen];
        sprintf(labs[0], "%s.%d", s, y);
        strcpy(labs[1], lab);
        if (strcmp(labs[0], labs[1])==0) nlab = 1; // <--
    } else if (node == 0) {
        nlab = 2;
        strcpy(labs[0], lab);
        sprintf(labs[1], "%d kits", kdp->nkit);
    } else {
        nlab = 1;
        strcpy(labs[0], lab);
    }
    if ((n < isbg) || kdp->nkid[n] || kdp->nlft[n]) {
        c = (node >= ngen);
        svg_node(ofp, xp[node] + xshf, yp[node], labs, nlab, c);
        yplft = yp[node] + hg;
    } else {
        yplft = yp[node];
    }
    nlft = kdp->nlft[node];
    if (nlft > 0) {
        for (k = 0; k < nlft; k++) {
            kitk = kdp->lfts[node][k];
            kits = kdp->kits[kitk];
            sprintf(labs[k], "%s", kits);
            nsvs[k] = kdp->nsvs[kitk];
        }
        kitk = kdp->lfts[node][0];
        g = kdp->kgrp[kitk];
        gn = g ? kdp->grps[g] : "";
        svg_left(ofp, xp[node], yplft, labs, nsvs, nlft, gn);
    }
}

void
plot_strs(struct KITDAT *kdp, FILE *ofp, int *xp, int *yp, int xshf, 
    int node, int left)
{
    char *sn, s[8][20];
    int i, ns, pv, mv, cv, lc, x1, y1, yy, nsig, nkid;
    static int ch = 15;

    if (left) {
        if (kdp->nkid[node] == 0) return;
        if (kdp->nlft[node] == 0) return;
    }
    nsig = kdp->sig[node].n;
    if (nsig == 0) return;
    ns = (nsig < 8) ? nsig : 8;
    for (i = 0; i < ns; i++) {
        sn = kdp->labs[kdp->sig[node].s[i]];
        pv = kdp->sig[node].p[i];
        mv = kdp->sig[node].m[i];
        sprintf(s[i], " %s=%d->%d", sn, pv, mv);
    }
    if (ns < nsig) strcat(s[ns - 1], "...");
    x1 = xp[node] + 25;
    y1 = yp[node] + 15 - ch * (ns - 1);
    if (left) { y1 += hg; } else { x1 += xshf; }
    nkid = kdp->nkid[node];
    yy = y1;
    for (i = 0; i < ns; i++) {
        lc = (kdp->sig[node].d[i] < kdp->sig[node].c[i]);
        if ((lc != left) || (nkid == 0)) {
            svg_text(ofp, x1, yy, s[i], 2, 0, 0);
            yy += ch;
        }
    }
    for (i = 0; i < ns; i++) {
        cv = left ? kdp->sig[node].d[i] : kdp->sig[node].c[i];
        if (!cv) cv = kdp->sig[node].d[i]; // <-- DEBUG
        sprintf(s[i], "%3d", cv);
    }
    x1 = xp[node] + wg + - 3;
    y1 = yp[node] + 15 - ch * (ns - 1);
    if (left) { y1 += hg; } else { x1 += xshf; }
    yy = y1;
    for (i = 0; i < ns; i++) {
        lc = (kdp->sig[node].d[i] < kdp->sig[node].c[i]);
        if ((lc != left) || (nkid == 0)) {
            svg_text(ofp, x1, yy, s[i], 1, 0, 0);
            yy += ch;
        }
    }
}

void
newext(char *fn, char *ext)
{
    while (*fn) {
        if (*fn++ == '.') {
            strcpy(fn, ext);
            break;
        }
    }
}

void
plot_tree(struct KITDAT *kdp, char *ifn, int group)
{
    char ofn[80], date[80], name[80], title[80], fmt[MAXGEN], lab[80];
    int xp[MAXGEN], yp[MAXGEN];
    int k, xshf, ngen, nsnp, nsbg, nnod;
    time_t t;
    FILE *ofp;
    static char *creator = "tremer";

    strcpy(name, ifn);
    strcpy(ofn, ifn);
    newext(ofn, "svg");
    ngen = kdp->ngen;
    nsnp = kdp->nsnp;
    nsbg = kdp->nsbg;
    nnod = ngen + nsnp + nsbg;
    position_node(kdp, fmt, xp, yp, 0);
    adjust_ypos(kdp, fmt, xp, yp);
    adjust_xpos(kdp, fmt, xp, yp);
    sprintf(title, "%s", name);
    if (group == 12) {
        sprintf(lab, "Groups %d+%d", group / 10, group % 10);
    } else {
        sprintf(lab, "Group %d", group);
    }
    time(&t);
    strcpy(date, ctime(&t));
    date[strlen(date) - 1] = 0;
    ofp = fopen(ofn, "w");
    svg_open(ofp, title, creator, date, wv, hv);
    plot_path(kdp, ofp, fmt, xp, yp);
    for (k = 0; k < nnod; k++) {
        if (fmt[k]) {
            xshf = xshift(kdp, fmt, xp, yp, k);
            if (fmt[k] != 3) plot_node(kdp, ofp, xp, yp, lab, xshf, k);
            plot_strs(kdp, ofp, xp, yp, xshf, k, 0);
            plot_strs(kdp, ofp, xp, yp, xshf, k, 1);
        }
    }
    svg_close(ofp);
    fclose(ofp);
}

