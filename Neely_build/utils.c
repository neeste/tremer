// utils.c - Utility functions for tremer

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "utils.h"

// Note: MAXSTR is defined in utils.h for function signature compatibility.

int
getval(char *line, char *str)
{
    char *s;
    int nsv = 0;

    nsv = 0;
    s = line;
    // create list of STR values
    while (*s && (nsv < MAXSTR)) {
        while(isspace(*s) || (*s=='-')) s++;
        if (*s == 'N') {
            str[nsv++] = 0;
        } else if ((*s >= '1') && (*s <= '9')) {
            str[nsv++] = (char) atoi(s);
        }
        while(*s && !(isspace(*s) || (*s=='-'))) s++;
    }
    return (nsv);
}

int
unique(char *v, char *uv, int n, int nz)
{
    int i, j, nuv;

    nuv = 0;
    for (i = 0; i < n; i++) {
        if (nz && (v[i] == 0)) continue; // only consider non-zero values
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
count_val(char *v, int n, int u)
{
    int i, cnt;

    cnt = 0;
    for (i = 0; i < n; i++) {
        if (v[i] == u) cnt++;
    }
    return (cnt);
}

