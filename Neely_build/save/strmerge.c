/* strmerge.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXFNL 256

int csv_read(char *fn, char **kitnam, int **kitstr) {
    char header[8192], data[8192], s[8192];
    char *b, *e;
    FILE *fp;

    // 1. Robust Kit Name Extraction (with fallback if _YDNA_ is missing)
    strcpy(s, fn);
    e = strstr(s, "_YDNA_");
    if (e) {
        *e = '\0';
    } else {
        e = strrchr(s, '.');
        if (e) *e = '\0';
    }
    for (b = s + strlen(s); b > s; b--) {
        if (b[-1] == '/' || b[-1] == '\\') break;
    }
    *kitnam = strdup(b);

    fp = fopen(fn, "r");
    if (!fp) return 0;

    // Read Header and Data lines
    if (!fgets(header, sizeof(header), fp)) { fclose(fp); return 0; }
    if (!fgets(data, sizeof(data), fp)) { fclose(fp); return 0; }
    fclose(fp);

    // 2. Parse Header to identify metadata columns to skip
    int skip_col[1024] = {0};
    char *hptr = header;
    int col_idx = 0;

    while (*hptr != '\0' && *hptr != '\r' && *hptr != '\n') {
        char tok[256] = {0};
        int i = 0;
        
        // CSV extraction respecting quotes
        if (*hptr == '"') {
            hptr++;
            while (*hptr != '"' && *hptr != '\0' && *hptr != '\r' && *hptr != '\n') {
                if (i < 255) tok[i++] = *hptr;
                hptr++;
            }
            if (*hptr == '"') hptr++;
        } else {
            while (*hptr != ',' && *hptr != '\0' && *hptr != '\r' && *hptr != '\n') {
                if (i < 255) tok[i++] = *hptr;
                hptr++;
            }
        }
        tok[i] = '\0';

        // Convert to lowercase for heuristic matching
        for (int k = 0; tok[k]; k++) {
            if (tok[k] >= 'A' && tok[k] <= 'Z') tok[k] += 32;
        }

        if (strstr(tok, "kit") || strstr(tok, "name") || strstr(tok, "haplogroup") || strstr(tok, "snp")) {
            skip_col[col_idx] = 1;
        }
        col_idx++;

        if (*hptr == ',') hptr++;
    }

    // 3. Parse Data columns
    int temp_str[2048];
    int nint = 0;
    char *dptr = data;
    col_idx = 0;

    while (*dptr != '\0' && *dptr != '\r' && *dptr != '\n') {
        char tok[256] = {0};
        int i = 0;
        
        // CSV extraction respecting quotes
        if (*dptr == '"') {
            dptr++;
            while (*dptr != '"' && *dptr != '\0' && *dptr != '\r' && *dptr != '\n') {
                if (i < 255) tok[i++] = *dptr;
                dptr++;
            }
            if (*dptr == '"') dptr++;
        } else {
            while (*dptr != ',' && *dptr != '\0' && *dptr != '\r' && *dptr != '\n') {
                if (i < 255) tok[i++] = *dptr;
                dptr++;
            }
        }
        tok[i] = '\0';

        // Process token only if it's not a skipped metadata column
        if (!skip_col[col_idx]) {
            // Clean token: remove internal spaces
            char ctok[256] = {0};
            int ci = 0;
            for(int k = 0; tok[k]; k++) {
                if (tok[k] != ' ') ctok[ci++] = tok[k];
            }
            ctok[ci] = '\0';

            // Handle Missing Values
            if (strlen(ctok) == 0 || strcmp(ctok, "-") == 0) {
                temp_str[nint++] = 0; 
            } 
            // Handle Multi-copy hyphenated markers
            else if (strchr(ctok, '-')) {
                char *sub = strtok(ctok, "-");
                while (sub) {
                    temp_str[nint++] = atoi(sub);
                    sub = strtok(NULL, "-");
                }
            } 
            // Handle standard integer
            else {
                temp_str[nint++] = atoi(ctok);
            }
        }

        col_idx++;
        if (*dptr == ',') dptr++;
    }

    // Allocate memory for the exact number of integers found
    *kitstr = calloc(nint, sizeof(int));
    for(int i = 0; i < nint; i++) {
        (*kitstr)[i] = temp_str[i];
    }

    return nint;
}

int main(int ac, char **av) {
    char *kitnam[512];
    int csvcnt, i, j, n, strcnt[512], *kitstr[512];
    FILE *ofp;
    static char ofn[MAXFNL] = "strmerge.txt";

    csvcnt = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'o') {
                if (ac > 2) {
                    strncpy(ofn, av[2], MAXFNL);
                    av++;
                    ac--;
                }
            }
        } else {
            strcnt[csvcnt] = csv_read(av[1], kitnam + csvcnt, kitstr + csvcnt);
            csvcnt++;
        }
        ac--;
        av++;
    }
    
    printf("Found %d CSV files.\n", csvcnt);
    
    // Write download file
    ofp = fopen(ofn, "w");
    if (!ofp) {
        fprintf(stderr, "Error opening output file %s\n", ofn);
        return 1;
    }
    
    for (i = 0; i < csvcnt; i++) {
        fprintf(ofp, "* Kit %s has %d STRs.\n", kitnam[i], strcnt[i]);
        fprintf(ofp, "%s", kitnam[i]);
        
        for (j = 0; j < strcnt[i]; j++) {
            n = kitstr[i][j];
            if (n) {
                fprintf(ofp, " %d", n);
            } else {
                fprintf(ofp, " N"); // Missing values output as N
            }
        }
        fprintf(ofp, "\n");
        free(kitnam[i]);
        free(kitstr[i]);
    }
    fclose(ofp);
    return 0;
}
