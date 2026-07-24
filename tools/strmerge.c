/* strmerge.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "../src/fallback_labels.h"

#define MAXFNL 256

// Helper to normalize marker names (remove DYS/DYF and make lowercase)
void normalize_marker_name(const char* raw, char* norm) {
    int i = 0, j = 0;
    while (raw[i]) {
        norm[j++] = tolower(raw[i]);
        i++;
    }
    norm[j] = '\0';
    if (strncmp(norm, "dys", 3) == 0) {
        memmove(norm, norm + 3, strlen(norm) - 2);
    } else if (strncmp(norm, "dyf", 3) == 0) {
        memmove(norm, norm + 3, strlen(norm) - 2);
    }
}

int csv_read(char *fn, char **kitnam, int **kitstr, int *max_idx_out) {
    char header[1048576], data[1048576], s[8192];
    char *b, *e;
    FILE *fp;

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

    if (!fgets(header, sizeof(header), fp)) { fclose(fp); return 0; }
    if (!fgets(data, sizeof(data), fp)) { fclose(fp); return 0; }
    fclose(fp);

    int col_mapping[4096]; // maps CSV column index to base fallback_labels index (-1 if skip)
    for(int i = 0; i < 4096; i++) col_mapping[i] = -1;
    
    char *hptr = header;
    int col_idx = 0;

    // Build lookup for fallback labels
    char norm_fallback[838][64];
    for (int i = 0; i < n_lab; i++) {
        normalize_marker_name(str_lab[i], norm_fallback[i]);
    }

    while (*hptr != '\0' && *hptr != '\r' && *hptr != '\n') {
        char tok[256] = {0};
        int i = 0;
        
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

        char norm_tok[256];
        normalize_marker_name(tok, norm_tok);

        if (strstr(norm_tok, "kit") || strstr(norm_tok, "name") || strstr(norm_tok, "haplogroup") || strstr(norm_tok, "snp")) {
            col_mapping[col_idx] = -1;
        } else {
            // Find base index. If multi-copy, match the 'a' variant or the base name
            int mapped = -1;
            char search_a[256];
            sprintf(search_a, "%sa", norm_tok);
            
            for (int k = 0; k < n_lab; k++) {
                if (strcmp(norm_tok, norm_fallback[k]) == 0 || strcmp(search_a, norm_fallback[k]) == 0) {
                    mapped = k;
                    break;
                }
            }
            // For edge cases like "y-gata-h4" -> "gatah4"
            if (mapped == -1 && strstr(norm_tok, "gata-h4")) {
                for(int k=0; k<n_lab; k++) { if(strcmp(norm_fallback[k], "gatah4") == 0) { mapped = k; break; } }
            }
            if (mapped == -1 && strstr(norm_tok, "gata-a10")) {
                for(int k=0; k<n_lab; k++) { if(strcmp(norm_fallback[k], "yga10") == 0) { mapped = k; break; } }
            }
            if (mapped == -1 && strstr(norm_tok, "ggaat-1b07")) {
                for(int k=0; k<n_lab; k++) { if(strcmp(norm_fallback[k], "y1b07") == 0) { mapped = k; break; } }
            }
            if (mapped == -1 && strstr(norm_tok, "ycaii")) {
                for(int k=0; k<n_lab; k++) { if(strcmp(norm_fallback[k], "ycaiia") == 0) { mapped = k; break; } }
            }
            if (mapped == -1 && strcmp(norm_tok, "389ii") == 0) {
                for(int k=0; k<n_lab; k++) { if(strcmp(norm_fallback[k], "389ii-i") == 0) { mapped = k; break; } }
            }
            
            if (col_idx < 4096) col_mapping[col_idx] = mapped;
        }
        col_idx++;
        if (*hptr == ',') hptr++;
    }

    int temp_str[838] = {0}; // 0 implies 'N'
    char *dptr = data;
    col_idx = 0;
    int max_mapped_idx = -1;

    while (*dptr != '\0' && *dptr != '\r' && *dptr != '\n') {
        char tok[256] = {0};
        int i = 0;
        
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

        int mapped = (col_idx < 4096) ? col_mapping[col_idx] : -1;
        if (mapped != -1) {
            char ctok[256] = {0};
            int ci = 0;
            for(int k = 0; tok[k]; k++) {
                if (tok[k] != ' ') ctok[ci++] = tok[k];
            }
            ctok[ci] = '\0';

            if (strlen(ctok) == 0 || strcmp(ctok, "-") == 0) {
                // leave as 0
            } 
            else if (strchr(ctok, '-')) {
                int multi_idx = mapped;
                char *sub = strtok(ctok, "-");
                while (sub && multi_idx < n_lab) { // Assign to mapped, mapped+1, etc.
                    temp_str[multi_idx] = atoi(sub);
                    if (multi_idx > max_mapped_idx) max_mapped_idx = multi_idx;
                    multi_idx++;
                    sub = strtok(NULL, "-");
                }
            } 
            else {
                temp_str[mapped] = atoi(ctok);
                if (mapped > max_mapped_idx) max_mapped_idx = mapped;
            }
        }

        col_idx++;
        if (*dptr == ',') dptr++;
    }

    *kitstr = calloc(max_mapped_idx + 1, sizeof(int));
    for(int i = 0; i <= max_mapped_idx; i++) {
        (*kitstr)[i] = temp_str[i];
    }
    *max_idx_out = max_mapped_idx;

    return max_mapped_idx + 1; // Number of markers to output
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
            strcnt[csvcnt] = csv_read(av[1], kitnam + csvcnt, kitstr + csvcnt, &n);
            csvcnt++;
        }
        ac--;
        av++;
    }
    
    printf("Found %d CSV files.\n", csvcnt);
    
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
                fprintf(ofp, " N");
            }
        }
        fprintf(ofp, "\n");
        free(kitnam[i]);
        free(kitstr[i]);
    }
    fclose(ofp);
    return 0;
}
