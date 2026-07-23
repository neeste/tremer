#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <time.h>

#define MAX_SNPS 10000
#define MAX_STR 128

typedef struct {
    char name[MAX_STR];
    int age;
    char parent[MAX_STR];
    char filename[MAX_STR]; 
} TreeEntry;

TreeEntry* tree_list;
int tree_count = 0;

int is_numeric(const char* str) {
    if (!str || *str == '\0') return 0;
    for (int i = 0; str[i] != '\0'; i++) {
        if (!isdigit(str[i]) && str[i] != '-') return 0;
    }
    return 1;
}

void extract_basename(const char* path, char* output) {
    const char* base = strrchr(path, '/');
    base = base ? base + 1 : path;
    strncpy(output, base, MAX_STR - 1);
    output[MAX_STR - 1] = '\0';
}

void clean_snp_name(char* snp_name) {
    if (strlen(snp_name) >= 3 && snp_name[1] == '-') {
        if ((snp_name[0] >= 'A' && snp_name[0] <= 'Z') || (snp_name[0] >= 'a' && snp_name[0] <= 'z')) {
            memmove(snp_name, snp_name + 2, strlen(snp_name) - 1);
        }
    }
}

// --- ALTERNATE MODE: Validation & Recommendations ---
void run_alternate_mode(const char* out_filename, char** input_files, int input_count, int force_update) {
    printf("\n[Alternate Mode] '%s' is newer than the input files.\n", out_filename);
    printf("Validating input files against the reference output...\n\n");

    TreeEntry* ref_tree = malloc(sizeof(TreeEntry) * MAX_SNPS);
    int ref_count = 0;
    FILE* f_out = fopen(out_filename, "r");
    if (f_out) {
        char line[512];
        while (fgets(line, sizeof(line), f_out)) {
            char token1[MAX_STR] = {0}, token2[MAX_STR] = {0}, token3[MAX_STR] = {0};
            int items = sscanf(line, "%127s %127s %127s", token1, token2, token3);
            if (items >= 1) {
                strcpy(ref_tree[ref_count].name, token1);
                ref_tree[ref_count].age = 0;
                strcpy(ref_tree[ref_count].parent, "");
                
                if (items == 3) {
                    ref_tree[ref_count].age = atoi(token2);
                    strcpy(ref_tree[ref_count].parent, token3);
                } else if (items == 2) {
                    if (is_numeric(token2)) ref_tree[ref_count].age = atoi(token2);
                    else strcpy(ref_tree[ref_count].parent, token2);
                }
                ref_count++;
            }
        }
        fclose(f_out);
    } else {
        printf("Error: Could not read reference file %s\n", out_filename);
        free(ref_tree);
        return;
    }

    int total_recommendations = 0;
    int files_with_snptree = 0;

    for (int i = 0; i < input_count; i++) {
        FILE* file = fopen(input_files[i], "r");
        if (!file) continue;

        char temp_filename[MAX_STR];
        snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", input_files[i]);
        FILE* temp_file = NULL;
        if (force_update) temp_file = fopen(temp_filename, "w");

        int in_snptree = 0;
        int has_snptree = 0; 
        int file_changed = 0;
        char line[1024];
        int recommendations = 0;

        while (fgets(line, sizeof(line), file)) {
            char first_token[MAX_STR];
            if (sscanf(line, "%127s", first_token) != 1) {
                if (temp_file) fputs(line, temp_file);
                continue;
            }

            if (first_token[0] == '/') {
                in_snptree = (strcmp(first_token, "/SNPTREE") == 0);
                if (in_snptree) has_snptree = 1;
                if (temp_file) fputs(line, temp_file);
                continue;
            }

            if (in_snptree && first_token[0] != '*') {
                char tokens[50][MAX_STR];
                int t_count = 0;
                char line_copy[1024];
                strcpy(line_copy, line);
                
                char* tok = strtok(line_copy, " \t\r\n");
                while (tok && t_count < 50) { strcpy(tokens[t_count++], tok); tok = strtok(NULL, " \t\r\n"); }

                if (t_count < 2) {
                    if (temp_file) fputs(line, temp_file);
                    continue;
                }

                int is_tremer = 1;
                for (int k = 0; tokens[1][k] != '\0'; k++) {
                    if (!isdigit(tokens[1][k]) && tokens[1][k] != '-') { is_tremer = 0; break; }
                }

                if (is_tremer) {
                    char name[MAX_STR] = {0}; strcpy(name, tokens[0]); clean_snp_name(name);
                    int age = atoi(tokens[1]);
                    char parent[MAX_STR] = {0}; 
                    if (t_count >= 3) { strcpy(parent, tokens[2]); clean_snp_name(parent); }

                    int found_ref = 0;
                    for (int r = 0; r < ref_count; r++) {
                        if (strcmp(ref_tree[r].name, name) == 0) {
                            found_ref = 1;
                            if (ref_tree[r].age != age || strcmp(ref_tree[r].parent, parent) != 0) {
                                if (recommendations == 0) printf("--- Updates %s in %s ---\n", force_update ? "applied" : "recommended", input_files[i]);
                                printf("  Current : %-12s %-5d %s\n", name, age, parent);
                                printf("  Expected: %-12s %-5d %s\n\n", ref_tree[r].name, ref_tree[r].age, ref_tree[r].parent);
                                recommendations++;
                                file_changed = 1;

                                if (temp_file) {
                                    if (strlen(ref_tree[r].parent) > 0) fprintf(temp_file, "%-15s %-5d %s\n", ref_tree[r].name, ref_tree[r].age, ref_tree[r].parent);
                                    else fprintf(temp_file, "%-15s %-5d\n", ref_tree[r].name, ref_tree[r].age);
                                }
                            } else {
                                if (temp_file) fputs(line, temp_file); // No changes needed
                            }
                            break;
                        }
                    }
                    if (!found_ref && temp_file) fputs(line, temp_file);

                } else {
                    // SAPP formatting detected: evaluate, and if -f is passed, always standardise to TREMER format
                    char new_sapp_lines[4096] = {0};
                    int sapp_mismatch = 0;

                    for (int j = 1; j < t_count; j++) {
                        char child_name[MAX_STR], parent_name[MAX_STR];
                        strcpy(child_name, tokens[j]); clean_snp_name(child_name);
                        strcpy(parent_name, tokens[j-1]); clean_snp_name(parent_name);

                        int found_ref = 0;
                        for (int r = 0; r < ref_count; r++) {
                            if (strcmp(ref_tree[r].name, child_name) == 0) {
                                found_ref = 1;
                                if (strcmp(ref_tree[r].parent, parent_name) != 0 || ref_tree[r].age != 0) {
                                    if (recommendations == 0) printf("--- Updates %s in %s ---\n", force_update ? "applied" : "recommended", input_files[i]);
                                    printf("  Current (SAPP)  : %s -> %s\n", parent_name, child_name);
                                    printf("  Expected (TREMER): %-12s %-5d %s\n\n", ref_tree[r].name, ref_tree[r].age, ref_tree[r].parent);
                                    recommendations++;
                                    sapp_mismatch = 1;
                                }

                                char temp_buf[256];
                                if (strlen(ref_tree[r].parent) > 0) snprintf(temp_buf, sizeof(temp_buf), "%-15s %-5d %s\n", ref_tree[r].name, ref_tree[r].age, ref_tree[r].parent);
                                else snprintf(temp_buf, sizeof(temp_buf), "%-15s %-5d\n", ref_tree[r].name, ref_tree[r].age);
                                strcat(new_sapp_lines, temp_buf);
                                break;
                            }
                        }
                        if (!found_ref) {
                            char temp_buf[256];
                            snprintf(temp_buf, sizeof(temp_buf), "%-15s %-5d %s\n", child_name, 0, parent_name);
                            strcat(new_sapp_lines, temp_buf);
                        }
                    }

                    if (sapp_mismatch || force_update) file_changed = 1;
                    if (temp_file) {
                        if (force_update) fputs(new_sapp_lines, temp_file);
                        else fputs(line, temp_file);
                    }
                }
            } else {
                if (temp_file) fputs(line, temp_file); // Outside SNPTREE block or is a comment
            }
        }
        fclose(file);
        
        if (temp_file) {
            fclose(temp_file);
            if (file_changed && force_update) {
                remove(input_files[i]);
                rename(temp_filename, input_files[i]);
            } else {
                remove(temp_filename);
            }
        }

        if (has_snptree) files_with_snptree++;
        total_recommendations += recommendations;
    }

    if (files_with_snptree > 0 && total_recommendations == 0) {
        printf("--- No updates needed in any input file ---\n\n");
    } else if (total_recommendations > 0 && force_update) {
        printf("--- All expected updates have been successfully applied to the input files. ---\n\n");
    }

    free(ref_tree);
}

// --- NORMAL MODE: Standard Tree Aggregation ---
int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s [-o output_file] [-f] <input_file1> [input_file2...]\n", argv[0]);
        return 1;
    }

    char out_filename[MAX_STR] = "snptree.txt"; 
    char* input_files[100];
    int input_count = 0;
    int force_update = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            strncpy(out_filename, argv[++i], MAX_STR - 1);
        } else if (strcmp(argv[i], "-f") == 0) {
            force_update = 1;
        } else {
            input_files[input_count++] = argv[i];
        }
    }

    if (input_count == 0) {
        printf("Error: No input files specified.\n");
        return 1;
    }

    struct stat out_stat;
    int out_exists = (stat(out_filename, &out_stat) == 0);
    time_t out_mtime = out_exists ? out_stat.st_mtime : 0;
    
    time_t max_in_mtime = 0;
    for (int i = 0; i < input_count; i++) {
        struct stat in_stat;
        if (stat(input_files[i], &in_stat) == 0) {
            if (in_stat.st_mtime > max_in_mtime) {
                max_in_mtime = in_stat.st_mtime;
            }
        }
    }

    if (out_exists && out_mtime > max_in_mtime) {
        run_alternate_mode(out_filename, input_files, input_count, force_update);
        return 0;
    }

    tree_list = malloc(sizeof(TreeEntry) * MAX_SNPS);

    for (int i = 0; i < input_count; i++) {
        FILE* file = fopen(input_files[i], "r");
        if (!file) continue;

        char base_fn[MAX_STR];
        extract_basename(input_files[i], base_fn);

        int in_snptree = 0;
        char line[1024];

        while (fgets(line, sizeof(line), file)) {
            char first_token[MAX_STR];
            if (sscanf(line, "%127s", first_token) != 1) continue;

            if (first_token[0] == '/') {
                in_snptree = (strcmp(first_token, "/SNPTREE") == 0);
                continue;
            }

            if (in_snptree && first_token[0] != '*') {
                char tokens[50][MAX_STR];
                int t_count = 0;
                char line_copy[1024];
                strcpy(line_copy, line);
                
                char* tok = strtok(line_copy, " \t\r\n");
                while (tok && t_count < 50) { strcpy(tokens[t_count++], tok); tok = strtok(NULL, " \t\r\n"); }

                if (t_count < 2) continue;

                int is_tremer = 1;
                for (int k = 0; tokens[1][k] != '\0'; k++) {
                    if (!isdigit(tokens[1][k]) && tokens[1][k] != '-') { is_tremer = 0; break; }
                }

                if (is_tremer) {
                    char name[MAX_STR] = {0}; strcpy(name, tokens[0]); clean_snp_name(name);
                    int age = atoi(tokens[1]);
                    char parent[MAX_STR] = {0}; 
                    if (t_count >= 3) { strcpy(parent, tokens[2]); clean_snp_name(parent); }

                    int is_dup = 0;
                    for (int j = 0; j < tree_count; j++) {
                        if (strcmp(tree_list[j].name, name) == 0) { is_dup = 1; break; }
                    }
                    if (!is_dup && tree_count < MAX_SNPS) {
                        strcpy(tree_list[tree_count].name, name);
                        tree_list[tree_count].age = age;
                        strcpy(tree_list[tree_count].parent, parent);
                        strcpy(tree_list[tree_count].filename, base_fn);
                        tree_count++;
                    }
                } else {
                    for (int j = 1; j < t_count; j++) {
                        char child_name[MAX_STR], parent_name[MAX_STR];
                        strcpy(child_name, tokens[j]); clean_snp_name(child_name);
                        strcpy(parent_name, tokens[j-1]); clean_snp_name(parent_name);

                        int is_dup = 0;
                        for (int x = 0; x < tree_count; x++) {
                            if (strcmp(tree_list[x].name, child_name) == 0) {
                                if (strlen(tree_list[x].parent) == 0) strcpy(tree_list[x].parent, parent_name);
                                is_dup = 1; break;
                            }
                        }
                        if (!is_dup && tree_count < MAX_SNPS) {
                            strcpy(tree_list[tree_count].name, child_name);
                            tree_list[tree_count].age = 0;
                            strcpy(tree_list[tree_count].parent, parent_name);
                            strcpy(tree_list[tree_count].filename, base_fn);
                            tree_count++;
                        }
                    }
                }
            }
        }
        fclose(file);
    }

    FILE* out = fopen(out_filename, "w");
    if (!out) {
        printf("Error: Could not open output file %s\n", out_filename);
        free(tree_list);
        return 1;
    }
    
    int max_name_len = 0;
    int max_age_len = 0;
    for (int i = 0; i < tree_count; i++) {
        int n_len = strlen(tree_list[i].name);
        if (n_len > max_name_len) max_name_len = n_len;
        
        if (tree_list[i].age > 0) {
            char age_str[32];
            sprintf(age_str, "%d", tree_list[i].age);
            int a_len = strlen(age_str);
            if (a_len > max_age_len) max_age_len = a_len;
        }
    }

    for (int i = 0; i < tree_count; i++) {
        if (strlen(tree_list[i].parent) > 0 && tree_list[i].age > 0) {
            fprintf(out, "%-*s %-*d %s\n", max_name_len, tree_list[i].name, max_age_len, tree_list[i].age, tree_list[i].parent);
        } else if (tree_list[i].age > 0) {
            fprintf(out, "%-*s %d\n", max_name_len, tree_list[i].name, tree_list[i].age);
        } else if (strlen(tree_list[i].parent) > 0) {
            fprintf(out, "%-*s %-*s %s\n", max_name_len, tree_list[i].name, max_age_len, "0", tree_list[i].parent);
        } else {
            fprintf(out, "%s\n", tree_list[i].name);
        }
    }

    fclose(out);
    printf("Aggregated SNP tree successfully written to %s\n", out_filename);
    
    free(tree_list);
    return 0;
}
