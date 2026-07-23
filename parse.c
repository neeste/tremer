#include "tremer.h"

int peek_next_char(FILE* file) {
    int c;
    while ((c = fgetc(file)) != EOF) {
        if (c == '*') {
            while ((c = fgetc(file)) != EOF && c != '\n') { }
        }
        else if (c > 32) { ungetc(c, file); return c; }
    }
    return EOF;
}

int get_next_token(FILE* file, char* token) {
    int c; int i = 0;
    while ((c = fgetc(file)) != EOF) {
        if (c == '*') { 
            while ((c = fgetc(file)) != EOF && c != '\n') { }
        }
        else if (c > 32) break;
    }
    if (c == EOF) return 0;
    token[i++] = c;
    while ((c = fgetc(file)) != EOF && c > 32) {
        if (c == '*') { ungetc(c, file); break; }
        if (i < MAX_STRING_LEN - 1) token[i++] = c;
    }
    token[i] = '\0';
    return 1;
}

int get_line_rest(FILE* file, char* line, int max_size) {
    int c; int i = 0;
    while ((c = fgetc(file)) != EOF) {
        if (c == '*') { 
            while ((c = fgetc(file)) != EOF && c != '\n') { }
            break; 
        }
        if (c == '\n' || c == '\r') break;
        if (i < max_size - 1) line[i++] = c;
    }
    line[i] = '\0';
    return (i > 0);
}

int find_kit(const char* id) {
    for (int i = 0; i < kit_count; i++) {
        if (strcmp(kits[i].id, id) == 0) return i;
    }
    return -1;
}

int get_or_create_kit(const char* id) {
    int idx = find_kit(id);
    if (idx != -1) return idx;
    if (kit_count >= MAX_KITS) {
        fprintf(stderr, "Error: Maximum kit capacity reached. Increase MAX_KITS and recompile.\n");
        return -1; // Gracefully exit
    }
    if (kit_count >= MAX_KITS - 1) return 0;
    strcpy(kits[kit_count].id, id);
    for (int m = 0; m < MAX_MARKERS; m++) kits[kit_count].str_values[m] = STR_MISSING;
    return kit_count++;
}

int is_kit_id(const char* token) {
    if (strcmp(token, "N") == 0) return 0;
    if (strncmp(token, "00", 2) == 0) return 0;
    for (int i = 0; token[i] != '\0'; i++) if ((token[i] >= 'A' && token[i] <= 'Z') || (token[i] >= 'a' && token[i] <= 'z')) return 1;
    if (atoi(token) > 200) return 1;
    return 0;
}

void strip_snp_affixes(char* token) {
    if (!token) return;
    int len = strlen(token);
    while (len > 0 && (token[len-1] == '+' || token[len-1] == '-' || token[len-1] == '?' || token[len-1] == '*')) {
        token[len-1] = '\0';
        len--;
    }
    int start = 0;
    while (token[start] == '~' || token[start] == '^' || token[start] == '*') {
        start++;
    }
    if (start > 0 && start <= len) {
        memmove(token, token + start, len - start + 1);
    }
}

void clean_snp_name(char* snp_name) {
    if (strlen(snp_name) >= 3 && snp_name[1] == '-') {
        if ((snp_name[0] >= 'A' && snp_name[0] <= 'Z') || (snp_name[0] >= 'a' && snp_name[0] <= 'z')) {
            if (root_haplogroup_letter == '\0') {
                root_haplogroup_letter = snp_name[0];
            }
            memmove(snp_name, snp_name + 2, strlen(snp_name) - 1);
        }
    }
}

void parse_markers(FILE* file) {
    char line[4096]; 
    int m_idx = 0;
    
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            char delim = ' ';
            if (strchr(line, ',')) delim = ',';
            else if (strchr(line, '\t')) delim = '\t';
            
            char* ptr = line; char* token;
            
            while (*ptr != '\0') {
                if (delim == ',' || delim == '\t') {
                    token = ptr;
                    char* next_delim = strchr(ptr, delim);
                    if (next_delim) {
                        *next_delim = '\0';
                        ptr = next_delim + 1;
                    } else {
                        ptr += strlen(ptr);
                    }
                } else {
                    while (*ptr == ' ') ptr++;
                    if (*ptr == '\0') break;
                    token = ptr;
                    while (*ptr != ' ' && *ptr != '\0') ptr++;
                    if (*ptr != '\0') { *ptr = '\0'; ptr++; }
                }
                
                char clean_tok[MAX_STRING_LEN] = {0};
                int c_idx = 0;
                for (int i = 0; token[i] != '\0' && c_idx < MAX_STRING_LEN - 1; i++) {
                    if (token[i] != '"' && token[i] != ' ' && token[i] != '\r' && token[i] != '\n') {
                        clean_tok[c_idx++] = token[i];
                    }
                }
                clean_tok[c_idx] = '\0';
                
                if (strlen(clean_tok) == 0) continue;
                
                char lower_tok[MAX_STRING_LEN];
                for (int i = 0; clean_tok[i] != '\0'; i++) {
                    lower_tok[i] = (clean_tok[i] >= 'A' && clean_tok[i] <= 'Z') ? clean_tok[i] + 32 : clean_tok[i];
                }
                lower_tok[strlen(clean_tok)] = '\0';
                
                if (strstr(lower_tok, "kit") || strstr(lower_tok, "name") || 
                    strstr(lower_tok, "group") || strstr(lower_tok, "ancestor") || 
                    strstr(lower_tok, "country") || strstr(lower_tok, "snp") ||
                    strcmp(lower_tok, "number") == 0) {
                    continue; 
                }
                
                if (m_idx < MAX_MARKERS) {
                    strcpy(marker_names[m_idx++], clean_tok);
                }
            }
        }
    }
    if (marker_count < m_idx) marker_count = m_idx;
}

void parse_modal(FILE* file) {
    char line[4096]; 
    int m_idx = 0;
    
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            char delim = ' ';
            if (strchr(line, ',')) delim = ',';
            else if (strchr(line, '\t')) delim = '\t';
            
            char* ptr = line; char* token;
            
            while (*ptr != '\0') {
                if (delim == ',' || delim == '\t') {
                    token = ptr;
                    char* next_delim = strchr(ptr, delim);
                    if (next_delim) {
                        *next_delim = '\0';
                        ptr = next_delim + 1;
                    } else {
                        ptr += strlen(ptr);
                    }
                } else {
                    while (*ptr == ' ') ptr++;
                    if (*ptr == '\0') break;
                    token = ptr;
                    while (*ptr != ' ' && *ptr != '\0') ptr++;
                    if (*ptr != '\0') { *ptr = '\0'; ptr++; }
                }
                
                char clean_tok[MAX_STRING_LEN] = {0};
                int c_idx = 0;
                for (int i = 0; token[i] != '\0' && c_idx < MAX_STRING_LEN - 1; i++) {
                    if (token[i] != '"' && token[i] != ' ' && token[i] != '\r' && token[i] != '\n') {
                        clean_tok[c_idx++] = token[i];
                    }
                }
                clean_tok[c_idx] = '\0';
                
                if (strlen(clean_tok) == 0) {
                    if (m_idx < MAX_MARKERS) modal_values[m_idx++] = STR_MISSING;
                    continue;
                }
                
                char lower_tok[MAX_STRING_LEN];
                for (int i = 0; clean_tok[i] != '\0'; i++) {
                    lower_tok[i] = (clean_tok[i] >= 'A' && clean_tok[i] <= 'Z') ? clean_tok[i] + 32 : clean_tok[i];
                }
                lower_tok[strlen(clean_tok)] = '\0';
                
                if (strstr(lower_tok, "modal") || strstr(lower_tok, "mode") || strstr(lower_tok, "group")) {
                    continue; 
                }
                
                if (m_idx < MAX_MARKERS) {
                    if (clean_tok[0] == 'N' || clean_tok[0] == 'n' || clean_tok[0] == 'X' || clean_tok[0] == 'x') {
                        modal_values[m_idx] = STR_MISSING;
                    } else {
                        int val = atoi(clean_tok);
                        modal_values[m_idx] = (val == 0) ? STR_MISSING : val;
                    }
                    m_idx++;
                }
            }
        }
    }
    if (marker_count < m_idx) marker_count = m_idx;
}

void parse_snptree(FILE* file) {
    char line[1024]; 
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            char tokens[50][MAX_STRING_LEN];
            int token_count = 0;
            
            char line_copy[1024];
            strcpy(line_copy, line);
            
            char* tok = strtok(line_copy, " \t\r\n");
            while (tok != NULL && token_count < 50) {
                strcpy(tokens[token_count++], tok);
                tok = strtok(NULL, " \t\r\n");
            }
            
            if (token_count < 2) continue; 
            
            int is_tremer = 1;
            for (int i = 0; tokens[1][i] != '\0'; i++) {
                if (tokens[1][i] < '0' || tokens[1][i] > '9') {
                    if (i == 0 && tokens[1][i] == '-') continue; 
                    is_tremer = 0; 
                    break;
                }
            }
            
            if (is_tremer) {
                char name[MAX_STRING_LEN] = {0};
                char parent_name[MAX_STRING_LEN] = {0};
                int date = atoi(tokens[1]);
                
                strcpy(name, tokens[0]);
                clean_snp_name(name);
                
                if (token_count >= 3) {
                    strcpy(parent_name, tokens[2]);
                }
                
                if (strlen(parent_name) == 1) {
                    if (root_haplogroup_letter == '\0') root_haplogroup_letter = parent_name[0];
                    parent_name[0] = '\0';
                } else if (strlen(parent_name) > 1) {
                    clean_snp_name(parent_name);
                }
                
                int is_dup = 0;
                for (int i = 0; i < snp_hierarchy_count; i++) {
                    if (strcmp(snp_hierarchy[i].name, name) == 0) { is_dup = 1; break; }
                }
                
                if (!is_dup && snp_hierarchy_count < MAX_SNPS) {
                    strcpy(snp_hierarchy[snp_hierarchy_count].parent_name_str, parent_name);
                    strcpy(snp_hierarchy[snp_hierarchy_count].name, name);
                    snp_hierarchy[snp_hierarchy_count].mrca_date = date;
                    snp_hierarchy[snp_hierarchy_count].parent = NULL;
                    snp_hierarchy_count++;
                }
                
            } else {
                for (int i = 1; i < token_count; i++) {
                    char parent_name[MAX_STRING_LEN];
                    strcpy(parent_name, tokens[i-1]);
                    clean_snp_name(parent_name);
                    
                    char child_name[MAX_STRING_LEN];
                    strcpy(child_name, tokens[i]);
                    clean_snp_name(child_name);
                    
                    if (strlen(parent_name) == 1 && i == 1) {
                        if (root_haplogroup_letter == '\0') root_haplogroup_letter = parent_name[0];
                        parent_name[0] = '\0';
                    }
                    
                    int is_dup = 0;
                    for (int j = 0; j < snp_hierarchy_count; j++) {
                        if (strcmp(snp_hierarchy[j].name, child_name) == 0) {
                            if (strlen(snp_hierarchy[j].parent_name_str) == 0) {
                                strcpy(snp_hierarchy[j].parent_name_str, parent_name);
                            }
                            is_dup = 1;
                            break;
                        }
                    }
                    
                    if (!is_dup && snp_hierarchy_count < MAX_SNPS) {
                        strcpy(snp_hierarchy[snp_hierarchy_count].parent_name_str, parent_name);
                        strcpy(snp_hierarchy[snp_hierarchy_count].name, child_name);
                        snp_hierarchy[snp_hierarchy_count].mrca_date = 0; 
                        snp_hierarchy[snp_hierarchy_count].parent = NULL;
                        snp_hierarchy_count++;
                    }
                }
            }
        }
    }
}

void parse_gendata(FILE* file) {
    char token[MAX_STRING_LEN]; int current_gen_idx = -1;
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_next_token(file, token)) {
            if (current_gen_idx == -1) {
                if (gen_hierarchy_count >= MAX_GEN_GROUPS - 1) break;
                char ancestor[MAX_STRING_LEN]; int date = 0;
                if (strchr(token, '.') != NULL) {
                    sscanf(token, "%63[^.].%d", ancestor, &date);
                } else {
                    strcpy(ancestor, token);
                    date = 0;
                }
                current_gen_idx = gen_hierarchy_count++;
                strcpy(gen_hierarchy[current_gen_idx].ancestor_name, ancestor);
                gen_hierarchy[current_gen_idx].date = date;
            } else {
                char* ptr = token; if (*ptr == '(') ptr++;
                if (*ptr == ')') current_gen_idx = -1;
                else if (*ptr != '\0') {
                    char kit_id[MAX_STRING_LEN]; int i = 0;
                    while (*ptr != '\0' && *ptr != ')') { if (i < MAX_STRING_LEN - 1) kit_id[i++] = *ptr; ptr++; }
                    kit_id[i] = '\0';
                    
                    char status = '+'; 
                    if (i > 0 && (kit_id[i-1] == '+' || kit_id[i-1] == '-' || kit_id[i-1] == '?')) {
                        status = kit_id[i-1];
                        kit_id[i-1] = '\0';
                    }
                    
                    if (strlen(kit_id) > 0) {
                        int index = find_kit(kit_id);
                        
                        // THE FIX: Add ALL kits to the genealogical group array, 
                        // including explicitly excluded ('-') kits, so downstream logic 
                        // can enforce the negative constraints.
                        if (index != -1) {
                            if (gen_hierarchy[current_gen_idx].kit_count < MAX_KITS) {
                                gen_hierarchy[current_gen_idx].kit_indices[gen_hierarchy[current_gen_idx].kit_count] = index;
                                gen_hierarchy[current_gen_idx].kit_statuses[gen_hierarchy[current_gen_idx].kit_count] = status;
                                gen_hierarchy[current_gen_idx].kit_count++;
                            }
                        }
                    }
                    if (*ptr == ')') current_gen_idx = -1;
                }
            }
        }
    }
}

void parse_strdata(FILE* file) {
    char line[2048];
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            
            char prefix[12]; strncpy(prefix, line, 11); prefix[11] = '\0';
            for(int i = 0; prefix[i]; i++) if(prefix[i] >= 'A' && prefix[i] <= 'Z') prefix[i] += 32;
            if (strncmp(prefix, "kit number", 10) == 0 || strncmp(prefix, "kit,", 4) == 0 || strncmp(prefix, "kit\t", 4) == 0 || strncmp(prefix, "kit ", 4) == 0) {
                continue; 
            }
            
            char delim = ' ';
            if (strchr(line, ',')) delim = ',';
            else if (strchr(line, '\t')) delim = '\t';
            
            char* ptr = line; char* token;
            int current_kit_idx = -1; int m_idx = 0; int token_count = 0;
            
            while (*ptr != '\0') {
                if (delim == ',' || delim == '\t') {
                    token = ptr;
                    char* next_delim = strchr(ptr, delim);
                    if (next_delim) {
                        *next_delim = '\0';
                        ptr = next_delim + 1;
                    } else {
                        ptr += strlen(ptr);
                    }
                } else {
                    while (*ptr == ' ') ptr++;
                    if (*ptr == '\0') break;
                    token = ptr;
                    while (*ptr != ' ' && *ptr != '\0') ptr++;
                    if (*ptr != '\0') { *ptr = '\0'; ptr++; }
                }
                
                token_count++;
                
                char clean_tok[MAX_STRING_LEN] = {0};
                int c_idx = 0;
                for (int i = 0; token[i] != '\0' && c_idx < MAX_STRING_LEN - 1; i++) {
                    if (token[i] != '"' && token[i] != ' ' && token[i] != '\r' && token[i] != '\n') {
                        clean_tok[c_idx++] = token[i];
                    }
                }
                clean_tok[c_idx] = '\0';
                
                if (token_count == 1) {
                    if (strlen(clean_tok) > 0 && is_kit_id(clean_tok)) {
                        current_kit_idx = get_or_create_kit(clean_tok);
                    } else { break; } 
                } 
                else {
                    if (token_count == 2) {
                        int is_pure_str = 1;
                        if (strlen(clean_tok) == 0) is_pure_str = 0; 
                        for(int i=0; clean_tok[i]!='\0'; i++) {
                            if (!((clean_tok[i] >= '0' && clean_tok[i] <= '9') || clean_tok[i] == '-' || 
                                  clean_tok[i] == 'N' || clean_tok[i] == 'n' || 
                                  clean_tok[i] == 'X' || clean_tok[i] == 'x')) {
                                is_pure_str = 0; break;
                            }
                        }
                        if (!is_pure_str) continue; 
                    }
                    
                    if (current_kit_idx != -1 && m_idx < MAX_MARKERS) {
                        if (strlen(clean_tok) == 0) {
                            kits[current_kit_idx].str_values[m_idx++] = STR_MISSING;
                            continue;
                        }
                        
                        char* sub_token = strtok(clean_tok, "-");
                        while (sub_token != NULL && m_idx < MAX_MARKERS) {
                            if (sub_token[0] == 'N' || sub_token[0] == 'n' || sub_token[0] == 'X' || sub_token[0] == 'x') {
                                kits[current_kit_idx].str_values[m_idx++] = STR_MISSING;
                            } else {
                                int val = atoi(sub_token);
                                kits[current_kit_idx].str_values[m_idx++] = (val == 0) ? STR_MISSING : val;     
                            }
                            sub_token = strtok(NULL, "-");
                        }
                    }
                }
            }
            if (current_kit_idx != -1 && m_idx > marker_count) {
                marker_count = m_idx;
            }
        }
    }
}

void parse_snpdata(FILE* file) {
    char line[4096]; 
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            char prefix[12]; strncpy(prefix, line, 11); prefix[11] = '\0';
            for(int i = 0; prefix[i]; i++) if(prefix[i] >= 'A' && prefix[i] <= 'Z') prefix[i] += 32;
            if (strncmp(prefix, "kit number", 10) == 0 || strncmp(prefix, "kit,", 4) == 0 || strncmp(prefix, "kit\t", 4) == 0 || strncmp(prefix, "kit ", 4) == 0) {
                continue; 
            }
            
            char delim = ' ';
            if (strchr(line, ',')) delim = ',';
            else if (strchr(line, '\t')) delim = '\t';
            
            char* ptr = line; char* token;
            int current_kit_idx = -1; int token_count = 0;
            
            while (*ptr != '\0') {
                if (delim == ',' || delim == '\t') {
                    token = ptr;
                    char* next_delim = strchr(ptr, delim);
                    if (next_delim) {
                        *next_delim = '\0';
                        ptr = next_delim + 1;
                    } else {
                        ptr += strlen(ptr);
                    }
                } else {
                    while (*ptr == ' ') ptr++;
                    if (*ptr == '\0') break;
                    token = ptr;
                    while (*ptr != ' ' && *ptr != '\0') ptr++;
                    if (*ptr != '\0') { *ptr = '\0'; ptr++; }
                }
                
                token_count++;
                
                char clean_tok[MAX_STRING_LEN] = {0};
                int c_idx = 0;
                for (int i = 0; token[i] != '\0' && c_idx < MAX_STRING_LEN - 1; i++) {
                    if (token[i] != '"' && token[i] != '(' && token[i] != ')' && token[i] != ' ' && token[i] != '\r' && token[i] != '\n') {
                        clean_tok[c_idx++] = token[i];
                    }
                }
                clean_tok[c_idx] = '\0';
                
                if (strlen(clean_tok) == 0) continue;
                
                if (token_count == 1) {
                    if (is_kit_id(clean_tok)) {
                        current_kit_idx = get_or_create_kit(clean_tok);
                    } else { break; } 
                } 
                else if (current_kit_idx != -1) {
                    int len = strlen(clean_tok);
                    char last_char = clean_tok[len-1];
                    
                    if (last_char == '+' || last_char == '-' || last_char == '?' || last_char == '*') {
                        SnpStatus status = SNP_POSITIVE;
                        if (last_char == '+') status = SNP_POSITIVE;
                        else if (last_char == '-' || last_char == '*') status = SNP_NEGATIVE;
                        else if (last_char == '?') status = SNP_UNTESTED;
                        
                        strip_snp_affixes(clean_tok);
                        clean_snp_name(clean_tok);
                        
                        if (strlen(clean_tok) > 0) {
                            SnpNode* new_snp = malloc(sizeof(SnpNode));
                            strcpy(new_snp->name, clean_tok);
                            new_snp->status = status;
                            
                            new_snp->next = kits[current_kit_idx].snps;
                            kits[current_kit_idx].snps = new_snp;
                        }
                    }
                }
            }
        }
    }
}

void parse_kitdata(FILE* file) {
    char token[MAX_STRING_LEN];
    while (peek_next_char(file) != EOF && peek_next_char(file) != '/') {
        if (get_next_token(file, token)) {
            char* dot_ptr = strchr(token, '.');
            if (dot_ptr != NULL) {
                *dot_ptr = '\0'; // Split kit ID from year
                int year = atoi(dot_ptr + 1);
                // 1. Strip absolutely everything but letters and numbers from the token
                char clean_token[MAX_STRING_LEN];
                int c = 0;
                for (int i = 0; token[i] != '\0'; i++) {
                    if ((token[i] >= 'A' && token[i] <= 'Z') || 
                        (token[i] >= 'a' && token[i] <= 'z') || 
                        (token[i] >= '0' && token[i] <= '9')) {
                        clean_token[c++] = token[i];
                    }
                }
                clean_token[c] = '\0';
                // 2. Forgiving Match against EXISTING kits
                for (int i = 0; i < kit_count; i++) {
                    
                    // Clean the existing kit ID the exact same way
                    char clean_existing[MAX_STRING_LEN];
                    int ec = 0;
                    for (int j = 0; kits[i].id[j] != '\0'; j++) {
                        if ((kits[i].id[j] >= 'A' && kits[i].id[j] <= 'Z') || 
                            (kits[i].id[j] >= 'a' && kits[i].id[j] <= 'z') || 
                            (kits[i].id[j] >= '0' && kits[i].id[j] <= '9')) {
                            clean_existing[ec++] = kits[i].id[j];
                        }
                    }
                    clean_existing[ec] = '\0';
                    // 3. Substring match: If one string safely exists inside the other, it's a match!
                    if (strstr(clean_existing, clean_token) != NULL || strstr(clean_token, clean_existing) != NULL) {
                        kits[i].birth_year = year;
                        break; 
                    }
                }
            }
        }
    }
}

void parse_groups(FILE* file) {
    char line[MAX_NODE_NAME_LEN];
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        if (get_line_rest(file, line, sizeof(line))) {
            char* open_paren = strchr(line, '(');
            if (open_paren != NULL) {
                char label[MAX_STRING_LEN] = {0};
                int len = open_paren - line;
                if (len >= MAX_STRING_LEN) len = MAX_STRING_LEN - 1;
                strncpy(label, line, len);
                label[len] = '\0';
                
                for (int i = len - 1; i >= 0 && (label[i] == ' ' || label[i] == '\t'); i--) {
                    label[i] = '\0';
                }
                
                char* label_start = label;
                while (*label_start == ' ' || *label_start == '\t') label_start++;
                
                if (strlen(label_start) > 0) {
                    char* ptr = open_paren + 1;
                    while (*ptr != '\0' && *ptr != ')') {
                        while (*ptr == ' ' || *ptr == '\t') ptr++;
                        if (*ptr == '\0' || *ptr == ')') break;
                        
                        char kit_id[MAX_STRING_LEN];
                        int i = 0;
                        while (*ptr != ' ' && *ptr != '\t' && *ptr != ')' && *ptr != '\0') {
                            if (i < MAX_STRING_LEN - 1) kit_id[i++] = *ptr;
                            ptr++;
                        }
                        kit_id[i] = '\0';
                        
                        int kit_len = strlen(kit_id);
                        if (kit_len > 0 && (kit_id[kit_len-1] == '+' || kit_id[kit_len-1] == '-' || kit_id[kit_len-1] == '?')) {
                            kit_id[kit_len-1] = '\0';
                        }
                        
                        if (strlen(kit_id) > 0) {
                            int index = find_kit(kit_id);
                            if (index != -1) {
                                strncpy(kits[index].group_label, label_start, MAX_STRING_LEN - 1);
                                kits[index].group_label[MAX_STRING_LEN - 1] = '\0';
                            }
                        }
                    }
                }
            }
        }
    }
}

// ==========================================
// NEW: SAFE BLOCK IGNORER
// ==========================================
void parse_ignore(FILE* file) {
    char line[4096];
    // Consume and discard lines until the next block header (starting with '/') or End of File
    while (peek_next_char(file) != '/' && peek_next_char(file) != EOF) {
        get_line_rest(file, line, sizeof(line));
    }
}
