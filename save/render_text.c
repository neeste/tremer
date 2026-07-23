#include "tremer.h"

int print_indent(FILE* out_file, int depth) {
    for (int i = 0; i < depth; i++) fprintf(out_file, "  "); 
    return depth * 2;
}

void print_item_wrapped(FILE* out_file, const char* text, int* current_col, int wrap_indent_col, int max_width) {
    if (text == NULL || strlen(text) == 0) return;
    int len = strlen(text);
    if (*current_col + len + 1 > max_width) {
        fprintf(out_file, "\n");
        for (int i = 0; i < wrap_indent_col; i++) fprintf(out_file, " ");
        *current_col = wrap_indent_col;
    }
    fprintf(out_file, "%s ", text);
    *current_col += len + 1;
}

// ==========================================
// ORPHANED KIT SWEEPER
// ==========================================
void text_find_assigned_kits(TreeNode* node, int assigned[]) {
    if (!node) return;
    for (int i = 0; i < node->kit_count; i++) {
        assigned[node->kit_indices[i]] = 1;
    }
    TreeNode* c = node->first_child;
    while (c) {
        text_find_assigned_kits(c, assigned);
        c = c->next_sibling;
    }
}

// ==========================================
// TEXT RENDERER
// ==========================================

void collect_all_html_kits(TreeNode* node, int* kit_array, int* count);

void collect_ancestral_str_kits(TreeNode* n, TreeNode* exclude_path, int* kit_array, int* count) {
    if (!n) return;
    for (int i = 0; i < n->kit_count; i++) {
        kit_array[(*count)++] = n->kit_indices[i];
    }
    TreeNode* c = n->first_child;
    while (c) {
        if (c != exclude_path && (c->type == NODE_STR || c->type == NODE_STR_BRANCH || c->type == NODE_DISTANCE)) {
            collect_all_html_kits(c, kit_array, count);
        }
        c = c->next_sibling;
    }
}

void collect_all_html_kits(TreeNode* node, int* kit_array, int* count) {
    if (!node) return;
    for (int i = 0; i < node->kit_count; i++) {
        kit_array[(*count)++] = node->kit_indices[i];
    }
    TreeNode* c = node->first_child;
    while (c) {
        collect_all_html_kits(c, kit_array, count);
        c = c->next_sibling;
    }
}

void generate_text_output_internal(FILE* out_file, TreeNode* node, int depth, int is_top) {
    if (node == NULL) return; 
    
    if (is_top) {
        TreeNode* ancestors[100];
        int ancestor_count = 0;
        TreeNode* curr = node;
        while (curr->parent != NULL && ancestor_count < 100) {
            ancestors[ancestor_count++] = curr->parent;
            curr = curr->parent;
        }
        
        for (int a = ancestor_count - 1; a >= 0; a--) {
            TreeNode* anc = ancestors[a];
            TreeNode* path_to_skip = (a == 0) ? node : ancestors[a-1];
            
            int anc_kits[MAX_KITS];
            int anc_count = 0;
            collect_ancestral_str_kits(anc, path_to_skip, anc_kits, &anc_count);
            
            if (anc_count > 0) {
                int current_col = print_indent(out_file, depth);
                char buffer[MAX_NODE_NAME_LEN + 32];
                if (anc->type == NODE_ROOT || strcmp(anc->name, "ROOT") == 0) {
                    sprintf(buffer, "Global Basal ");
                } else {
                    sprintf(buffer, "%s ", anc->name);
                }
                
                current_col += fprintf(out_file, "%s", buffer);
                int wrap_indent_col = (depth * 2) + 4; // Locked to structural depth
                for (int i = 0; i < anc_count; i++) {
                    print_item_wrapped(out_file, kits[anc_kits[i]].id, &current_col, wrap_indent_col, 80);
                }
                fprintf(out_file, "\n");
            }
        }
    }
    
    if (node->type == NODE_ROOT) { 
        if (node->kit_count > 0) {
            int current_col = print_indent(out_file, depth);
            current_col += fprintf(out_file, "ROOT ");
            int wrap_indent_col = (depth * 2) + 4; // Locked to structural depth
            for (int i = 0; i < node->kit_count; i++) {
                print_item_wrapped(out_file, kits[node->kit_indices[i]].id, &current_col, wrap_indent_col, 80);
            }
            fprintf(out_file, "\n");
        }
        TreeNode* child = node->first_child; 
        while (child != NULL) { 
            generate_text_output_internal(out_file, child, depth, 0); 
            child = child->next_sibling; 
        } 
        return; 
    }
    
    int current_col = print_indent(out_file, depth);
    int wrap_indent_col = (depth * 2) + 4;  // Locked to structural depth
    
    if (node->type == NODE_STR || node->type == NODE_STR_BRANCH) {
        char prefix_buf[MAX_NODE_NAME_LEN];
        // Allow STR nodes to display their dynamically calculated ages
        if (node->date > 0) {
            sprintf(prefix_buf, ". %s.%d", node->name, node->date);
        } else {
            sprintf(prefix_buf, ". %s", node->name);
        }
        print_item_wrapped(out_file, prefix_buf, &current_col, wrap_indent_col, 80);
        
        for (int i = 0; i < node->mutation_count; i++) {
            char mut_buf[64];
            sprintf(mut_buf, "%s=%d->%d", 
                    marker_names[node->mutations[i].marker_index], 
                    node->mutations[i].old_val, 
                    node->mutations[i].new_val);
            print_item_wrapped(out_file, mut_buf, &current_col, wrap_indent_col, 80);
        }
    } else if (node->type != NODE_ROOT && node->type != NODE_DISTANCE) {
        char name_buf[MAX_NODE_NAME_LEN + 32];
        if (node->date > 0) {
            sprintf(name_buf, ". %s.%d", node->name, node->date);
        } else {
            sprintf(name_buf, ". %s", node->name);
        }
        print_item_wrapped(out_file, name_buf, &current_col, wrap_indent_col, 80);
    }
    
    for (int i = 0; i < node->kit_count; i++) {
        print_item_wrapped(out_file, kits[node->kit_indices[i]].id, &current_col, wrap_indent_col, 80);
    }
    
    if (node->parent != NULL && node->type != NODE_STR && node->type != NODE_STR_BRANCH) {
        int subtree_kits[MAX_KITS]; int subtree_kit_count = 0; 
        collect_subtree_kits(node, subtree_kits, &subtree_kit_count);
        for (int m = 0; m < marker_count; m++) {
            if (node->local_modal[m] > 0 && 
                node->parent->local_modal[m] > 0 && 
                node->local_modal[m] != node->parent->local_modal[m]) {
                
                int occurrence_count = 0; 
                for (int k = 0; k < subtree_kit_count; k++) {
                    if (kits[subtree_kits[k]].str_values[m] == node->local_modal[m]) occurrence_count++;
                }

                if (occurrence_count >= 2) {
                    char buffer[MAX_NODE_NAME_LEN];
                    sprintf(buffer, "%s=%d->%d", marker_names[m], node->parent->local_modal[m], node->local_modal[m]);
                    print_item_wrapped(out_file, buffer, &current_col, wrap_indent_col, 80);
                }
            }
        }
    }
    
    fprintf(out_file, "\n");
    int child_depth = (node->type == NODE_DISTANCE) ? depth : depth + 1;
    TreeNode* child = node->first_child; 
    while (child != NULL) { 
        generate_text_output_internal(out_file, child, child_depth, 0); 
        child = child->next_sibling; 
    }
}

void generate_text_output(FILE* out_file, TreeNode* node, int depth) {
    generate_text_output_internal(out_file, node, depth, 1);
}

// ==========================================
// HTML TABLE OUTPUT STAGE
// ==========================================

void generate_html_table_output(const char* filename, TreeNode* root) {
    FILE* f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n<title>STR Marker Table</title>\n");
    fprintf(f, "<style>\n");
    fprintf(f, "body { font-family: sans-serif; margin: 20px; }\n");
    fprintf(f, "table { border-collapse: collapse; font-size: 12px; white-space: nowrap; }\n");
    fprintf(f, "th, td { border: 1px solid #ddd; text-align: center; }\n");
    fprintf(f, "td { padding: 4px; width: 20px; min-width: 20px; }\n");
    fprintf(f, "th.marker { background-color: #d9e2f3; position: sticky; top: 0; box-shadow: 0 2px 2px -1px rgba(0,0,0,0.1); height: 100px; padding: 4px; vertical-align: bottom; }\n");
    fprintf(f, "th.marker span { writing-mode: vertical-rl; transform: rotate(180deg); display: inline-block; }\n");
    fprintf(f, "td.kitid { text-align: left; padding: 4px 8px; font-weight: bold; background-color: #fafafa; position: sticky; left: 0; box-shadow: 2px 0 2px -1px rgba(0,0,0,0.1); }\n");
    fprintf(f, "th.kitid { background-color: #d9e2f3; position: sticky; left: 0; top: 0; z-index: 2; padding: 4px 8px; }\n");
    fprintf(f, "</style>\n</head>\n<body>\n");

    fprintf(f, "<h2>Group %d - STR Marker Values</h2>\n", global_group_num);
    fprintf(f, "<table>\n<thead>\n<tr><th class=\"kitid\">Kit ID</th>");

    int num_markers = (marker_count < 111) ? marker_count : 111;
    for (int m = 0; m < num_markers; m++) {
        fprintf(f, "<th class=\"marker\"><span>%s</span></th>", marker_names[m]);
    }
    fprintf(f, "</tr>\n</thead>\n<tbody>\n");

    int ordered_kits[MAX_KITS];
    int ordered_count = 0;
    
    collect_all_html_kits(root, ordered_kits, &ordered_count);

    for (int i = 0; i < ordered_count; i++) {
        int k_idx = ordered_kits[i];
        
        int is_basal = 0;
        
        for (int b = 0; b < root->kit_count; b++) {
            if (root->kit_indices[b] == k_idx) { is_basal = 1; break; }
        }
        
        if (is_basal && (root->type == NODE_ROOT || strcmp(root->name, "ROOT") == 0)) {
            fprintf(f, "<tr><td class=\"kitid\">%s <span style=\"color:#aaa; font-weight:normal; font-size:0.75em; font-style:italic;\">(ROOT)</span></td>", kits[k_idx].id);
        } else {
            fprintf(f, "<tr><td class=\"kitid\">%s</td>", kits[k_idx].id);
        }

        for (int m = 0; m < num_markers; m++) {
            int val = kits[k_idx].str_values[m];
            int modal = root->local_modal[m]; 
            if (val == STR_MISSING || val <= 0) {
                fprintf(f, "<td></td>");
            } else {
                if (modal != STR_MISSING && modal > 0) {
                    if (val > modal) {
                        fprintf(f, "<td style=\"background-color: #add8e6;\">%d</td>", val);
                    } else if (val < modal) {
                        fprintf(f, "<td style=\"background-color: #ffb6c1;\">%d</td>", val);
                    } else {
                        fprintf(f, "<td>%d</td>", val);
                    }
                } else {
                    fprintf(f, "<td>%d</td>", val);
                }
            }
        }
        fprintf(f, "</tr>\n");
    }

    fprintf(f, "<tr><td class=\"kitid\" style=\"background-color: #e8e8e8; border-top: 2px solid #333;\">MODAL</td>");
    for (int m = 0; m < num_markers; m++) {
        int g_mod = modal_values[m]; 
        if (g_mod == STR_MISSING || g_mod <= 0) {
            fprintf(f, "<td style=\"background-color: #e8e8e8; border-top: 2px solid #333; font-weight: bold;\">-</td>");
        } else {
            fprintf(f, "<td style=\"background-color: #e8e8e8; border-top: 2px solid #333; font-weight: bold;\">%d</td>", g_mod);
        }
    }
    fprintf(f, "</tr>\n");

    fprintf(f, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose(f);
}

void generate_html_filtered_table_output(const char* filename, TreeNode* root) {
    FILE* f = fopen(filename, "w");
    if (!f) return;

    // 1. We ONLY collect kits from the target node and its direct descendants.
    int ordered_kits[MAX_KITS];
    int ordered_count = 0;
    collect_all_html_kits(root, ordered_kits, &ordered_count);

    // 2. Establish Baseline Modal
    int baseline_modal[MAX_MARKERS];
    for (int m = 0; m < marker_count; m++) {
        if (root->parent != NULL) {
            baseline_modal[m] = root->parent->local_modal[m];
        } else {
            baseline_modal[m] = root->local_modal[m];
        }
    }

    // 3. Find which markers have mutations to display
    int kept_markers[MAX_MARKERS];
    int kept_count = 0;

    for (int m = 0; m < marker_count; m++) {
        int diff_count = 0;
        int modal = baseline_modal[m]; 
        if (modal == STR_MISSING || modal <= 0) continue;
        
        for (int i = 0; i < ordered_count; i++) {
            int val = kits[ordered_kits[i]].str_values[m];
            if (val > 0 && val != modal) diff_count++;
        }
        if (diff_count >= 2) kept_markers[kept_count++] = m;
    }

    if (kept_count < 4) {
        kept_count = 0; 
        for (int m = 0; m < marker_count; m++) {
            int diff_count = 0;
            int modal = baseline_modal[m]; 
            if (modal == STR_MISSING || modal <= 0) continue;
            
            for (int i = 0; i < ordered_count; i++) {
                int val = kits[ordered_kits[i]].str_values[m];
                if (val > 0 && val != modal) diff_count++;
            }
            if (diff_count >= 1) kept_markers[kept_count++] = m;
        }
    }

    // 4. Print HTML Header
    fprintf(f, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n<title>STR Marker Table (Variant)</title>\n");
    fprintf(f, "<style>\n");
    fprintf(f, "body { font-family: sans-serif; margin: 20px; }\n");
    fprintf(f, "table { border-collapse: collapse; font-size: 12px; white-space: nowrap; }\n");
    fprintf(f, "th, td { border: 1px solid #ddd; text-align: center; }\n");
    fprintf(f, "td { padding: 4px; width: 20px; min-width: 20px; }\n");
    fprintf(f, "th.marker { background-color: #d9e2f3; position: sticky; top: 0; box-shadow: 0 2px 2px -1px rgba(0,0,0,0.1); height: 100px; padding: 4px; vertical-align: bottom; }\n");
    fprintf(f, "th.marker span { writing-mode: vertical-rl; transform: rotate(180deg); display: inline-block; }\n");
    fprintf(f, "td.kitid { text-align: left; padding: 4px 8px; font-weight: bold; background-color: #fafafa; position: sticky; left: 0; box-shadow: 2px 0 2px -1px rgba(0,0,0,0.1); }\n");
    fprintf(f, "th.kitid { background-color: #d9e2f3; position: sticky; left: 0; top: 0; z-index: 2; padding: 4px 8px; }\n");
    fprintf(f, "</style>\n</head>\n<body>\n");

    // 5. Build the display name dynamically using the live, raw struct variables
    char display_name[MAX_NODE_NAME_LEN * 2] = "";
    
    if (root->type == NODE_MERGED && strlen(root->gen_name) > 0) {
        // Format both the SNP half and the GEN half safely
        char snp_part[MAX_NODE_NAME_LEN];
        char gen_part[MAX_NODE_NAME_LEN];
        
        if (root->date > 0) sprintf(snp_part, "%s.%d", root->name, root->date);
        else strcpy(snp_part, root->name);
        
        if (root->gen_date > 0) sprintf(gen_part, "%s.%d", root->gen_name, root->gen_date);
        else strcpy(gen_part, root->gen_name);
        
        sprintf(display_name, "%s %s", snp_part, gen_part);
    } else {
        // Standard unmerged node formatting
        if (root->date > 0) {
            sprintf(display_name, "%s.%d", root->name, root->date);
        } else {
            strcpy(display_name, root->name);
        }
    }

    // ===============================================
    // HTML TABLE GENERATION STATEMENTS
    // ===============================================
    if (root->type == NODE_ROOT || strcmp(root->name, "ROOT") == 0) {
        fprintf(f, "<h2>Group %d - Variant STR Marker Values</h2>\n", global_group_num);
    } else {
        fprintf(f, "<h2>Group %d - Variant STR Marker Values (Node: %s)</h2>\n", global_group_num, display_name);
    }
    
    fprintf(f, "<table>\n<thead>\n<tr><th class=\"kitid\">Kit ID</th>");

    for (int m = 0; m < kept_count; m++) {
        fprintf(f, "<th class=\"marker\"><span>%s</span></th>", marker_names[kept_markers[m]]);
    }
    fprintf(f, "</tr>\n</thead>\n<tbody>\n");

    // 6. Print Kit Rows
    for (int i = 0; i < ordered_count; i++) {
        int k_idx = ordered_kits[i];
        
        // Print the plain Kit ID (No tags!)
        fprintf(f, "<tr><td class=\"kitid\">%s</td>", kits[k_idx].id);
        
        for (int j = 0; j < kept_count; j++) {
            int m = kept_markers[j];
            int val = kits[k_idx].str_values[m];
            int modal = baseline_modal[m]; 
            if (val == STR_MISSING || val <= 0) {
                fprintf(f, "<td></td>");
            } else {
                if (val > modal) {
                    fprintf(f, "<td style=\"background-color: #add8e6;\">%d</td>", val);
                } else if (val < modal) {
                    fprintf(f, "<td style=\"background-color: #ffb6c1;\">%d</td>", val);
                } else {
                    fprintf(f, "<td>%d</td>", val);
                }
            }
        }
        fprintf(f, "</tr>\n");
    }

    // 7. Print Modal Row
    fprintf(f, "<tr><td class=\"kitid\" style=\"background-color: #e8e8e8; border-top: 2px solid #333;\">MODAL</td>");
    for (int j = 0; j < kept_count; j++) {
        int m = kept_markers[j];
        int local_mod = root->local_modal[m];
        if (local_mod == STR_MISSING || local_mod <= 0) {
            fprintf(f, "<td style=\"background-color: #e8e8e8; border-top: 2px solid #333; font-weight: bold;\">-</td>");
        } else {
            fprintf(f, "<td style=\"background-color: #e8e8e8; border-top: 2px solid #333; font-weight: bold;\">%d</td>", local_mod);
        }
    }
    fprintf(f, "</tr>\n");
    fprintf(f, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose(f);
}

void generate_html_gen_table_output(const char* filename, TreeNode* root) {
    FILE* f = fopen(filename, "w");
    if (!f) return;
    fprintf(f, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n<title>GEN and SNP Status Table</title>\n");
    fprintf(f, "<style>\n");
    fprintf(f, "body { font-family: sans-serif; margin: 20px; }\n");
    fprintf(f, "table { border-collapse: collapse; font-size: 12px; white-space: nowrap; }\n");
    fprintf(f, "th, td { border: 1px solid #ddd; text-align: center; }\n");
    fprintf(f, "td { padding: 4px; width: 20px; min-width: 20px; }\n");
    fprintf(f, "th.marker { background-color: #d9e2f3; position: sticky; top: 0; box-shadow: 0 2px 2px -1px rgba(0,0,0,0.1); height: 150px; padding: 4px; vertical-align: bottom; }\n");
    fprintf(f, "th.marker span { writing-mode: vertical-rl; transform: rotate(180deg); display: inline-block; }\n");
    fprintf(f, "td.kitid { text-align: left; padding: 4px 8px; font-weight: bold; background-color: #fafafa; position: sticky; left: 0; box-shadow: 2px 0 2px -1px rgba(0,0,0,0.1); }\n");
    fprintf(f, "th.kitid { background-color: #d9e2f3; position: sticky; left: 0; top: 0; z-index: 2; padding: 4px 8px; }\n");
    fprintf(f, "</style>\n</head>\n<body>\n");
    fprintf(f, "<h2>Group %d - GEN & SNP Status Values</h2>\n", global_group_num);
    fprintf(f, "<table>\n<thead>\n<tr><th class=\"kitid\">Kit ID</th>");
    for (int g = 0; g < gen_hierarchy_count; g++) {
        if (gen_hierarchy[g].date > 0) {
            fprintf(f, "<th class=\"marker\"><span>%s.%d</span></th>", gen_hierarchy[g].ancestor_name, gen_hierarchy[g].date);
        } else {
            fprintf(f, "<th class=\"marker\"><span>%s</span></th>", gen_hierarchy[g].ancestor_name);
        }
    }
    
    for (int s = 0; s < snp_hierarchy_count; s++) {
        if (snp_hierarchy[s].mrca_date > 0) {
            fprintf(f, "<th class=\"marker\"><span>%s.%d</span></th>", snp_hierarchy[s].name, snp_hierarchy[s].mrca_date);
        } else {
            fprintf(f, "<th class=\"marker\"><span>%s</span></th>", snp_hierarchy[s].name);
        }
    }
    
    fprintf(f, "</tr>\n</thead>\n<tbody>\n");
    int ordered_kits[MAX_KITS];
    int ordered_count = 0;
    
    collect_all_html_kits(root, ordered_kits, &ordered_count);
    for (int i = 0; i < ordered_count; i++) {
        int k_idx = ordered_kits[i];
        
        int is_basal = 0;
        
        for (int b = 0; b < root->kit_count; b++) {
            if (root->kit_indices[b] == k_idx) { is_basal = 1; break; }
        }
        
        if (is_basal && (root->type == NODE_ROOT || strcmp(root->name, "ROOT") == 0)) {
            fprintf(f, "<tr><td class=\"kitid\">%s <span style=\"color:#aaa; font-weight:normal; font-size:0.75em; font-style:italic;\">(ROOT)</span></td>", kits[k_idx].id);
        } else {
            fprintf(f, "<tr><td class=\"kitid\">%s</td>", kits[k_idx].id);
        }
        for (int g = 0; g < gen_hierarchy_count; g++) {
            char status = ' ';
            for (int j = 0; j < gen_hierarchy[g].kit_count; j++) {
                if (gen_hierarchy[g].kit_indices[j] == k_idx) {
                    status = gen_hierarchy[g].kit_statuses[j];
                    break;
                }
            }
            
            if (status == '+') {
                fprintf(f, "<td style=\"background-color: #add8e6;\">+</td>");
            } else if (status == '-') {
                fprintf(f, "<td style=\"background-color: #ffb6c1;\">-</td>");
            } else if (status == '?') {
                fprintf(f, "<td>?</td>");
            } else {
                fprintf(f, "<td></td>");
            }
        }
        
        for (int s = 0; s < snp_hierarchy_count; s++) {
            SnpNode* curr = kits[k_idx].snps;
            char status = ' ';
            
            while(curr != NULL) {
                if (strcmp(curr->name, snp_hierarchy[s].name) == 0) {
                    if (curr->status == SNP_POSITIVE) status = '+';
                    else if (curr->status == SNP_NEGATIVE) status = '-';
                    else if (curr->status == SNP_UNTESTED) status = '?';
                    break;
                }
                curr = curr->next;
            }
            
            if (status == '+') {
                fprintf(f, "<td style=\"background-color: #add8e6;\">+</td>");
            } else if (status == '-') {
                fprintf(f, "<td style=\"background-color: #ffb6c1;\">-</td>");
            } else if (status == '?') {
                fprintf(f, "<td>?</td>");
            } else {
                fprintf(f, "<td></td>");
            }
        }
        
        fprintf(f, "</tr>\n");
    }
    fprintf(f, "</tbody>\n</table>\n</body>\n</html>\n");
    fclose(f);
}

void generate_modal_output(const char* filename, TreeNode* root) {
    if (!root) return;
    FILE* f = fopen(filename, "w");
    if (!f) return;
    // Line 1: The Header
    fprintf(f, "/MODAL\n");
    // Line 2: The sequentially ordered baseline values
    for (int m = 0; m < marker_count; m++) {
        int val = root->local_modal[m];
        
        if (val == STR_MISSING || val <= 0) {
            fprintf(f, "0"); // Output 0 for missing/unresolved markers
        } else {
            fprintf(f, "%d", val);
        }
        
        // Tab-delimit the values for clean formatting
        if (m < marker_count - 1) {
            fprintf(f, "\t");
        }
    }
    fprintf(f, "\n");
    
    fclose(f);
}

