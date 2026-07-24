#include "tremer.h"
#include "fallback_labels.h" // Assuming you still have this
#include <stdlib.h> 

// Explicitly set global pointers to NULL
Kit* kits = NULL;
int kit_count = 0;

int modal_values[MAX_MARKERS];
int marker_count = 0;
char (*marker_names)[MAX_STRING_LEN] = NULL; // This is a pointer to an array!

SnpTreeNode* snp_hierarchy = NULL;
int snp_hierarchy_count = 0;

GenGroup* gen_hierarchy = NULL;
int gen_hierarchy_count = 0;

TreeNode* tree_nodes = NULL; 
int tree_node_count = 0;
TreeNode* root_node = NULL;
int str_node_counter = 1; 

char project_name[MAX_STRING_LEN] = "Project";
int global_group_num = 0;
int relaxed_mode = 0;
char root_haplogroup_letter = '\0';

// =========================================================
// UPDATED SEARCH FUNCTION: Uses new dedicated GEN variables
// =========================================================
TreeNode* search_tree_for_target(TreeNode* node, const char* target) {
    if (!node) return NULL;

    // 1. Check SNP side (Base name and Date-appended name)
    char snp_display[MAX_NODE_NAME_LEN + 32];
    if (node->date > 0) {
        sprintf(snp_display, "%s.%d", node->name, node->date);
    } else {
        strcpy(snp_display, node->name);
    }

    if (strcmp(node->name, target) == 0 || strcmp(snp_display, target) == 0) return node;

    // 2. Check GEN side if this is a merged node
    if (node->type == NODE_MERGED && strlen(node->gen_name) > 0) {
        char gen_display[MAX_NODE_NAME_LEN + 32];
        if (node->gen_date > 0) {
            sprintf(gen_display, "%s.%d", node->gen_name, node->gen_date);
        } else {
            strcpy(gen_display, node->gen_name);
        }

        // If the user typed the exact genealogy alias or its dated version, we have a match
        if (strcmp(node->gen_name, target) == 0 || strcmp(gen_display, target) == 0) return node;
    }

    // 3. Match against specific kit IDs
    for (int i = 0; i < node->kit_count; i++) {
        if (strcmp(kits[node->kit_indices[i]].id, target) == 0) {
            return node; 
        }
    }

    // 4. Recurse into children
    TreeNode* c = node->first_child;
    while (c) {
        TreeNode* found = search_tree_for_target(c, target);
        if (found) return found;
        c = c->next_sibling;
    }

    return NULL;
}

void run_syntax_checks_and_export_json(const char* proj_name, int grp_num) {
    char json_filename[256];
    sprintf(json_filename, "%s_report_g%d.json", proj_name, grp_num);

    FILE* jf = fopen(json_filename, "w");
    if (!jf) return;

    char warnings[100][256];
    int w_count = 0;
    char errors[100][256];
    int e_count = 0;

    if (marker_count > 0) {
        for (int i = 0; i < kit_count; i++) {
            if (strlen(kits[i].id) == 0) continue;
            
            int valid_strs = 0;
            for (int m = 0; m < marker_count; m++) {
                if (kits[i].str_values[m] > 0) valid_strs++;
            }
            
            if (valid_strs == 0) {
                if (w_count < 100) sprintf(warnings[w_count++], "DATA: Kit '%s' has no STR data and will be ignored by the tree builder.", kits[i].id);
            }
        }
    }

    for (int i = 0; i < snp_hierarchy_count; i++) {
        if (strlen(snp_hierarchy[i].parent_name_str) > 0) {
            int found = 0;
            for (int j = 0; j < snp_hierarchy_count; j++) {
                if (strcmp(snp_hierarchy[j].name, snp_hierarchy[i].parent_name_str) == 0) {
                    found = 1; break;
                }
            }
            if (!found) {
                if (w_count < 100) sprintf(warnings[w_count++], "SNPTREE: Parent '%s' for node '%s' is not defined in the tree.", snp_hierarchy[i].parent_name_str, snp_hierarchy[i].name);
            }
        }
    }

    fprintf(jf, "{\n");
    fprintf(jf, "  \"status\": \"%s\",\n", (e_count > 0) ? "error" : (w_count > 0 ? "warning" : "success"));
    
    fprintf(jf, "  \"errors\": [\n");
    for (int i = 0; i < e_count; i++) {
        fprintf(jf, "    \"%s\"%s\n", errors[i], (i < e_count - 1) ? "," : "");
    }
    fprintf(jf, "  ],\n");

    fprintf(jf, "  \"warnings\": [\n");
    for (int i = 0; i < w_count; i++) {
        fprintf(jf, "    \"%s\"%s\n", warnings[i], (i < w_count - 1) ? "," : "");
    }
    fprintf(jf, "  ]\n");

    fprintf(jf, "}\n");
    fclose(jf);
}

/* --- Memory Management --- */

void initialize_memory() {
    kits = calloc(MAX_KITS, sizeof(Kit));
    snp_hierarchy = calloc(MAX_SNPS, sizeof(SnpTreeNode));
    gen_hierarchy = calloc(MAX_GEN_GROUPS, sizeof(GenGroup));
    tree_nodes = calloc(MAX_TREE_NODES, sizeof(TreeNode));
    marker_names = calloc(MAX_MARKERS, sizeof(char[MAX_STRING_LEN]));
    
    if (!kits || !snp_hierarchy || !gen_hierarchy || !tree_nodes || !marker_names) {
        fprintf(stderr, "FATAL ERROR: Failed to allocate heap memory. Exiting.\n");
        exit(1);
    }
}

void cleanup_memory() {
}

int main(int argc, char* argv[]) {
    initialize_memory();

    kit_count = 0;
    marker_count = 0;
    snp_hierarchy_count = 0;
    gen_hierarchy_count = 0;
    tree_node_count = 0;
    root_node = NULL;
    str_node_counter = 1;
    global_group_num = 0;
    root_haplogroup_letter = '\0';
    relaxed_mode = 0;
    allele_counts_initialized = 0;
    strcpy(project_name, "Project");

    char target_filter[MAX_STRING_LEN] = "";
    char cmd_project_name[MAX_STRING_LEN] = "";
    int cmd_group_num = -1;
    int has_cmd_modal = 0;
    int has_filter = 0;
    int has_cmd_project = 0;
    int has_cmd_group = 0;
    int input_file_idx = -1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            strncpy(target_filter, argv[i+1], MAX_STRING_LEN - 1);
            target_filter[MAX_STRING_LEN - 1] = '\0';
            has_filter = 1;
            i++; 
        } else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            strncpy(cmd_project_name, argv[i+1], MAX_STRING_LEN - 1);
            cmd_project_name[MAX_STRING_LEN - 1] = '\0';
            has_cmd_project = 1;
            i++;
        } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            cmd_group_num = atoi(argv[i+1]);
            has_cmd_group = 1;
            i++;
        } else if (strcmp(argv[i], "-m") == 0) {
            has_cmd_modal = 1;
        } else if (strcmp(argv[i], "-relaxed") == 0) {
            relaxed_mode = 1;
        } else if (argv[i][0] != '-') {
            input_file_idx = i;
            break;
        }
    }

    if (input_file_idx == -1) {
        printf("Usage: %s [Options] <input_file.txt>\n", argv[0]);
        printf("Options:\n");
        printf("  -m    output modal.txt\n");
        printf("  -n    node_or_kit_name\n");
        printf("  -g    group_name\n");
        printf("  -p    project_name\n");
        return 1;
    }

    // Secondary memory initialization (legacy safety block)
    kits = calloc(MAX_KITS, sizeof(Kit)); 
    tree_nodes = calloc(MAX_TREE_NODES, sizeof(TreeNode)); 
    gen_hierarchy = calloc(MAX_GEN_GROUPS, sizeof(GenGroup)); 
    marker_names = calloc(MAX_MARKERS, sizeof(*marker_names)); 
    snp_hierarchy = calloc(MAX_SNPS, sizeof(SnpTreeNode));
    
    if (!kits || !tree_nodes || !gen_hierarchy || !marker_names || !snp_hierarchy) {
        printf("Memory allocation failed. Please verify system RAM availability.\n");
        return 1;
    }

    FILE* file = fopen(argv[input_file_idx], "r"); 
    if (!file) return 1;

    for (int i = 0; i < MAX_MARKERS; i++) modal_values[i] = STR_MISSING;

    char* filename_only = strrchr(argv[input_file_idx], '/');
    filename_only = filename_only ? filename_only + 1 : argv[input_file_idx];

    char temp_fn[256];
    strncpy(temp_fn, filename_only, 255);
    temp_fn[255] = '\0';

    char* token_proj = strtok(temp_fn, "-_.");
    if (token_proj) {
        strncpy(project_name, token_proj, MAX_STRING_LEN - 1);
        project_name[MAX_STRING_LEN - 1] = '\0';
        char* token_grp = strtok(NULL, "-_.");
        if (token_grp) {
            int gnum = 0, has_digits = 0;
            for (int j = 0; token_grp[j] != '\0'; j++) {
                if (token_grp[j] >= '0' && token_grp[j] <= '9') {
                    gnum = gnum * 10 + (token_grp[j] - '0');
                    has_digits = 1;
                }
            }
            if (has_digits) global_group_num = gnum;
        }
    }

    char tag[MAX_STRING_LEN]; 
    while (peek_next_char(file) != EOF) { 
        if (get_next_token(file, tag)) { 
            if (strcmp(tag, "/MARKERS") == 0) parse_markers(file); 
            else if (strcmp(tag, "/MODAL") == 0) parse_modal(file); 
            else if (strcmp(tag, "/SNPTREE") == 0) parse_snptree(file); 
            else if (strcmp(tag, "/GENDATA") == 0) parse_gendata(file); 
            else if (strcmp(tag, "/STRDATA") == 0) parse_strdata(file); 
            else if (strcmp(tag, "/SNPDATA") == 0) parse_snpdata(file); 
            else if (strcmp(tag, "/KITDATA") == 0) parse_kitdata(file); 
            else if (strcmp(tag, "/GROUPS") == 0) parse_groups(file); 
            else if (strcmp(tag, "/PROJECT") == 0) {
                char line[256];
                if (get_line_rest(file, line, sizeof(line))) {
                    char* p = line;
                    while (*p == ' ' || *p == '\t') p++;
                    if (strlen(p) > 0) {
                        strncpy(project_name, p, MAX_STRING_LEN - 1);
                        project_name[MAX_STRING_LEN - 1] = '\0';
                        char* end = project_name + strlen(project_name) - 1;
                        while (end >= project_name && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '\r')) {
                            *end = '\0'; end--;
                        }
                    }
                }
            } else if (strcmp(tag, "/GROUPNUM") == 0) {
                char line[256];
                if (get_line_rest(file, line, sizeof(line))) {
                    global_group_num = atoi(line);
                }
            }
        } 
    }
    fclose(file);

    if (has_cmd_project) {
        strncpy(project_name, cmd_project_name, MAX_STRING_LEN - 1);
        project_name[MAX_STRING_LEN - 1] = '\0';
    }
    if (has_cmd_group) {
        global_group_num = cmd_group_num;
    }

    run_syntax_checks_and_export_json(project_name, global_group_num);

    compute_missing_global_modals(marker_count);
    
    for (int m = 0; m < marker_count; m++) {
        if (marker_names[m][0] == '\0') {
            if (m < n_lab && strlen(str_lab[m]) > 0) strcpy(marker_names[m], str_lab[m]);
            else sprintf(marker_names[m], "%d", m + 1);
        }
    }
    
    // Core Engine Execution
    impute_snp_statuses();
    build_skeleton_and_bucket_kits();
    initialize_modals_top_down(root_node, modal_values, marker_count); 
    refine_modals_bottom_up(root_node, marker_count);
    
    group_shared_str_mutations(root_node, marker_count); 
    apply_rule_of_two(root_node, marker_count); 
    
    reassign_root_orphans(marker_count);
    
    group_shared_str_mutations(root_node, marker_count); 
    apply_rule_of_two(root_node, marker_count); 
    
    check_kits_and_reassign(marker_count);
    
    // --- Final Quality Check Phase ---
    post_process_parsimony_scrub(root_node, marker_count);
    
    reevaluate_mutation_labels(root_node, marker_count); 
    
    void prune_empty_str_nodes_with_kits(TreeNode* node);
    collapse_empty_str_nodes(root_node);
    prune_empty_str_nodes_with_kits(root_node);
    finalize_str_node_numbering(root_node);

    sweep_orphaned_kits(root_node);

    // ==========================================================
    // SECOND PARSIMONY SCRUB
    // ==========================================================
    post_process_parsimony_scrub(root_node, marker_count);
    
    // ==========================================================
    // STR AGE ESTIMATION (Explicit 2-Pass Logic)
    // ==========================================================
    estimate_str_node_ages(root_node);

    // ==========================================================
    // EXPORT THE CALCULATED MODAL BLOCK IF REQUESTED (-m)
    // ==========================================================
    if (has_cmd_modal && root_node) {
        char modal_filename[256];
        sprintf(modal_filename, "%s_modal_g%d.txt", project_name, global_group_num);
        generate_modal_output(modal_filename, root_node);
        printf("Generated suggested MODAL block: %s\n", modal_filename);
    }
    
    if (has_filter) {
        TreeNode* target_node = search_tree_for_target(root_node, target_filter);
        char html_b_filename[512];
        
        sprintf(html_b_filename, "%s_b%d.html", project_name, global_group_num);
        
        if (target_node) {
            generate_html_filtered_table_output(html_b_filename, target_node);
            printf("Generated isolated variant table for '%s': %s\n", target_filter, html_b_filename);
        } else {
            printf("Warning: Target '%s' not found. Defaulting to full tree.\n", target_filter);
            generate_html_filtered_table_output(html_b_filename, root_node);
        }
    } else {
        char txt_filename[256]; 
        sprintf(txt_filename, "tree%d.txt", global_group_num);
        FILE* out_file = fopen(txt_filename, "w");
        if (out_file) { 
            fprintf(out_file, "Group%d - %d kits\n", global_group_num, kit_count); 
            generate_text_output(out_file, root_node, 0); 
            fclose(out_file); 
        }
        
        char html_filename[256]; 
        sprintf(html_filename, "%s_s%d.html", project_name, global_group_num);
        generate_html_table_output(html_filename, root_node);

        char html_b_filename[256]; 
        sprintf(html_b_filename, "%s_b%d.html", project_name, global_group_num);
        generate_html_filtered_table_output(html_b_filename, root_node);

        char html_g_filename[256]; 
        sprintf(html_g_filename, "%s_g%d.html", project_name, global_group_num);
        generate_html_gen_table_output(html_g_filename, root_node);

        char svg_filename[256]; 
        sprintf(svg_filename, "%s_s%d.svg", project_name, global_group_num);
        generate_svg_output(svg_filename, root_node);
        
        // printf("Full tree output generated.\n");
    }
    
    for (int i = 0; i < kit_count; i++) {
        SnpNode* current = kits[i].snps;
        while (current != NULL) { SnpNode* next = current->next; free(current); current = next; }
    }
    free(kits); free(tree_nodes); free(gen_hierarchy); free(marker_names); free(snp_hierarchy);
    cleanup_memory();
    return 0;
}

