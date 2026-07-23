#include "tremer.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MIN_SHARED            2
#define MIN_KITS_CLUSTER_ROOT 2
#define MIN_KITS_CLUSTER_NODE 4

int calculate_str_distance(int strA[], int strB[], int total_markers);
void branch_by_shared_mutations(TreeNode* node, int total_markers);
void cluster_inferred_str_nodes(TreeNode* node, int total_markers);

/* --- 1. Utility and Search Functions --- */

SnpTreeNode* find_snp_by_name(const char* name) {
    for (int i = 0; i < snp_hierarchy_count; i++) if (strcmp(snp_hierarchy[i].name, name) == 0) return &snp_hierarchy[i];
    return NULL;
}

int get_snp_depth(SnpTreeNode* snp) {
    int depth = 0; while (snp != NULL) { depth++; snp = snp->parent; } return depth;
}

SnpTreeNode* find_raw_terminal_snp(Kit* kit) {
    SnpTreeNode* terminal = NULL; int max_depth = -1;
    SnpNode* current_snp = kit->snps;
    while (current_snp != NULL) {
        if (current_snp->status == SNP_POSITIVE) {
            SnpTreeNode* snp_ref = find_snp_by_name(current_snp->name);
            if (snp_ref != NULL) { int depth = get_snp_depth(snp_ref); if (depth > max_depth) { max_depth = depth; terminal = snp_ref; } }
        }
        current_snp = current_snp->next;
    }
    return terminal;
}

int is_snp_ancestor(SnpTreeNode* ancestor, SnpTreeNode* descendant) {
    if (ancestor == NULL) return 0;
    SnpTreeNode* curr = descendant;
    while (curr != NULL) { if (curr == ancestor) return 1; curr = curr->parent; }
    return 0;
}

SnpTreeNode* get_gen_group_mrca_snp(int gen_idx) {
    SnpTreeNode* deepest = NULL;
    int max_depth = -1;
    for (int k = 0; k < gen_hierarchy[gen_idx].kit_count; k++) {
        char st = gen_hierarchy[gen_idx].kit_statuses[k];
        if (st == '-') continue; 
        int kit_idx = gen_hierarchy[gen_idx].kit_indices[k];
        SnpTreeNode* term = find_raw_terminal_snp(&kits[kit_idx]);
        if (term != NULL) {
            int depth = get_snp_depth(term);
            if (depth > max_depth) { max_depth = depth; deepest = term; }
        }
    }
    return deepest;
}

SnpTreeNode* find_imputed_terminal_snp(Kit* kit, int kit_idx) {
    SnpTreeNode* terminal = find_raw_terminal_snp(kit);
    if (terminal != NULL) return terminal;
    int best_gen_idx = -1; int min_size = 999999;
    for (int g = 0; g < gen_hierarchy_count; g++) {
        int in_gen = 0; char status = ' ';
        for (int k = 0; k < gen_hierarchy[g].kit_count; k++) {
            if (gen_hierarchy[g].kit_indices[k] == kit_idx) { in_gen = 1; status = gen_hierarchy[g].kit_statuses[k]; break; }
        }
        if (in_gen && (status == '+' || status == ' ') && gen_hierarchy[g].kit_count < min_size) {
            SnpTreeNode* mrca = get_gen_group_mrca_snp(g);
            if (mrca != NULL) { min_size = gen_hierarchy[g].kit_count; best_gen_idx = g; }
        }
    }
    if (best_gen_idx != -1) {
        SnpTreeNode* mrca = get_gen_group_mrca_snp(best_gen_idx);
        while (mrca != NULL) {
            int conflict = 0; SnpTreeNode* check_node = mrca;
            while (check_node) {
                SnpNode* sn = kit->snps;
                while (sn) {
                    if (strcmp(sn->name, check_node->name) == 0 && sn->status == SNP_NEGATIVE) { conflict = 1; break; }
                    sn = sn->next;
                }
                if (conflict) break;
                check_node = check_node->parent;
            }
            if (!conflict) return mrca; 
            mrca = mrca->parent; 
        }
    }
    return NULL;
}

void impute_snp_statuses() {
    for (int i = 0; i < snp_hierarchy_count; i++) {
        if (strlen(snp_hierarchy[i].parent_name_str) > 0) {
            snp_hierarchy[i].parent = find_snp_by_name(snp_hierarchy[i].parent_name_str);
        }
    }

    for (int k = 0; k < kit_count; k++) {
        if (kits[k].id[0] == '\0') continue;
        
        for (int s = 0; s < snp_hierarchy_count; s++) {
            SnpTreeNode* snp_node = &snp_hierarchy[s];
            SnpTreeNode* curr = snp_node->parent;
            int is_descendant_of_neg = 0;
            
            while (curr != NULL) {
                SnpNode* sn = kits[k].snps;
                int found_neg = 0;
                while (sn) {
                    if (strcmp(sn->name, curr->name) == 0 && sn->status == SNP_NEGATIVE) { found_neg = 1; break; }
                    sn = sn->next;
                }
                if (found_neg) { is_descendant_of_neg = 1; break; }
                curr = curr->parent;
            }
            
            if (is_descendant_of_neg) {
                SnpNode* sn = kits[k].snps;
                int found = 0;
                while(sn) {
                    if (strcmp(sn->name, snp_node->name) == 0) { sn->status = SNP_NEGATIVE; found = 1; break; }
                    sn = sn->next;
                }
                if (!found) {
                    SnpNode* new_sn = calloc(1, sizeof(SnpNode));
                    strcpy(new_sn->name, snp_node->name);
                    new_sn->status = SNP_NEGATIVE;
                    new_sn->next = kits[k].snps;
                    kits[k].snps = new_sn;
                }
            }
        }

        SnpTreeNode* term = find_imputed_terminal_snp(&kits[k], k);
        if (term != NULL) {
            for (int s = 0; s < snp_hierarchy_count; s++) {
                SnpTreeNode* snp_node = &snp_hierarchy[s];
                int is_ancestor = is_snp_ancestor(snp_node, term);
                int is_descendant = is_snp_ancestor(term, snp_node);
                
                SnpStatus implied_status = SNP_UNTESTED;
                if (is_ancestor) implied_status = SNP_POSITIVE; 
                else if (!is_descendant) implied_status = SNP_NEGATIVE; 

                if (implied_status != SNP_UNTESTED) {
                    SnpNode* sn = kits[k].snps;
                    int found = 0;
                    while(sn) {
                        if (strcmp(sn->name, snp_node->name) == 0) {
                            if (sn->status == SNP_UNTESTED) sn->status = implied_status;
                            found = 1; break;
                        }
                        sn = sn->next;
                    }
                    if (!found) {
                        SnpNode* new_sn = calloc(1, sizeof(SnpNode));
                        strcpy(new_sn->name, snp_node->name);
                        new_sn->status = implied_status;
                        new_sn->next = kits[k].snps;
                        kits[k].snps = new_sn;
                    }
                }
            }
        }
    }
}

/* --- 2. Shield and Skeletal Rules --- */

int is_kit_negative_for_specific_node(int kit_idx, TreeNode* node) {
    if (!node || node->type == NODE_ROOT) return 0;
    char temp[MAX_NODE_NAME_LEN]; strcpy(temp, node->name);
    char* tok = strtok(temp, " ");
    while (tok) {
        char base[MAX_NODE_NAME_LEN]; strcpy(base, tok);
        char* dot = strchr(base, '.'); if (dot) *dot = '\0';
        SnpNode* snp = kits[kit_idx].snps;
        while (snp) {
            if (strcmp(snp->name, base) == 0 && snp->status == SNP_NEGATIVE) return 1; 
            snp = snp->next;
        }
        tok = strtok(NULL, " ");
    }
    return 0;
}

int is_kit_negative_for_node(int kit_idx, TreeNode* node) {
    if (!node || node->type == NODE_ROOT) return 0;
    TreeNode* curr = node;
    while (curr != NULL && curr->type != NODE_ROOT) {
        if (is_kit_negative_for_specific_node(kit_idx, curr)) return 1;
        curr = curr->parent; 
    }
    return 0;
}

int is_tree_node_ancestor(TreeNode* ancestor, TreeNode* descendant) {
    if (!ancestor || !descendant) return 0;
    TreeNode* curr = descendant; while (curr) { if (curr == ancestor) return 1; curr = curr->parent; }
    return 0;
}

void find_skeletal_node_by_name_and_date_recursive(TreeNode* node, const char* clean_base, int target_date, TreeNode** found) {
    if (!node || *found) return;
    if (node->type == NODE_SNP || node->type == NODE_GEN || node->type == NODE_MERGED) {
        char temp[MAX_NODE_NAME_LEN]; strcpy(temp, node->name);
        char* tok = strtok(temp, " ");
        while (tok) {
            char base[MAX_NODE_NAME_LEN]; strcpy(base, tok);
            char* dot = strchr(base, '.'); if (dot) *dot = '\0';
            if (strcmp(base, clean_base) == 0) {
                if (node->type == NODE_GEN) {
                    if (target_date <= 0 || node->date == target_date) { *found = node; return; }
                } else { *found = node; return; }
            }
            tok = strtok(NULL, " ");
        }
    }
    TreeNode* c = node->first_child; while (c) { find_skeletal_node_by_name_and_date_recursive(c, clean_base, target_date, found); c = c->next_sibling; }
}

TreeNode* find_skeletal_node_by_name_and_date(const char* base_name, int target_date) {
    char clean_base[MAX_NODE_NAME_LEN]; strcpy(clean_base, base_name);
    char* query_dot = strchr(clean_base, '.'); if (query_dot) *query_dot = '\0';
    TreeNode* found = NULL; find_skeletal_node_by_name_and_date_recursive(root_node, clean_base, target_date, &found);
    return found;
}

int is_acceptable_target(int kit_idx, TreeNode* target_node) {
    if (!target_node) return 0;
    
    // 1. Find the deepest explicit GEN assignment for this kit
    TreeNode* explicit_gen_node = NULL;
    for (int g = 0; g < gen_hierarchy_count; g++) {
        for (int i = 0; i < gen_hierarchy[g].kit_count; i++) {
            if (gen_hierarchy[g].kit_indices[i] == kit_idx) {
                char status = gen_hierarchy[g].kit_statuses[i];
                TreeNode* gen_anchor = find_skeletal_node_by_name_and_date(gen_hierarchy[g].ancestor_name, gen_hierarchy[g].date);
                if (gen_anchor) {
                    if (status == '+' || status == ' ') {
                        // FIX 1: Only update if the new anchor is deeper in the tree
                        if (!explicit_gen_node || is_tree_node_ancestor(explicit_gen_node, gen_anchor)) {
                            explicit_gen_node = gen_anchor;
                        }
                    }
                    else if (status == '-' && is_tree_node_ancestor(gen_anchor, target_node)) {
                        return 0;
                    }
                }
            }
        }
    }
    
    // 2. Paper-Trail Lock: A kit MUST sit on its explicit GEN node or a descendant of it
    if (explicit_gen_node && !is_tree_node_ancestor(explicit_gen_node, target_node)) return 0;
    
    // 3. Strict Genealogical Boundary: A kit cannot cross into a deeper explicit GEN node's territory.
    if (!relaxed_mode && explicit_gen_node && explicit_gen_node != target_node) {
        TreeNode* curr = target_node;
        while (curr != NULL && curr != explicit_gen_node) {
            if (curr->type == NODE_GEN || curr->type == NODE_SNP) {
                return 0;
            }
            curr = curr->parent;
        }
    }
    
    // 3. Evaluate SNP constraints
    SnpNode* sn = kits[kit_idx].snps;
    while (sn) {
        TreeNode* anchor = find_skeletal_node_by_name_and_date(sn->name, -1);
        if (anchor) {
            // Negative SNP constraint is absolute
            if (sn->status == SNP_NEGATIVE && is_tree_node_ancestor(anchor, target_node)) return 0;
            
            // Positive SNP constraint
            if (sn->status == SNP_POSITIVE && !is_tree_node_ancestor(anchor, target_node)) {
                
                // FIX 2: If the positive SNP is parallel to the explicit GEN node, 
                // the paper-trail overrides the SNP constraint.
                int bypass = 0;
                if (explicit_gen_node) {
                    int is_parallel = !is_tree_node_ancestor(explicit_gen_node, anchor) && !is_tree_node_ancestor(anchor, explicit_gen_node);
                    if (is_parallel) {
                        bypass = 1;
                    }
                }
                if (!bypass) return 0;
            }
        }
        sn = sn->next;
    }
    return 1;
}

/* --- 3. Date Helpers --- */

int get_min_known_descendant_date(TreeNode* node) {
    if (!node) return 999999;
    int min_date = 999999;
    TreeNode* c = node->first_child;
    while (c) {
        if ((c->type == NODE_SNP || c->type == NODE_GEN || c->type == NODE_MERGED) && c->date > 0 && c->date < min_date) {
            min_date = c->date;
        }
        int child_min = get_min_known_descendant_date(c);
        if (child_min > 0 && child_min < min_date) min_date = child_min;
        c = c->next_sibling;
    }
    return min_date;
}

void impute_missing_snp_dates_proportional(TreeNode* node, int total_markers) {
    if (!node) return;
    if (node->type == NODE_ROOT) {
        TreeNode* c = node->first_child;
        while (c) { impute_missing_snp_dates_proportional(c, total_markers); c = c->next_sibling; }
        return;
    }
    if ((node->type == NODE_SNP || node->type == NODE_GEN || node->type == NODE_MERGED) && node->date <= 0 && node->parent) {
        int parent_date = node->parent->date;
        if (parent_date > 0) {
            int child_date = get_min_known_descendant_date(node);
            if (child_date == 999999 || child_date <= parent_date) node->date = parent_date + 30;
            else node->date = parent_date + ((child_date - parent_date) / 2);
        }
    }
    TreeNode* c = node->first_child;
    while (c) { impute_missing_snp_dates_proportional(c, total_markers); c = c->next_sibling; }
}

void calculate_str_node_ages_proportional(TreeNode* node, int total_markers) {
    if (!node) return;
    if ((node->type == NODE_STR || node->type == NODE_STR_BRANCH) && node->parent && node->parent->date > 0) {
        int parent_date = node->parent->date;
        float m_up = 0;
        for (int m = 0; m < total_markers; m++) {
            if (node->parent->local_modal[m] > 0 && node->local_modal[m] > 0 && node->parent->local_modal[m] != node->local_modal[m]) m_up++;
        }
        int child_date = get_min_known_descendant_date(node);
        if (child_date == 999999 || child_date <= parent_date) node->date = parent_date + (int)(m_up * 30);
        else node->date = parent_date + ((child_date - parent_date) / 2);
    }
    TreeNode* c = node->first_child;
    while (c) { calculate_str_node_ages_proportional(c, total_markers); c = c->next_sibling; }
}

/* --- 4. Basic Tree Building --- */

int get_node_priority(NodeType type) {
    if (type == NODE_DISTANCE) return 0; 
    if (type == NODE_GEN) return 1;
    if (type == NODE_SNP || type == NODE_MERGED) return 2; 
    if (type == NODE_STR_BRANCH) return 3; 
    if (type == NODE_STR) return 4;
    return 5; 
}

void add_child_to_node(TreeNode* parent, TreeNode* child) {
    child->parent = parent; child->next_sibling = NULL;
    if (parent->first_child == NULL) { parent->first_child = child; return; }
    int child_priority = get_node_priority(child->type);
    if (get_node_priority(parent->first_child->type) > child_priority) { child->next_sibling = parent->first_child; parent->first_child = child; return; }
    TreeNode* current = parent->first_child;
    while (current->next_sibling != NULL && get_node_priority(current->next_sibling->type) <= child_priority) current = current->next_sibling;
    child->next_sibling = current->next_sibling; current->next_sibling = child;
}

int count_intersection(int* setA, int countA, int* setB, int countB) {
    int matches = 0;
    for (int i = 0; i < countA; i++) for (int j = 0; j < countB; j++) if (setA[i] == setB[j]) { matches++; break; }
    return matches;
}

void collect_subtree_kits(TreeNode* node, int* collected_indices, int* count) {
    if (node == NULL) return;
    for (int i = 0; i < node->kit_count; i++) {
        if (*count < MAX_KITS) collected_indices[(*count)++] = node->kit_indices[i];
    }
    TreeNode* child = node->first_child;
    while (child != NULL) { collect_subtree_kits(child, collected_indices, count); child = child->next_sibling; }
}

void post_process_age_correction(TreeNode* node, int parent_date) {
    if (!node) return;
    if (node->type == NODE_SNP && node->date > 0 && parent_date > 0 && node->date < parent_date) node->date = parent_date;
    int effective_date = (node->date > 0) ? node->date : parent_date;
    TreeNode* child = node->first_child;
    while (child) { post_process_age_correction(child, effective_date); child = child->next_sibling; }
}

int force_snp_ages_bottom_up(TreeNode* node) {
    if (!node) return 999999;
    int min_child_date = 999999;
    TreeNode* c = node->first_child;
    while (c) {
        int c_date = force_snp_ages_bottom_up(c);
        if (c_date > 0 && c_date < min_child_date) min_child_date = c_date;
        c = c->next_sibling;
    }
    if (node->type == NODE_SNP && node->date > 0 && min_child_date != 999999 && node->date > min_child_date) node->date = min_child_date; 
    return (node->date > 0) ? node->date : min_child_date;
}

void post_process_merges(TreeNode* node) {
    if (!node) return;
    int merged_any = 1;
    while (merged_any) {
        merged_any = 0; TreeNode* prev = NULL; TreeNode* curr = node->first_child;
        while (curr) {
            if (node->type == NODE_SNP && curr->type == NODE_GEN) {
                int intersect = count_intersection(node->defined_kits, node->defined_kit_count, curr->defined_kits, curr->defined_kit_count);
                
                // 1. Calculate Date Gap
                int age_gap = 0;
                if (node->date > 0 && curr->date > 0) {
                    age_gap = curr->date - node->date;
                    if (age_gap < 0) age_gap = -age_gap; // Get absolute difference
                }

                // 2. Overlap vs Subset Logic
                int meets_overlap = ((intersect * 4 >= node->defined_kit_count * 3) && (intersect * 4 >= curr->defined_kit_count * 3));
                int is_pure_subset = (intersect == node->defined_kit_count && intersect > 0); 

                if (meets_overlap || is_pure_subset) {
                    
                    // 3. Conflict Checks
                    int conflict = 0;
                    for (int k = 0; k < node->defined_kit_count; k++) {
                        if (is_kit_negative_for_specific_node(node->defined_kits[k], curr)) { conflict = 1; break; }
                    }
                    
                    if (!conflict && !is_pure_subset) {
                        for (int k = 0; k < curr->defined_kit_count; k++) {
                            if (is_kit_negative_for_specific_node(curr->defined_kits[k], node)) { conflict = 1; break; }
                        }
                    }

                    // 4. Final Decision
                    if (!conflict && (node->date <= 0 || curr->date <= 0 || age_gap <= 30)) {
                        node->type = NODE_MERGED;
                        
                        // NEW LOGIC: Store the GEN data in its own dedicated fields!
                        strcpy(node->gen_name, curr->name);
                        node->gen_date = curr->date;
                        
                        // The SNP data (node->name and node->date) remains untouched!
                        
                        for (int k = 0; k < curr->defined_kit_count; k++) {
                            int exists = 0; for (int c = 0; c < node->defined_kit_count; c++) if (node->defined_kits[c] == curr->defined_kits[k]) { exists = 1; break; }
                            if (!exists && node->defined_kit_count < MAX_KITS) node->defined_kits[node->defined_kit_count++] = curr->defined_kits[k];
                        }
                        
                        if (prev == NULL) node->first_child = curr->next_sibling; else prev->next_sibling = curr->next_sibling;
                        TreeNode* ch = curr->first_child;
                        while(ch) { TreeNode* nxt = ch->next_sibling; ch->next_sibling = NULL; add_child_to_node(node, ch); ch = nxt; }
                        merged_any = 1; break; 
                    }
                }
            }
            prev = curr; curr = curr->next_sibling;
        }
    }
    TreeNode* child = node->first_child; while (child) { post_process_merges(child); child = child->next_sibling; }
}

void build_skeleton_and_bucket_kits() {
    for (int k = 0; k < kit_count; k++) {
        int valid = 0; for (int m = 0; m < marker_count; m++) if (kits[k].str_values[m] > 0) valid++;
        if (valid == 0) kits[k].id[0] = '\0';
    }
    
    TreeNode** groups = calloc(MAX_SNPS + MAX_GEN_GROUPS, sizeof(TreeNode*)); int group_count = 0;
    
    for (int i = 0; i < snp_hierarchy_count; i++) {
        TreeNode* n = &tree_nodes[tree_node_count++]; n->type = NODE_SNP; strcpy(n->name, snp_hierarchy[i].name); n->date = snp_hierarchy[i].mrca_date;
        for (int k = 0; k < kit_count; k++) {
            if (kits[k].id[0] == '\0') continue; 
            SnpTreeNode* curr = find_imputed_terminal_snp(&kits[k], k);
            while(curr) { if (strcmp(curr->name, snp_hierarchy[i].name) == 0) { if (n->defined_kit_count < MAX_KITS) n->defined_kits[n->defined_kit_count++] = k; break; } curr = curr->parent; }
        }
        groups[group_count++] = n;
    }
    
    for (int i = 0; i < gen_hierarchy_count; i++) {
        TreeNode* n = &tree_nodes[tree_node_count++]; n->type = NODE_GEN; strcpy(n->name, gen_hierarchy[i].ancestor_name); n->date = gen_hierarchy[i].date;
        for (int k = 0; k < gen_hierarchy[i].kit_count; k++) {
            char status = gen_hierarchy[i].kit_statuses[k];
            if (status == '+' || status == ' ') if (n->defined_kit_count < MAX_KITS) n->defined_kits[n->defined_kit_count++] = gen_hierarchy[i].kit_indices[k];
        }
        groups[group_count++] = n;
    }
    
    for (int i = 0; i < group_count - 1; i++) {
        for (int j = i + 1; j < group_count; j++) {
            int swap = 0;
            if (groups[j]->defined_kit_count > groups[i]->defined_kit_count) {
                swap = 1;
            }
            else if (groups[j]->defined_kit_count == groups[i]->defined_kit_count) {
                if (groups[j]->type == NODE_SNP && groups[i]->type == NODE_GEN) {
                    int gen_has_neg = 0;
                    for (int k = 0; k < groups[i]->defined_kit_count; k++) {
                        if (is_kit_negative_for_specific_node(groups[i]->defined_kits[k], groups[j])) { gen_has_neg = 1; break; }
                    }
                    if (!gen_has_neg) swap = 1;
                }
                // THE FIX: Add a tie-breaker for phylogenetically equivalent SNP blocks
                else if (groups[j]->type == NODE_SNP && groups[i]->type == NODE_SNP) {
                    SnpTreeNode* snp_i = find_snp_by_name(groups[i]->name);
                    SnpTreeNode* snp_j = find_snp_by_name(groups[j]->name);
                    
                    // If 'j' is the true ancestor of 'i', force it to sort first!
                    if (is_snp_ancestor(snp_j, snp_i)) {
                        swap = 1; 
                    }
                }
            }
            if (swap) { TreeNode* t = groups[i]; groups[i] = groups[j]; groups[j] = t; }
        }
    }
    
    root_node = &tree_nodes[tree_node_count++]; root_node->type = NODE_ROOT; strcpy(root_node->name, "ROOT");
    for (int k = 0; k < kit_count; k++) if (kits[k].id[0] != '\0') root_node->defined_kits[root_node->defined_kit_count++] = k;
    
    for (int i = 0; i < group_count; i++) {
        TreeNode* n = groups[i]; if (n->defined_kit_count == 0) continue;
        TreeNode* parent = root_node;
        while (1) {
            TreeNode* best_child = NULL; int max_intersect = 0; TreeNode* child = parent->first_child;
            while (child) {
                int intersect = count_intersection(n->defined_kits, n->defined_kit_count, child->defined_kits, child->defined_kit_count);
                int conflict = 0;
                
                if (child->type == NODE_SNP || child->type == NODE_MERGED) {
                    for (int k = 0; k < n->defined_kit_count; k++) {
                        if (is_kit_negative_for_specific_node(n->defined_kits[k], child)) { conflict = 1; break; }
                    }
                }

                if (!conflict && intersect > max_intersect) { max_intersect = intersect; best_child = child; }
                child = child->next_sibling;
            }
            
            if (best_child) {
                int threshold = n->defined_kit_count / 2;
                if (max_intersect > threshold) parent = best_child;
                else break;
            } else break;
        }
        if (n != NULL) add_child_to_node(parent, n);
    }
    
    post_process_age_correction(root_node, 0); force_snp_ages_bottom_up(root_node); post_process_merges(root_node);
    
    for (int k = 0; k < kit_count; k++) {
        if (kits[k].id[0] == '\0') continue; 
        TreeNode* curr = root_node;
        while(1) {
            TreeNode* found = NULL; TreeNode* c = curr->first_child;
            while(c) { 
                for(int i = 0; i < c->defined_kit_count; i++) if (c->defined_kits[i] == k && !is_kit_negative_for_node(k, c)) { found = c; break; } 
                if (found) break; c = c->next_sibling; 
            }
            if (found) curr = found; else break;
        }
        if (curr->kit_count < MAX_KITS) curr->kit_indices[curr->kit_count++] = k;
    }
    free(groups);
}

/* --- 5. Modals and STR Handling --- */

void compute_missing_global_modals(int total_markers) {
    for (int m = 0; m < total_markers; m++) {
        if (modal_values[m] <= 0) {
            int freqs[256] = {0}, max_freq = 0, best_val = STR_MISSING;
            for (int i = 0; i < kit_count; i++) { 
                if (kits[i].id[0] == '\0') continue; 
                int val = kits[i].str_values[m]; 
                if (val > 0 && val < 256) { freqs[val]++; if (freqs[val] > max_freq) { max_freq = freqs[val]; best_val = val; } } 
            }
            modal_values[m] = best_val;
        }
    }
}

void estimate_node_modal(TreeNode* node, int total_markers) {
    if (node->type == NODE_ROOT) { for (int m = 0; m < total_markers; m++) node->local_modal[m] = modal_values[m]; return; }
    for (int m = 0; m < total_markers; m++) {
        int freq[256] = {0}, max_freq = 0;
        for (int i = 0; i < node->kit_count; i++) { 
            int val = kits[node->kit_indices[i]].str_values[m]; if (val > 0 && val < 256) { freq[val]++; if (freq[val] > max_freq) max_freq = freq[val]; } 
        }
        TreeNode* c = node->first_child;
        while (c) { if (c->type != NODE_DISTANCE) { int val = c->local_modal[m]; if (val > 0 && val < 256) { freq[val]++; if (freq[val] > max_freq) max_freq = freq[val]; } } c = c->next_sibling; }
        int p_modal = (node->parent != NULL) ? node->parent->local_modal[m] : modal_values[m]; 
        int mode_val = STR_MISSING;
        if (max_freq > 0) { 
            int retain = (node->type == NODE_STR_BRANCH || node->type == NODE_STR) ? 1 : 2;
            if (p_modal > 0 && p_modal < 256 && freq[p_modal] >= retain) mode_val = p_modal; 
            else { for (int v = 0; v < 256; v++) if (freq[v] == max_freq) { mode_val = v; break; } } 
        }
        node->local_modal[m] = (mode_val > 0) ? mode_val : p_modal;
    }
}

void initialize_modals_top_down(TreeNode* node, int* parent_modal, int total_markers) {
    if (node == NULL) return;
    for (int m = 0; m < total_markers; m++) node->local_modal[m] = parent_modal[m];
    TreeNode* child = node->first_child; while (child != NULL) { initialize_modals_top_down(child, node->local_modal, total_markers); child = child->next_sibling; }
}

void refine_modals_bottom_up(TreeNode* node, int marker_count) {
    if (!node) return;
    TreeNode* child = node->first_child; while (child) { refine_modals_bottom_up(child, marker_count); child = child->next_sibling; }
    if (node->kit_count == 0 && node->first_child == NULL) return;
    for (int m = 0; m < marker_count; m++) {
        int counts[256] = {0}, max_count = 0, best = STR_MISSING;
        for (int i = 0; i < node->kit_count; i++) {
            int val = kits[node->kit_indices[i]].str_values[m]; if (val > 0 && val < 256) { counts[val]++; if (counts[val] > max_count) { max_count = counts[val]; best = val; } }
        }
        TreeNode* c = node->first_child;
        while (c) { if (c->type != NODE_DISTANCE) { int val = c->local_modal[m]; if (val > 0 && val < 256) { counts[val]++; if (counts[val] > max_count) { max_count = counts[val]; best = val; } } } c = c->next_sibling; }
        if (node->local_modal[m] > 0 && counts[node->local_modal[m]] == max_count) best = node->local_modal[m];
        if (best != STR_MISSING && max_count > 0) node->local_modal[m] = best;
    }
}

void group_shared_str_mutations(TreeNode* node, int total_markers) {
    if (node == NULL) return;
    
    // THE SHIELD: Prevent STR branches from forming inside established genealogical or biological nodes
    

    int grouped = 1;
    while (grouped) {
        grouped = 0; int child_node_count = 0; TreeNode* c = node->first_child; while (c) { child_node_count++; c = c->next_sibling; }
        if (child_node_count < 2) break; 
        int best_m = -1, best_val = -1, max_shared = 1;
        for (int m = 0; m < total_markers; m++) {
            int p_val = node->local_modal[m]; 
            if (p_val <= 0) continue;
            int freqs[256] = {0}, child_node_matches[256] = {0}; 
            c = node->first_child; while (c) { int cv = c->local_modal[m]; if (cv > 0 && cv != p_val && cv < 256) { freqs[cv]++; child_node_matches[cv]++; } c = c->next_sibling; }
            for (int i = 0; i < node->kit_count; i++) { 
                if (relaxed_mode || (node->type != NODE_GEN && node->type != NODE_SNP)) {
                    int kv = kits[node->kit_indices[i]].str_values[m]; if (kv > 0 && kv != p_val && kv < 256) freqs[kv]++; 
                }
            }
            for (int v = 0; v < 256; v++) {
                if (child_node_matches[v] >= 2 && freqs[v] > max_shared) { best_m = m; best_val = v; max_shared = freqs[v]; }
            }
        }
        if (best_m != -1) {
            TreeNode* strn = &tree_nodes[tree_node_count++]; strn->type = NODE_STR_BRANCH; sprintf(strn->name, "STR%02d", str_node_counter++);
            strn->parent = node; for (int i = 0; i < MAX_MARKERS; i++) strn->local_modal[i] = node->local_modal[i]; strn->local_modal[best_m] = best_val;
            strn->mutations[0].marker_index = best_m;
            strn->mutations[0].old_val = node->local_modal[best_m];
            strn->mutations[0].new_val = best_val;
            strn->mutation_count = 1;
            TreeNode* prev = NULL; c = node->first_child;
            int g_cnt = 0;
            while (c != NULL) { 
                TreeNode* nxt = c->next_sibling; 
                if (c->local_modal[best_m] == best_val) { if (prev == NULL) node->first_child = nxt; else prev->next_sibling = nxt; c->next_sibling = NULL; add_child_to_node(strn, c); g_cnt++; } else prev = c;
                c = nxt; 
            }
            if (strcmp(node->name, "WilliamNeely.1591") == 0) {
                printf("DEBUG: WilliamNeely grouped %d children into %s for marker %s\n", g_cnt, strn->name, marker_names[best_m]);
            }
            if (relaxed_mode || (node->type != NODE_GEN && node->type != NODE_SNP)) {
                int k_k[MAX_KITS], k_c = 0;
                for (int i = 0; i < node->kit_count; i++) { 
                    int k_i = node->kit_indices[i]; if (kits[k_i].str_values[best_m] == best_val) { if (strn->kit_count < MAX_KITS) strn->kit_indices[strn->kit_count++] = k_i; } 
                    else { if (k_c < MAX_KITS) k_k[k_c++] = k_i; } 
                }
                node->kit_count = k_c; for (int i = 0; i < k_c; i++) node->kit_indices[i] = k_k[i];
            }
            add_child_to_node(node, strn); estimate_node_modal(strn, total_markers); grouped = 1;
        }
    }
    TreeNode* child = node->first_child; while (child != NULL) { group_shared_str_mutations(child, total_markers); child = child->next_sibling; }
}

void apply_rule_of_two(TreeNode* node, int total_markers) {
    if (node == NULL) return;
    
    // THE SHIELD: Prevent STR branches from forming inside established genealogical or biological nodes
    

    if (node->kit_count >= 2) {
        int best_m = -1, best_v = -1, max_s = 1;
        for (int m = 0; m < total_markers; m++) {
            int p_val = node->local_modal[m]; if (p_val <= 0) continue;
            int freqs[256] = {0}; for (int i = 0; i < node->kit_count; i++) { int v = kits[node->kit_indices[i]].str_values[m]; if (v > 0 && v != p_val && v < 256) freqs[v]++; } 
            for (int v = 0; v < 256; v++) if (freqs[v] > max_s) { max_s = freqs[v]; best_m = m; best_v = v; }
        }
        if (best_m != -1) {
            TreeNode* str_n = &tree_nodes[tree_node_count++]; str_n->type = NODE_STR; sprintf(str_n->name, "STR%02d", str_node_counter++);
            for (int i = 0; i < MAX_MARKERS; i++) str_n->local_modal[i] = node->local_modal[i]; str_n->local_modal[best_m] = best_v;
            add_child_to_node(node, str_n); int k_k[MAX_KITS], k_c = 0;
            for (int i = 0; i < node->kit_count; i++) { 
                int k_i = node->kit_indices[i]; if (kits[k_i].str_values[best_m] == best_v) { if (str_n->kit_count < MAX_KITS) str_n->kit_indices[str_n->kit_count++] = k_i; } 
                else { if (k_c < MAX_KITS) k_k[k_c++] = k_i; } 
            }
            node->kit_count = k_c; for (int i = 0; i < k_c; i++) node->kit_indices[i] = k_k[i];
            apply_rule_of_two(str_n, total_markers); apply_rule_of_two(node, total_markers); return;
        }
    }
    TreeNode* child = node->first_child; while (child != NULL) { apply_rule_of_two(child, total_markers); child = child->next_sibling; }
}

void post_process_parsimony_scrub(TreeNode* node, int total_markers) {
    if (!node) return;
    // 1. Recurse bottom-up (process children first), passing total_markers down
    TreeNode* c = node->first_child;
    while (c) {
        TreeNode* next = c->next_sibling;
        post_process_parsimony_scrub(c, total_markers);
        c = next;
    }
    // 2. Check if THIS node is a useless pass-through
    // Condition: Node has exactly 1 child, and 0 kits of its own.
    if (node->kit_count == 0 && node->first_child != NULL && node->first_child->next_sibling == NULL) {
        TreeNode* single_child = node->first_child;
        // Only compress if the child is an inferred STR node. 
        // (We don't want to accidentally erase user-defined SNP nodes!)
        if (single_child->type == NODE_STR || single_child->type == NODE_STR_BRANCH) {
            // A. Move the kits up to the parent
            node->kit_count = single_child->kit_count;
            for (int i = 0; i < single_child->kit_count; i++) {
                node->kit_indices[i] = single_child->kit_indices[i];
            }
            node->defined_kit_count = single_child->defined_kit_count;
            for(int i = 0; i < single_child->defined_kit_count; i++){
                node->defined_kits[i] = single_child->defined_kits[i];
            }
            // B. Combine mutations (The STR mutations are biologically synonymous with the SNP at this point)
            for (int i = 0; i < single_child->mutation_count; i++) {
                if (node->mutation_count < MAX_STR_STACK) {
                    node->mutations[node->mutation_count++] = single_child->mutations[i];
                }
            }
            // C. Move the grandchildren up
            node->first_child = single_child->first_child;
            TreeNode* sub = node->first_child;
            while (sub) {
                sub->parent = node;
                sub = sub->next_sibling;
            }
            // D. We leave node->name, node->date, node->gen_name, and node->gen_date untouched!
            // The redundant child is now fully bypassed and structurally removed.
        }
    }
}

void reevaluate_mutation_labels(TreeNode* node, int total_markers) {
    if (!node) return; node->mutation_count = 0; 
    if (node->parent != NULL) {
        int sub_kits[MAX_KITS], sub_c = 0; collect_subtree_kits(node, sub_kits, &sub_c);
        for (int m = 0; m < total_markers; m++) {
            int pv = node->parent->local_modal[m], cv = node->local_modal[m];
            if (pv > 0 && cv > 0 && pv != cv) {
                int occ = 0; for (int k = 0; k < sub_c; k++) if (kits[sub_kits[k]].str_values[m] == cv) occ++;
                if (occ > 0 && (node->type == NODE_STR || node->type == NODE_STR_BRANCH || occ >= 2)) {
                    if (node->mutation_count < 20) { node->mutations[node->mutation_count].marker_index = m; node->mutations[node->mutation_count].old_val = pv; node->mutations[node->mutation_count].new_val = cv; node->mutation_count++; }
                }
            }
        }
    }
    TreeNode* c = node->first_child; while (c) { reevaluate_mutation_labels(c, total_markers); c = c->next_sibling; }
}

void finalize_str_node_numbering_recursive(TreeNode* node, int* counter) {
    if (!node) return; if (node->type == NODE_STR || node->type == NODE_STR_BRANCH) snprintf(node->name, MAX_NODE_NAME_LEN, "STR%02d", (*counter)++);
    TreeNode* c = node->first_child; while (c) { finalize_str_node_numbering_recursive(c, counter); c = c->next_sibling; }
}

void finalize_str_node_numbering(TreeNode* root) { 
    int counter = 1; 
    finalize_str_node_numbering_recursive(root, &counter); 
}

/* --- 6. Reassignment Helpers --- */

void find_kit_in_active_tree(TreeNode* node, int kit_idx, TreeNode** found_node, int* found_idx) {
    if (!node || *found_node) return;
    for (int i = 0; i < node->kit_count; i++) { if (node->kit_indices[i] == kit_idx) { *found_node = node; *found_idx = i; return; } }
    TreeNode* c = node->first_child; while (c) { find_kit_in_active_tree(c, kit_idx, found_node, found_idx); c = c->next_sibling; }
}

TreeNode* find_node_with_mutation(TreeNode* n, int m_idx, int old_val, int new_val, TreeNode* exclude) {
    if (n == NULL || n == exclude) return NULL;
    if (n != root_node && n->parent != NULL) {
        if (n->parent->local_modal[m_idx] == old_val && n->local_modal[m_idx] == new_val) {
            if (n->type != NODE_DISTANCE) return n;
        }
    }
    TreeNode* c = n->first_child;
    while (c) {
        TreeNode* found = find_node_with_mutation(c, m_idx, old_val, new_val, exclude);
        if (found) return found;
        c = c->next_sibling;
    }
    return NULL;
}

void evaluate_reassignment_targets_recursive(TreeNode* target, int kit_idx, int total_markers, TreeNode* current_node, TreeNode** best_node, int* min_dist) {
    if (!target) return;
    if (target->type != NODE_DISTANCE && target != current_node) {
        if (is_acceptable_target(kit_idx, target)) {
            int dist = 0;
            for (int m = 0; m < total_markers; m++) {
                int kv = kits[kit_idx].str_values[m], nv = target->local_modal[m];
                if (kv > 0 && nv > 0 && kv != nv) dist++;
            }
            if (dist < *min_dist || (relaxed_mode && dist == *min_dist && *best_node && is_tree_node_ancestor(*best_node, target))) { *min_dist = dist; *best_node = target; }
        }
    }
    TreeNode* c = target->first_child; while (c) { evaluate_reassignment_targets_recursive(c, kit_idx, total_markers, current_node, best_node, min_dist); c = c->next_sibling; }
}

void find_closest_node_recursive(TreeNode* node, int kit_idx, int total_markers, TreeNode** best_node, int* min_dist) {
    if (node == NULL) return;
    if (node != root_node && (node->type == NODE_SNP || node->type == NODE_GEN || node->type == NODE_MERGED || node->type == NODE_STR_BRANCH)) {
        if (is_acceptable_target(kit_idx, node)) { 
            int dist = 0;
            for (int m = 0; m < total_markers; m++) {
                int k_val = kits[kit_idx].str_values[m];
                int n_val = node->local_modal[m];
                if (k_val > 0 && n_val > 0 && k_val != n_val) dist++;
            }
            if (dist < *min_dist || (relaxed_mode && dist == *min_dist && *best_node && is_tree_node_ancestor(*best_node, node))) { *min_dist = dist; *best_node = node; }
        }
    }
    TreeNode* c = node->first_child;
    while (c) { find_closest_node_recursive(c, kit_idx, total_markers, best_node, min_dist); c = c->next_sibling; }
}

void reassign_root_orphans(int total_markers) {
    TreeNode* prev = NULL;
    TreeNode* c = root_node->first_child;
    while (c) {
        TreeNode* next = c->next_sibling;
        if (c->type == NODE_STR) {
            if (c->mutation_count > 0) {
                int m_idx = c->mutations[0].marker_index;
                int old_val = c->mutations[0].old_val;
                int new_val = c->mutations[0].new_val;
                TreeNode* target = find_node_with_mutation(root_node, m_idx, old_val, new_val, c);
                if (target) {
                    for (int k=0; k < c->kit_count; k++) {
                        if (target->kit_count < MAX_KITS) target->kit_indices[target->kit_count++] = c->kit_indices[k];
                    }
                    if (prev == NULL) root_node->first_child = next;
                    else prev->next_sibling = next;
                    c = next;
                    continue;
                }
            }
        }
        prev = c;
        c = next;
    }

    int orphan_kits[MAX_KITS]; int orphan_count = 0;
    for (int i = 0; i < root_node->kit_count; i++) if (orphan_count < MAX_KITS) orphan_kits[orphan_count++] = root_node->kit_indices[i];
    root_node->kit_count = 0;

    prev = NULL; c = root_node->first_child;
    while (c) {
        TreeNode* next = c->next_sibling;
        if (c->type == NODE_STR) {
            for (int i = 0; i < c->kit_count; i++) if (orphan_count < MAX_KITS) orphan_kits[orphan_count++] = c->kit_indices[i];
            if (prev == NULL) root_node->first_child = next;
            else prev->next_sibling = next;
            c = next;
            continue;
        }
        prev = c; c = next;
    }

    for (int i = 0; i < orphan_count; i++) {
        int kit_idx = orphan_kits[i];
        TreeNode* best_node = NULL; int min_dist = 999999;
        find_closest_node_recursive(root_node, kit_idx, total_markers, &best_node, &min_dist);
        if (best_node != NULL) { if (best_node->kit_count < MAX_KITS) best_node->kit_indices[best_node->kit_count++] = kit_idx; } 
        else { if (root_node->kit_count < MAX_KITS) root_node->kit_indices[root_node->kit_count++] = kit_idx; }
    }

    // --- NEW PREDICTIVE CLUSTERING (IMMUNE TO SCRUBBERS) ---
    // 1. Force the final root modal to establish the mathematical baseline
    compute_missing_global_modals(total_markers);
    initialize_modals_top_down(root_node, modal_values, total_markers);
    
    // 2. Execute the clustering on the final remaining orphans!
    branch_by_shared_mutations(root_node, total_markers);
    //cluster_inferred_str_nodes(root_node, total_markers);
    // -------------------------------------------------------
}

void prune_empty_str_nodes(TreeNode* node) {
    if (!node) return;
    TreeNode* prev = NULL; TreeNode* curr = node->first_child;
    while (curr) {
        TreeNode* next = curr->next_sibling; prune_empty_str_nodes(curr);
        if ((curr->type == NODE_STR || curr->type == NODE_STR_BRANCH) && curr->kit_count == 0 && curr->first_child == NULL) {
            if (prev == NULL) node->first_child = next; else prev->next_sibling = next;
            curr = next; continue;
        }
        prev = curr; curr = next;
    }
}

void collapse_empty_str_nodes(TreeNode* node) {
    if (!node) return; TreeNode* prev = NULL; TreeNode* curr = node->first_child;
    while (curr) {
        TreeNode* next = curr->next_sibling; collapse_empty_str_nodes(curr); 
        if ((curr->type == NODE_STR || curr->type == NODE_STR_BRANCH) && curr->kit_count == 0 && curr->mutation_count == 0) {
            TreeNode* c_child = curr->first_child;
            if (c_child) {
                TreeNode* child_iter = c_child; while (child_iter->next_sibling) { child_iter->parent = node; child_iter = child_iter->next_sibling; }
                child_iter->parent = node; if (prev == NULL) node->first_child = c_child; else prev->next_sibling = c_child;
                child_iter->next_sibling = next; prev = child_iter;
            } else { if (prev == NULL) node->first_child = next; else prev->next_sibling = next; }
            curr = next; continue;
        }
        prev = curr; curr = next;
    }
}

void prune_empty_str_nodes_with_kits(TreeNode* node) {
    if (!node) return; TreeNode* prev = NULL; TreeNode* curr = node->first_child;
    while (curr) {
        TreeNode* next = curr->next_sibling; prune_empty_str_nodes_with_kits(curr); 
        if ((curr->type == NODE_STR || curr->type == NODE_STR_BRANCH) && curr->mutation_count == 0) {
            for(int i = 0; i < curr->kit_count; i++) {
                if(node->kit_count < MAX_KITS) node->kit_indices[node->kit_count++] = curr->kit_indices[i];
            }
            for(int i = 0; i < curr->defined_kit_count; i++) {
                if(node->defined_kit_count < MAX_KITS) node->defined_kits[node->defined_kit_count++] = curr->defined_kits[i];
            }
            TreeNode* c_child = curr->first_child;
            if (c_child) {
                TreeNode* child_iter = c_child; while (child_iter->next_sibling) { child_iter->parent = node; child_iter = child_iter->next_sibling; }
                child_iter->parent = node; if (prev == NULL) node->first_child = c_child; else prev->next_sibling = c_child;
                child_iter->next_sibling = next; prev = child_iter;
            } else { if (prev == NULL) node->first_child = next; else prev->next_sibling = next; }
            curr = next; continue;
        }
        prev = curr; curr = next;
    }
}


/* --- 7. Validation and Reporter --- */

SnpTreeNode* get_static_snp_mrca(SnpTreeNode* a, SnpTreeNode* b) {
    if (!a || !b) return NULL;
    SnpTreeNode* curr_a = a;
    while (curr_a) {
        SnpTreeNode* curr_b = b;
        while (curr_b) {
            if (curr_a == curr_b) return curr_a;
            curr_b = curr_b->parent;
        }
        curr_a = curr_a->parent;
    }
    return NULL;
}

void report_gendata_inconsistencies(int total_markers) {
    int warning_count = 0;
    const int MAX_WARNINGS = 4;

    for (int g = 0; g < gen_hierarchy_count; g++) {
        char* gen_name = gen_hierarchy[g].ancestor_name;
        int g_date = gen_hierarchy[g].date;

        char full_gen_name[MAX_NODE_NAME_LEN];
        if (g_date > 0) snprintf(full_gen_name, sizeof(full_gen_name), "%s.%d", gen_name, g_date);
        else strcpy(full_gen_name, gen_name);

        TreeNode* gen_node = find_skeletal_node_by_name_and_date(gen_name, g_date);

        // 1. SNP Contradictions (Consolidated and Topologically Aware)
        int group_snp_warned = 0; 
        for (int i = 0; i < gen_hierarchy[g].kit_count && !group_snp_warned; i++) {
            if (gen_hierarchy[g].kit_statuses[i] == '+' || gen_hierarchy[g].kit_statuses[i] == ' ') {
                int kit_idx = gen_hierarchy[g].kit_indices[i];
                SnpNode* sn_neg = kits[kit_idx].snps;
                
                while (sn_neg && !group_snp_warned) {
                    if (sn_neg->status == SNP_NEGATIVE) {
                        int found_pos = 0;
                        int other_idx = -1;
                        char pos_kit_id[64] = "";
                        for (int j = 0; j < gen_hierarchy[g].kit_count; j++) {
                            if (i == j) continue;
                            if (gen_hierarchy[g].kit_statuses[j] == '+' || gen_hierarchy[g].kit_statuses[j] == ' ') {
                                other_idx = gen_hierarchy[g].kit_indices[j];
                                SnpNode* sn_pos = kits[other_idx].snps;
                                while (sn_pos) {
                                    if (strcmp(sn_pos->name, sn_neg->name) == 0 && sn_pos->status == SNP_POSITIVE) {
                                        found_pos = 1;
                                        strcpy(pos_kit_id, kits[other_idx].id);
                                        break;
                                    }
                                    sn_pos = sn_pos->next;
                                }
                            }
                            if (found_pos) break;
                        }
                        
                        if (found_pos) {
                            int is_true_conflict = 1; 
                            SnpTreeNode* snp_node = find_snp_by_name(sn_neg->name);
                            SnpTreeNode* term_neg = find_raw_terminal_snp(&kits[kit_idx]);
                            SnpTreeNode* term_pos = find_raw_terminal_snp(&kits[other_idx]);
                            
                            if (snp_node && term_neg && term_pos) {
                                SnpTreeNode* mrca = get_static_snp_mrca(term_neg, term_pos);
                                if (mrca) {
                                    if (!is_snp_ancestor(snp_node, mrca)) is_true_conflict = 0;
                                }
                            }

                            if (is_true_conflict) {
                                if (warning_count < MAX_WARNINGS) {
                                    fprintf(stderr, "WARNING: GENDATA conflict! Group '%s' has conflicting SNP statuses (e.g., Kit %s is NEGATIVE for '%s', while Kit %s is POSITIVE).\n",
                                        full_gen_name, kits[kit_idx].id, sn_neg->name, pos_kit_id);
                                    warning_count++;
                                } else if (warning_count == MAX_WARNINGS) {
                                    fprintf(stderr, "WARNING: Maximum warning limit (%d) reached. Further warnings suppressed.\n", MAX_WARNINGS);
                                    warning_count++;
                                    return;
                                }
                                group_snp_warned = 1; 
                            }
                        }
                    }
                    sn_neg = sn_neg->next;
                }
            }
        }

        // 2. STR Outliers
        if (gen_node) {
            for (int i = 0; i < gen_hierarchy[g].kit_count; i++) {
                if (gen_hierarchy[g].kit_statuses[i] == '+' || gen_hierarchy[g].kit_statuses[i] == ' ') {
                    int kit_idx = gen_hierarchy[g].kit_indices[i];
                    int dist = 0, valid = 0;
                    for (int m = 0; m < total_markers; m++) {
                        int kv = kits[kit_idx].str_values[m], nv = gen_node->local_modal[m];
                        if (kv > 0 && nv > 0) { valid++; if (kv != nv) dist++; }
                    }
                    if (valid > 0 && (float)dist / valid >= 0.15) {
                        if (warning_count < MAX_WARNINGS) {
                            fprintf(stderr, "WARNING: GENDATA outlier! Kit %s in group '%s' has an unusually high STR distance (%d mutations / %d markers) from the group modal.\n",
                                kits[kit_idx].id, full_gen_name, dist, valid);
                            warning_count++;
                        } else if (warning_count == MAX_WARNINGS) {
                            fprintf(stderr, "WARNING: Maximum warning limit (%d) reached. Further warnings suppressed.\n", MAX_WARNINGS);
                            warning_count++;
                            return;
                        }
                    }
                }
            }
        }

        // 3. Unjustified Exclusions (Consolidated to ONE per GEN node)
        if (gen_node) {
            int group_exclusion_warned = 0;
            for (int i = 0; i < gen_hierarchy[g].kit_count && !group_exclusion_warned; i++) {
                if (gen_hierarchy[g].kit_statuses[i] == '-') {
                    int kit_idx = gen_hierarchy[g].kit_indices[i];
                    int dist = 0, valid = 0;
                    for (int m = 0; m < total_markers; m++) {
                        int kv = kits[kit_idx].str_values[m], nv = gen_node->local_modal[m];
                        if (kv > 0 && nv > 0) { valid++; if (kv != nv) dist++; }
                    }
                    if (valid > 0 && dist < 1) {
                        if (warning_count < MAX_WARNINGS) {
                            fprintf(stderr, "WARNING: GENDATA exclusion! Group '%s' has unjustified exclusions (e.g., Kit %s is explicitly excluded (-) but is genetically extremely close with STR distance %d).\n",
                                full_gen_name, kits[kit_idx].id, dist);
                            warning_count++;
                        } else if (warning_count == MAX_WARNINGS) {
                            fprintf(stderr, "WARNING: Maximum warning limit (%d) reached. Further warnings suppressed.\n", MAX_WARNINGS);
                            warning_count++;
                            return;
                        }
                        group_exclusion_warned = 1; // Stop checking exclusions for this GEN group
                    }
                }
            }
        }
    }
}

void report_tree_topology_conflicts() {
    int warning_count = 0;
    const int MAX_WARNINGS = 10;

    for (int k = 0; k < kit_count; k++) {
        if (kits[k].id[0] == '\0') continue;

        TreeNode* current_node = NULL;
        int current_idx = -1;
        
        // Find where the kit actually ended up in the final tree
        find_kit_in_active_tree(root_node, k, &current_node, &current_idx);
        if (!current_node) continue;

        // Trace up the tree from the kit's assigned node to the root
        TreeNode* ancestor = current_node;
        while (ancestor != NULL && ancestor->type != NODE_ROOT) {
            
            // If the kit is explicitly negative for this ancestor node, we have a structural conflict
            if (is_kit_negative_for_specific_node(k, ancestor)) {
                if (warning_count < MAX_WARNINGS) {
                    fprintf(stderr, "WARNING: Topology conflict! Kit '%s' is placed under node '%s', but explicitly tested NEGATIVE for it.\n",
                            kits[k].id, ancestor->name);
                    warning_count++;
                } else if (warning_count == MAX_WARNINGS) {
                    fprintf(stderr, "WARNING: Maximum topology warnings reached. Further warnings suppressed.\n");
                    warning_count++;
                }
                
                // Break out of the while loop for this kit so we don't spam multiple warnings 
                // for the same kit if the whole upper branch is wrong
                break; 
            }
            ancestor = ancestor->parent;
        }
    }
}

/* --- 8. Main Reassignment Loop --- */

void find_all_assigned_kits(TreeNode* node, int assigned[]) {
    if (!node) return;
    for (int i = 0; i < node->kit_count; i++) {
        assigned[node->kit_indices[i]] = 1;
    }
    TreeNode* c = node->first_child;
    while (c) {
        find_all_assigned_kits(c, assigned);
        c = c->next_sibling;
    }
}

void sweep_orphaned_kits(TreeNode* root) {
    if (!root) return;
    int assigned[MAX_KITS] = {0};
    find_all_assigned_kits(root, assigned);
    
    for (int i = 0; i < MAX_KITS; i++) {
        if (kits[i].id[0] == '\0') break; 
        if (!assigned[i]) {
            if (root->kit_count < MAX_KITS) {
                root->kit_indices[root->kit_count++] = i;
                assigned[i] = 1; 
            }
        }
    }
}

void compress_linear_str_chains(TreeNode* node) {
    if (!node) return;
    
    // Process bottom-up so deep chains zip up entirely in one pass
    TreeNode* c = node->first_child;
    while (c) {
        compress_linear_str_chains(c);
        c = c->next_sibling;
    }
    
    // If this is an empty STR node...
    if ((node->type == NODE_STR || node->type == NODE_STR_BRANCH) && node->kit_count == 0) {
        int child_count = 0;
        TreeNode* only_child = NULL;
        TreeNode* iter = node->first_child;
        while (iter) {
            child_count++;
            only_child = iter;
            iter = iter->next_sibling;
        }
        
        // ...and it has exactly ONE child that is also an STR node, stack them!
        if (child_count == 1 && (only_child->type == NODE_STR || only_child->type == NODE_STR_BRANCH)) {
            
            // 1. Stack the mutations
            for (int i = 0; i < only_child->mutation_count; i++) {
                if (node->mutation_count < 50) { // MAX_MUTATIONS_PER_NODE limit
                    node->mutations[node->mutation_count++] = only_child->mutations[i];
                }
            }
            
            // 2. Pull the kits up
            for (int i = 0; i < only_child->kit_count; i++) {
                node->kit_indices[node->kit_count++] = only_child->kit_indices[i];
            }
            
            // 3. Bypass the child and connect directly to grandchildren
            node->first_child = only_child->first_child;
            TreeNode* oc_child = node->first_child;
            while (oc_child) {
                oc_child->parent = node;
                oc_child = oc_child->next_sibling;
            }
            
            // 4. Inherit the deeper local modal
            for (int m = 0; m < marker_count; m++) {
                node->local_modal[m] = only_child->local_modal[m];
            }
        }
    }
}

void apply_predictive_clustering_recursive(TreeNode* node, int total_markers) {
    if (!node) return;

    // 1. Process this specific node ONLY if it meets the threshold
    int min_kits = relaxed_mode ? 2 : MIN_KITS_CLUSTER_NODE;
    if (node->kit_count >= min_kits) {
        branch_by_shared_mutations(node, total_markers);
        //cluster_inferred_str_nodes(node, total_markers);
    }

    // 2. Safely collect the original children before recursing.
    int orig_child_count = 0;
    TreeNode* orig_children[MAX_KITS];
    TreeNode* c = node->first_child;
    while (c) {
        orig_children[orig_child_count++] = c;
        c = c->next_sibling;
    }

    // 3. Recurse down the tree into the sub-branches (ALWAYS DO THIS)
    for (int i = 0; i < orig_child_count; i++) {
        apply_predictive_clustering_recursive(orig_children[i], total_markers);
    }
}

void check_kits_and_reassign(int total_markers) {
    int moved = 1, pass = 0;
    while (moved && pass++ < 10) {
        moved = 0;
        for (int k = 0; k < kit_count; k++) {
            if (kits[k].id[0] == '\0') continue;
            TreeNode* current_node = NULL; int current_idx = -1; find_kit_in_active_tree(root_node, k, &current_node, &current_idx);
            if (!current_node) continue;
            int is_curr_ok = is_acceptable_target(k, current_node);
            int current_dist = 0;
            for (int m = 0; m < total_markers; m++) { int kv = kits[k].str_values[m], nv = current_node->local_modal[m]; if (kv > 0 && nv > 0 && kv != nv) current_dist++; }
            int min_dist = is_curr_ok ? current_dist : 999999; TreeNode* best_node = is_curr_ok ? current_node : NULL;
            evaluate_reassignment_targets_recursive(root_node, k, total_markers, current_node, &best_node, &min_dist);
            if (best_node && best_node != current_node) {
                for (int j = current_idx; j < current_node->kit_count - 1; j++) current_node->kit_indices[j] = current_node->kit_indices[j+1];
                current_node->kit_count--; if (best_node->kit_count < MAX_KITS) best_node->kit_indices[best_node->kit_count++] = k;
                moved = 1;
            }
        }
        if (moved) { initialize_modals_top_down(root_node, modal_values, total_markers); refine_modals_bottom_up(root_node, total_markers); }
    }
    
    // 1. Run the strict biological scrubbers
    collapse_empty_str_nodes(root_node);
    impute_missing_snp_dates_proportional(root_node, total_markers); 
    calculate_str_node_ages_proportional(root_node, total_markers);

    // 2. Architecturally sound sweep: Rescue lost kits to the ROOT
    sweep_orphaned_kits(root_node);

    // 3. Establish mathematical baseline for the newly populated ROOT
    compute_missing_global_modals(total_markers);
    initialize_modals_top_down(root_node, modal_values, total_markers);
    
    // 4. Safely cluster basal kits into branches
    branch_by_shared_mutations(root_node, total_markers);
    //cluster_inferred_str_nodes(root_node, total_markers);
    apply_predictive_clustering_recursive(root_node, total_markers);
    compress_linear_str_chains(root_node);

    // 5. Force the modals to calculate 
    compute_missing_global_modals(total_markers);
    initialize_modals_top_down(root_node, modal_values, total_markers);
    refine_modals_bottom_up(root_node, total_markers);
    reevaluate_mutation_labels(root_node, total_markers);
    
    // 6. Final reports
    report_gendata_inconsistencies(total_markers);
    report_tree_topology_conflicts();
}

int calculate_str_distance(int strA[], int strB[], int total_markers) {
    int distance = 0;
    for (int i = 0; i < total_markers; i++) {
        if (strA[i] != STR_MISSING && strB[i] != STR_MISSING && strA[i] > 0 && strB[i] > 0) {
            distance += abs(strA[i] - strB[i]);
        }
    }
    return distance;
}

void branch_by_shared_mutations(TreeNode* node, int total_markers) {
    if (node->kit_count < MIN_SHARED) return;
    for (int m = 0; m < total_markers; m++) {
        int modal = node->local_modal[m];
        if (modal == STR_MISSING || modal <= 0) continue;

        int counts[256] = {0};
        for (int i = 0; i < node->kit_count; i++) {
            int val = kits[node->kit_indices[i]].str_values[m];
            if (val > 0 && val != modal && val < 256) counts[val]++;
        }

        for (int v = 0; v < 256; v++) {
            if (counts[v] >= 2) {
                TreeNode* newNode = &tree_nodes[tree_node_count++];
                newNode->type = NODE_STR_BRANCH;
                sprintf(newNode->name, "Branch_M%d_V%d", m, v);
                newNode->parent = node;
                newNode->kit_count = 0;
                newNode->first_child = NULL;
                newNode->next_sibling = NULL;
                newNode->date = 0;
                
                // Explicitly log the mutation so the HTML/Text renderer prints it!
                newNode->mutation_count = 1;
                newNode->mutations[0].marker_index = m;
                newNode->mutations[0].old_val = modal;
                newNode->mutations[0].new_val = v;

                for (int j = 0; j < total_markers; j++) newNode->local_modal[j] = node->local_modal[j];
                newNode->local_modal[m] = v;

                for (int i = 0; i < node->kit_count; i++) {
                    int k_idx = node->kit_indices[i];
                    if (kits[k_idx].str_values[m] == v) {
                        newNode->kit_indices[newNode->kit_count++] = k_idx;
                        node->kit_indices[i] = node->kit_indices[node->kit_count - 1];
                        node->kit_count--;
                        i--;
                    }
                }
                
                add_child_to_node(node, newNode);
                branch_by_shared_mutations(newNode, total_markers);
            }
        }
    }
}

void cluster_inferred_str_nodes(TreeNode* node, int total_markers) {
    if (node->kit_count < MIN_KITS_CLUSTER_ROOT) return;

    int active_kits[MAX_KITS];
    int active_count = node->kit_count;
    for (int i = 0; i < node->kit_count; i++) active_kits[i] = node->kit_indices[i];

    if (node->kit_count <= 2) return;

    while (active_count > 1) {
        int best_i = -1, best_j = -1, min_distance = 999999;
        for (int i = 0; i < active_count; i++) {
            for (int j = i + 1; j < active_count; j++) {
                int dist = calculate_str_distance(kits[active_kits[i]].str_values, kits[active_kits[j]].str_values, total_markers);
                if (dist < min_distance) { min_distance = dist; best_i = i; best_j = j; }
            }
        }

        if (best_i == -1 || best_j == -1) break;

        TreeNode* inferred_node = &tree_nodes[tree_node_count++];
        inferred_node->type = NODE_STR_BRANCH;
        sprintf(inferred_node->name, "Inferred_%d", str_node_counter++);
        inferred_node->parent = node;
        inferred_node->kit_count = 0;
        inferred_node->first_child = NULL;
        inferred_node->next_sibling = NULL;
        inferred_node->mutation_count = 0;
        inferred_node->date = 0;

        int kit1 = active_kits[best_i], kit2 = active_kits[best_j];
        inferred_node->kit_indices[inferred_node->kit_count++] = kit1;
        inferred_node->kit_indices[inferred_node->kit_count++] = kit2;

        for (int j = 0; j < total_markers; j++) inferred_node->local_modal[j] = node->local_modal[j];

        add_child_to_node(node, inferred_node);

        active_kits[best_j] = active_kits[active_count - 1]; active_count--;
        active_kits[best_i] = active_kits[active_count - 1]; active_count--;

        for (int k = 0; k < node->kit_count; k++) {
            if (node->kit_indices[k] == kit1 || node->kit_indices[k] == kit2) {
                node->kit_indices[k] = node->kit_indices[node->kit_count - 1];
                node->kit_count--;
                k--;
            }
        }
    }
}

#define DEFAULT_BIRTH_YEAR 1950

// Helper to get minimum birth year from direct kits
int get_min_kit_birth_year(TreeNode* node) {
    int min_y = 999999;
    for (int i = 0; i < node->kit_count; i++) {
        int k_idx = node->kit_indices[i];
        
        // NEW: If the kit has no birth year (0), default to 1950
        int effective_year = (kits[k_idx].birth_year > 0) ? kits[k_idx].birth_year : DEFAULT_BIRTH_YEAR;
        
        if (effective_year < min_y) {
            min_y = effective_year;
        }
    }
    return min_y;
}

// Helper to get minimum date from child nodes
int get_min_child_date(TreeNode* node) {
    int min_y = 999999;
    TreeNode* child = node->first_child;
    while (child) {
        if (child->type == NODE_MERGED && child->gen_date > 0 && child->gen_date < min_y) {
            min_y = child->gen_date;
        } else if (child->date > 0 && child->date < min_y) {
            min_y = child->date;
        }
        child = child->next_sibling;
    }
    return min_y;
}

// Helper to find the historical date of the immediate parent
int get_parent_date(TreeNode* node) {
    if (!node || !node->parent) return 0;
    
    // Always prefer the genealogical anchor if it's a merged node
    if (node->parent->type == NODE_MERGED && node->parent->gen_date > 0) {
        return node->parent->gen_date;
    }
    // Fallback to the standard SNP/Calculated date
    if (node->parent->date > 0) {
        return node->parent->date;
    }
    return 0;
}

// PASS 1: Locate STR leaves (no children)
void pass1_leaf_str_nodes(TreeNode* node) {
    if (!node) return;
    
    if ((node->type == NODE_STR_BRANCH || node->type == NODE_STR) && node->date <= 0 && node->first_child == NULL) {
        int oldest_kit = get_min_kit_birth_year(node);
        
        if (oldest_kit != 999999) {
            int calc_date = oldest_kit;
            
            // 1. Calculate base Genetic Distance
            if (node->mutation_count > 0) {
                calc_date = oldest_kit - (node->mutation_count * STR_GEN_LENGTH);
            }
            
            // 2. Apply Bounds Checking
            int parent_date = get_parent_date(node);
            
            if (parent_date > 0 && calc_date < parent_date) {
                calc_date = parent_date; // Floor: Cannot be older than the parent
            }
            if (calc_date > oldest_kit) {
                calc_date = oldest_kit;  // Ceiling: Cannot be younger than the kits
            }
            
            node->date = calc_date;
//            printf("MATH: Estimated Leaf %s -> %d (Anchor: %d, GD: %d, Parent Bound: %d)\n", 
//                   node->name, node->date, oldest_kit, node->mutation_count, parent_date);
        }
    }

    TreeNode* c = node->first_child;
    while (c) {
        pass1_leaf_str_nodes(c);
        c = c->next_sibling;
    }
}

// PASS 2: Locate STR parents (has children)
void pass2_parent_str_nodes(TreeNode* node) {
    if (!node) return;

    // Post-order traversal (bottom-up)
    TreeNode* c = node->first_child;
    while (c) {
        pass2_parent_str_nodes(c);
        c = c->next_sibling;
    }

    if ((node->type == NODE_STR_BRANCH || node->type == NODE_STR) && node->date <= 0 && node->first_child != NULL) {
        int oldest_overall = 999999;
        
        int oldest_kit = get_min_kit_birth_year(node);
        if (oldest_kit < oldest_overall) oldest_overall = oldest_kit;

        int oldest_child = get_min_child_date(node);
        if (oldest_child < oldest_overall) oldest_overall = oldest_child;

        if (oldest_overall != 999999) {
            int calc_date = oldest_overall;
            
            // 1. Calculate base Genetic Distance
            if (node->mutation_count > 0) {
                calc_date = oldest_overall - (node->mutation_count * STR_GEN_LENGTH);
            }
            
            // 2. Apply Bounds Checking
            int parent_date = get_parent_date(node);
            
            if (parent_date > 0 && calc_date < parent_date) {
                calc_date = parent_date; // Floor: Cannot be older than the parent
            }
            if (calc_date > oldest_overall) {
                calc_date = oldest_overall; // Ceiling: Cannot be younger than its descendants
            }

            node->date = calc_date;
//            printf("MATH: Estimated Parent %s -> %d (Anchor: %d, GD: %d, Parent Bound: %d)\n", 
//                   node->name, node->date, oldest_overall, node->mutation_count, parent_date);
        }
    }
}

void estimate_str_node_ages(TreeNode* root) {
    pass1_leaf_str_nodes(root);
    pass2_parent_str_nodes(root);
}

