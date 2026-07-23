#ifndef TREMER_H
#define TREMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_KITS 1000           
#define MAX_MARKERS 1000        
#define MAX_SNPS 3000           
#define MAX_GEN_GROUPS 1000     
#define MAX_TREE_NODES 5000     
#define MAX_STRING_LEN 128
#define MAX_NODE_NAME_LEN 2048
#define MAX_STR_STACK 8 
#define STR_MISSING 0
#define STR_GEN_LENGTH 34 // Average years per generation (adjustable!)

typedef enum {
    SNP_NEGATIVE = -1,
    SNP_UNTESTED = 0,
    SNP_POSITIVE = 1
} SnpStatus;

typedef enum {
    NODE_ROOT,
    NODE_SNP,
    NODE_GEN,
    NODE_MERGED,
    NODE_STR_BRANCH, 
    NODE_STR,        
    NODE_DISTANCE 
} NodeType;

typedef struct SnpNode {
    char name[MAX_STRING_LEN];
    SnpStatus status;
    struct SnpNode* next;
} SnpNode;

typedef struct {
    char id[MAX_NODE_NAME_LEN]; 
    char group_label[MAX_STRING_LEN];
    int str_values[MAX_MARKERS];
    SnpNode* snps;              
    char gen_name[MAX_STRING_LEN]; 
    int gen_date;                  
    int birth_year;                  
} Kit;

typedef struct SnpTreeNode {
    char name[MAX_STRING_LEN];
    int mrca_date;
    char parent_name_str[MAX_STRING_LEN]; 
    struct SnpTreeNode* parent; 
} SnpTreeNode;

typedef struct {
    char ancestor_name[MAX_STRING_LEN];
    int date;
    int kit_indices[MAX_KITS];
    char kit_statuses[MAX_KITS];
    int kit_count;
} GenGroup;

typedef struct {
    int marker_index;
    int old_val;
    int new_val;
} StrMutation;

typedef struct TreeNode {
    NodeType type;
    char name[MAX_NODE_NAME_LEN]; 
    int date;          
    char gen_name[MAX_STRING_LEN]; 
    int gen_date;                  
    float x;
    int y;
    struct TreeNode* parent;
    struct TreeNode* first_child;
    struct TreeNode* next_sibling;
    int defined_kits[MAX_KITS];
    int defined_kit_count;
    int kit_indices[MAX_KITS];
    int kit_count;
    int local_modal[MAX_MARKERS]; 
    
    StrMutation mutations[20];
    int mutation_count;
} TreeNode;

extern Kit* kits;
extern int kit_count;
extern int modal_values[MAX_MARKERS];
extern int marker_count;
extern char (*marker_names)[MAX_STRING_LEN]; 
extern SnpTreeNode* snp_hierarchy;
extern int snp_hierarchy_count;
extern GenGroup* gen_hierarchy;
extern int gen_hierarchy_count;
extern TreeNode* tree_nodes; 
extern int tree_node_count;
extern TreeNode* root_node;
extern int str_node_counter; 
extern char project_name[MAX_STRING_LEN];
extern int global_group_num;
extern char root_haplogroup_letter;
extern int n_lab;
extern const char *str_lab[838];

int calculate_str_distance(int strA[], int strB[], int total_markers);
void branch_by_shared_mutations(TreeNode* node, int total_markers);
void cluster_inferred_str_nodes(TreeNode* node, int total_markers);

int peek_next_char(FILE* file);
int get_next_token(FILE* file, char* token);
int get_line_rest(FILE* file, char* line, int max_size);
int find_kit(const char* id);
int get_or_create_kit(const char* id);
int is_kit_id(const char* token);
void clean_snp_name(char* snp_name);
void parse_markers(FILE* file);
void parse_modal(FILE* file);
void parse_snptree(FILE* file);
void parse_gendata(FILE* file);
void parse_strdata(FILE* file);
void parse_snpdata(FILE* file);
void parse_kitdata(FILE* file);
void parse_groups(FILE* file);

SnpTreeNode* find_snp_by_name(const char* name);
int get_snp_depth(SnpTreeNode* snp);
SnpTreeNode* find_raw_terminal_snp(Kit* kit);
int is_snp_ancestor(SnpTreeNode* ancestor, SnpTreeNode* descendant);
SnpTreeNode* get_gen_group_mrca_snp(int gen_idx);
SnpTreeNode* find_imputed_terminal_snp(Kit* kit, int kit_idx);
int get_node_priority(NodeType type);
void add_child_to_node(TreeNode* parent, TreeNode* child);
int count_intersection(int* setA, int countA, int* setB, int countB);
void collect_subtree_kits(TreeNode* node, int* collected_indices, int* count);

void post_process_age_correction(TreeNode* node, int parent_date);
void post_process_merges(TreeNode* node);
void estimate_str_node_ages(TreeNode* node);
void apply_top_down_str_ages(TreeNode* node, int parent_date);

void impute_snp_statuses();
void build_skeleton_and_bucket_kits();
void compute_missing_global_modals(int total_markers);
void estimate_node_modal(TreeNode* node, int total_markers);
void initialize_modals_top_down(TreeNode* node, int* parent_modal, int total_markers);
void refine_modals_bottom_up(TreeNode* node, int total_markers);
void group_shared_str_mutations(TreeNode* node, int total_markers);
void apply_rule_of_two(TreeNode* node, int total_markers);
TreeNode* find_node_with_mutation(TreeNode* n, int m_idx, int old_val, int new_val, TreeNode* exclude);
void find_closest_node_recursive(TreeNode* node, int kit_idx, int total_markers, TreeNode** best_node, int* min_dist);
void reassign_root_orphans(int total_markers);
void check_kits_and_reassign(int total_markers); 
void post_process_parsimony_scrub(TreeNode* node, int total_markers);

void collapse_empty_str_nodes(TreeNode* node);
void sweep_orphaned_kits(TreeNode* root);

void finalize_str_node_numbering(TreeNode* root);
void reevaluate_mutation_labels(TreeNode* node, int total_markers);

void prepare_tree_for_svg(TreeNode* node);
int get_max_depth(TreeNode* node, int depth);
float calculate_coords(TreeNode* node, int depth, float* leaf_x);
float get_box_height(TreeNode* node);
void compute_level_heights(TreeNode* node, int depth, float max_box_h[], float max_mut_h[]);
void draw_svg_edges(FILE* f, TreeNode* node);
void draw_svg_nodes(FILE* f, TreeNode* node);
void generate_svg_output(const char* filename, TreeNode* root);
void generate_modal_output(const char* filename, TreeNode* root);

int print_indent(FILE* out_file, int depth);
void print_item_wrapped(FILE* out_file, const char* text, int* current_col, int wrap_indent_col, int max_width);
void generate_text_output(FILE* out_file, TreeNode* node, int depth);
void generate_html_table_output(const char* filename, TreeNode* root);
void generate_html_filtered_table_output(const char* filename, TreeNode* root);
void generate_html_gen_table_output(const char* filename, TreeNode* root);

#endif
