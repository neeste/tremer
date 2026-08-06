#include "tremer.h"

#define SVG_X_SPAC 176
#define SVG_MARGIN 50
#define BOX_W 128

static float level_y_coords[5000] = {0};

void find_assigned_kits(TreeNode* node, int assigned[]) {
    if (!node) return;
    for (int i = 0; i < node->kit_count; i++) {
        assigned[node->kit_indices[i]] = 1;
    }
    TreeNode* c = node->first_child;
    while (c) {
        find_assigned_kits(c, assigned);
        c = c->next_sibling;
    }
}

void prepare_tree_for_svg(TreeNode* node) {
    if (!node) return;
    
    if (node->kit_count > 0) {
        if (tree_node_count < MAX_TREE_NODES) {
            TreeNode* kit_node = &tree_nodes[tree_node_count++];
            kit_node->type = NODE_DISTANCE; 
            strcpy(kit_node->name, ""); 
            kit_node->date = 0;
            kit_node->kit_count = node->kit_count;
            for(int i = 0; i < node->kit_count; i++) kit_node->kit_indices[i] = node->kit_indices[i];
            for(int m = 0; m < MAX_MARKERS; m++) kit_node->local_modal[m] = node->local_modal[m];
            
            kit_node->first_child = NULL;
            kit_node->next_sibling = NULL;
            
            node->kit_count = 0;
            add_child_to_node(node, kit_node);
        }
    }
    
    TreeNode* c = node->first_child;
    while (c) {
        TreeNode* next = c->next_sibling;
        if (c->type != NODE_DISTANCE) prepare_tree_for_svg(c);
        c = next;
    }
}

int get_max_depth(TreeNode* node, int depth) {
    if (!node) return depth;
    int max = depth;
    TreeNode* c = node->first_child;
    while (c) {
        int d = get_max_depth(c, depth + 1);
        if (d > max) max = d;
        c = c->next_sibling;
    }
    return max;
}

float calculate_coords(TreeNode* node, int depth, float* leaf_x) {
    if (!node) return 0;
    node->y = depth; 
    if (!node->first_child) {
        node->x = *leaf_x;
        *leaf_x += 1.0;
    } else {
        TreeNode* c = node->first_child;
        float first_x = -1, last_x = -1;
        while (c) {
            float cx = calculate_coords(c, depth + 1, leaf_x);
            if (first_x < 0) first_x = cx;
            last_x = cx;
            c = c->next_sibling;
        }
        node->x = (first_x + last_x) / 2.0;
    }
    return node->x;
}

float get_box_height(TreeNode* node) {
    float box_h = 10;
    if (node->type == NODE_DISTANCE) {
        box_h += node->kit_count * 16;
    } else {
        if (node->type == NODE_ROOT) {
            box_h += 40; 
            if (root_haplogroup_letter != '\0') box_h += 18;
        } else if (node->type == NODE_MERGED && strlen(node->gen_name) > 0) {
            box_h += 36; // 2 lines (18px each) for SNP and GEN
        } else {
            box_h += 18; // 1 line for standard nodes
        }
    }
    if (box_h < 30) box_h = 30;
    return box_h;
}

void compute_level_heights(TreeNode* node, int depth, float max_box_h[], float max_mut_h[]) {
    if (!node) return;
    float bh = get_box_height(node);
    if (bh > max_box_h[depth]) max_box_h[depth] = bh;
    
    if (node->parent != NULL && node->type != NODE_DISTANCE) {
        int occ_count = node->mutation_count;
        if (occ_count > MAX_STR_STACK) occ_count = MAX_STR_STACK;
        
        float mh = 0;
        if (occ_count > 0) {
            mh = (occ_count * 12) - 5; 
            if (mh < 0) mh = 0;
        }
        
        if (mh > max_mut_h[depth]) max_mut_h[depth] = mh;
    }
    
    TreeNode* c = node->first_child;
    while (c) {
        compute_level_heights(c, depth + 1, max_box_h, max_mut_h);
        c = c->next_sibling;
    }
}

void draw_svg_edges(FILE* f, TreeNode* node) {
    if (!node) return;
    float box_h = get_box_height(node);

    TreeNode* c = node->first_child;
    while (c) {
        float px = node->x * SVG_X_SPAC + SVG_MARGIN + BOX_W / 2.0;
        float py = level_y_coords[node->y] + box_h; 
        float cx = c->x * SVG_X_SPAC + SVG_MARGIN + BOX_W / 2.0;
        float cy = level_y_coords[c->y]; 
        float mid_y = py + 20; 
        
        fprintf(f, "<path d=\"M %.1f %.1f L %.1f %.1f L %.1f %.1f L %.1f %.1f\" fill=\"none\" stroke=\"#8ab4f8\" stroke-width=\"3\"/>\n", 
                px, py, px, mid_y, cx, mid_y, cx, cy);
                
        draw_svg_edges(f, c);
        c = c->next_sibling;
    }
}

void draw_svg_nodes(FILE* f, TreeNode* node) {
    if (!node) return;
    
    float x = node->x * SVG_X_SPAC + SVG_MARGIN;
    float y = level_y_coords[node->y];
    float box_h = get_box_height(node);
    
    const char* fill = (node->type == NODE_DISTANCE) ? "#FFFFE0" : "#ADD8E6";
    const char* stroke = (node->type == NODE_DISTANCE) ? "#FFC0CB" : "#87CEEB";

    fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%d\" height=\"%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"2\" rx=\"5\" ry=\"5\"/>\n",
            x, y, BOX_W, box_h, fill, stroke);
            
    float text_y;
    
    if (node->type == NODE_DISTANCE) {
        text_y = y + 16;
        for (int i = 0; i < node->kit_count; i++) {
            int kit_idx = node->kit_indices[i];
            int valid_strs = 0;
            for(int m = 0; m < marker_count; m++) if(kits[kit_idx].str_values[m] > 0) valid_strs++;
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"11\" fill=\"#000\" text-anchor=\"start\">%s</text>\n", x + 10.0, text_y, kits[kit_idx].id);
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"11\" fill=\"#000\" text-anchor=\"end\">%d</text>\n", x + BOX_W - 10.0, text_y, valid_strs);
            text_y += 16;
        }
        
        if (node->kit_count > 0 && strlen(kits[node->kit_indices[0]].group_label) > 0) {
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"11\" fill=\"red\" text-anchor=\"start\">%s</text>\n", x + 10.0, y + box_h + 14.0, kits[node->kit_indices[0]].group_label);
        }
        
    } else {
        int printed_lines = 0;
        
        if (node->type == NODE_ROOT) {
            text_y = y + 20;
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"middle\">Group %d</text>\n", x + BOX_W/2.0, text_y, global_group_num);
            text_y += 18;
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" fill=\"#333\" text-anchor=\"middle\">%d kits</text>\n", x + BOX_W/2.0, text_y, kit_count);
            if (root_haplogroup_letter != '\0') {
                text_y += 18;
                fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" fill=\"#333\" text-anchor=\"middle\">Haplogroup %c</text>\n", x + BOX_W/2.0, text_y, root_haplogroup_letter);
            }
            printed_lines = 1; 
            
        } else if (node->type == NODE_MERGED && strlen(node->gen_name) > 0) {
            // --- NEW STRUCTURAL PARSING FOR MERGED NODES ---
            text_y = y + 20;
            
            // Print Top Line: SNP Name and Date
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"middle\">", x + BOX_W/2.0, text_y);
            if (node->date > 0) fprintf(f, "%s.%d", node->name, node->date);
            else fprintf(f, "%s", node->name);
            fprintf(f, "</text>\n");
            
            text_y += 18;
            printed_lines++;
            
            // Print Bottom Line: GEN Name and Date
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"middle\">", x + BOX_W/2.0, text_y);
            if (node->gen_date > 0) fprintf(f, "%s.%d", node->gen_name, node->gen_date);
            else fprintf(f, "%s", node->gen_name);
            fprintf(f, "</text>\n");
            
            text_y += 18;
            printed_lines++;
            
        } else {
            // --- STANDARD UNMERGED NODES ---
            text_y = y + 20;
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"middle\">", x + BOX_W/2.0, text_y);
            if (node->date > 0) {
                fprintf(f, "%s.%d", node->name, node->date);
            } else {
                fprintf(f, "%s", node->name);
            }
            fprintf(f, "</text>\n");
            text_y += 18;
            printed_lines++;
        }
        
        if (printed_lines == 0) {
             fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#000\" text-anchor=\"middle\">STR</text>\n", x + BOX_W/2.0, y + 20);
        }

        if (node->parent != NULL) {
            int mut_y = y - 10; 
            for (int i = 0; i < node->mutation_count && i < MAX_STR_STACK; i++) {
                int m_idx = node->mutations[i].marker_index;

                int occ = 0; 
                int sub_kits[MAX_KITS]; int sub_c = 0;
                collect_subtree_kits(node, sub_kits, &sub_c);
                for (int k = 0; k < sub_c; k++) {
                    if (kits[sub_kits[k]].str_values[m_idx] == node->mutations[i].new_val) occ++;
                }

                fprintf(f, "<text x=\"%.1f\" y=\"%d\" font-family=\"monospace\" font-size=\"10\" fill=\"#555\" text-anchor=\"start\">%s=%d-&gt;%d</text>\n", 
                        x + 10, mut_y, marker_names[m_idx], node->mutations[i].old_val, node->mutations[i].new_val);
                fprintf(f, "<text x=\"%.1f\" y=\"%d\" font-family=\"monospace\" font-size=\"10\" fill=\"#555\" text-anchor=\"end\">%d</text>\n", 
                        x + BOX_W - 10, mut_y, occ);
                mut_y -= 12; 
            }
        }

        // --- NEW KIT COUNT LOGIC ---
        // 1. Collect all kits attached to this node AND its downstream branches
        int total_kits_for_node = 0;
        int temp_sub_kits[MAX_KITS]; // Using your existing MAX_KITS definition
        collect_subtree_kits(node, temp_sub_kits, &total_kits_for_node);

        // 2. Print the total count below the box, right-justified
        if (total_kits_for_node > 0) {
            fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"11\" fill=\"#555\" text-anchor=\"end\">%d</text>\n", 
                    x + BOX_W - 10.0, y + box_h + 14.0, total_kits_for_node);
        }
        // ---------------------------
        
    }
    
    TreeNode* c = node->first_child;
    while (c) { draw_svg_nodes(f, c); c = c->next_sibling; }
}

void generate_svg_output(const char* filename, TreeNode* root) {
    FILE* f = fopen(filename, "w");
    if (!f) return;
    
    int assigned_map[MAX_KITS] = {0};
    find_assigned_kits(root, assigned_map);
    
    for (int i = 0; i < kit_count; i++) {
        if (assigned_map[i] == 0) {
            if (root->kit_count < MAX_KITS) {
                root->kit_indices[root->kit_count++] = i;
            }
        }
    }
    
    prepare_tree_for_svg(root);
    
    float next_leaf_x = 0;
    calculate_coords(root, 0, &next_leaf_x);
    
    int max_depth = get_max_depth(root, 0);
    float max_box_height[5000] = {0};
    float max_mut_height[5000] = {0};
    
    compute_level_heights(root, 0, max_box_height, max_mut_height);
    
    level_y_coords[0] = SVG_MARGIN + max_mut_height[0];
    for (int i = 1; i <= max_depth; i++) {
        level_y_coords[i] = level_y_coords[i-1] + max_box_height[i-1] + max_mut_height[i] + 40;
    }
    
    float actual_leaves = next_leaf_x > 0 ? next_leaf_x - 1.0 : 0;
    float svg_w = actual_leaves * SVG_X_SPAC + SVG_MARGIN * 2 + BOX_W;
    float svg_h = level_y_coords[max_depth] + max_box_height[max_depth] + 20.0 + SVG_MARGIN;
    
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char date_str[64];
    strftime(date_str, sizeof(date_str), "%d-%b-%Y", &tm);
    
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" style=\"background-color: white;\">\n", svg_w, svg_h);
    
    fprintf(f, "<text x=\"%.1f\" y=\"%.1f\" font-family=\"Arial\" font-size=\"12\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"end\">%s</text>\n", svg_w - SVG_MARGIN, (float)(SVG_MARGIN - 20), date_str);
    
    draw_svg_edges(f, root);
    draw_svg_nodes(f, root);
    fprintf(f, "</svg>\n");
    fclose(f);
}
