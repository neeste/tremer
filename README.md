# TREMER

TREMER is a phylogenetic tree program designed to analyze genetic and genealogical data (such as Y-DNA STR and SNP data) and construct hierarchical trees detailing the relationships among various kits.

## Features
- **Phylogenetic Tree Construction:** Builds robust biological trees by evaluating shared STR mutations and established SNP hierarchies.
- **Genealogical Consistency:** Strictly enforces genealogical boundaries to ensure that predictive clustering does not pull kits away from their explicitly documented paper-trail ancestors.
- **Predictive Clustering:** Automatically infers unrecorded ancestral branches when multiple descendant lines share unique combinations of genetic mutations.
- **Modals and Distances:** Calculates ancestral modal haplotypes top-down and refines them bottom-up, reassigning orphans based on minimal genetic distance without violating known genealogy.

## Building the Project

The core engine is written in C for performance. You can build the main executable using `make`:

```bash
make
```

This will compile the source files (including `main.c`, `parse.c`, `tree_logic.c`, `render_svg.c`, and `render_text.c`) and produce the `tremer_rewrite` executable.

## Usage

You can run the program against a specific input file or a project directory. 

```bash
./tremer_rewrite <input_file_or_directory>
```

For example, to process a project folder:
```bash
./tremer_rewrite Neely_project
```

Or a specific input file:
```bash
./tremer_rewrite Neely_project/Neely_s2.txt
```

### Outputs
The program generates multiple outputs, including:
- `trees.txt`: A plain-text representation of the phylogenetic tree.
- `*.svg`: Scalable Vector Graphic visualizations of the trees.
- `*.html`: HTML wrappers for viewing the generated SVGs in a browser.

## Project Structure
- `main.c`: The core execution pipeline.
- `tree_logic.c`: The heavy-lifting algorithms for branching, distance clustering, rule evaluation, and modal computation.
- `parse.c`: Input parsers for STR/SNP genetic data and genealogical constraints.
- `render_text.c` / `render_svg.c`: Output generation modules.
- `wasm/`: Contains WebAssembly ports and JavaScript wrappers for running the engine in web environments.
