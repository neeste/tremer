# TREMER

TREMER is a blazing-fast, privacy-first phylogenetic tree program designed to analyze genetic and genealogical data (such as Y-DNA STR and SNP data) and construct hierarchical trees detailing the relationships among various kits. 

TREMER runs locally on your machine or directly in your browser via WebAssembly, ensuring your genetic data remains entirely private and is never uploaded to third-party clustering servers.

## Features
- **Strict Biological Skeletons:** Builds robust biological trees by evaluating shared STR mutations against established SNP hierarchies, preventing predictive clustering from pulling kits away from known paper-trail ancestors.
- **Heuristic Discovery Mode (Relaxed):** Optionally bypass strict paper-trail locks to cluster kits based on overall STR similarity, bridging the gap with visual "eyeball test" clustering tools like SAPP.
- **Automated Data Pipelines:** Includes utilities for merging raw FTDNA CSV exports and automatically fetching live SNP age estimates directly from FTDNA Discover.
- **WebAssembly UI:** An interactive web interface to generate and view trees directly in your browser without installing command-line tools.

## Getting Started (Web Interface)
The easiest way to use TREMER is via the included web interface.
1. Open `admin/index.html` in your web browser.
2. Upload your project `.txt` file.
3. Check the "Heuristic Discovery Mode" box if you want relaxed SAPP-style clustering.
4. Click "Generate Tree" to instantly calculate and visualize your phylogenetic tree!

## Getting Started (Command Line)
The core engine is written in C for maximum performance.

### 1. Compile the Project
```bash
make
```
This compiles the `tremer_rewrite` executable and the utilities in the `tools/` directory (`strmerge` and `strdata`).

### 2. Generate Trees
You can run the program against a specific input file or a project directory. 
```bash
./tremer_rewrite Neely_project
```
To run in **Relaxed Mode** (SAPP-style heuristic clustering):
```bash
./tremer_rewrite -relaxed Neely_project
```

### 3. Streamlined Workflow
For project maintainers, TREMER includes an automated pipeline:
- **`make build`**: Automatically merges new FTDNA STR `.csv` exports placed in your `../YSTR` directory and securely injects them into your project files without overwriting manual edits.
- **`make update`**: Runs the `tools/update_tremer_ages.py` script to fetch live SNP age updates from FTDNA for all SNPs listed in your project files. The script is protected against race conditions and will safely preserve any manual edits you make while it is running.

## Anonymous Sample Data
To see exactly how TREMER formats its inputs and handles multiple CSVs, check out the `samples/` directory! It contains mock 111-STR exports, merged STR data, and a full sample project input.

### Outputs
TREMER generates multiple outputs:
- `trees.txt`: A plain-text representation of the phylogenetic tree.
- `*.svg`: Scalable Vector Graphic visualizations of the trees.
- `*.html`: HTML wrappers for viewing the generated SVGs in a browser.
