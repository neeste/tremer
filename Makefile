# Compiler and flags
CC = gcc
CFLAGS = -Wall -O3

# Executable and source file names
TARGET = tremer_rewrite
SRC = main.c parse.c tree_logic.c render_svg.c render_text.c

# Automatically gather input files that start with "Neely_s" and end with ".txt"
Neely_INPUTS = $(wildcard Neely_project/Neely_s*.txt)

# Phony targets to prevent conflicts with files 
.PHONY: run Neely wasm clean update em

# Default target:

run : Neely 

YSTR = ../YSTR
PROJ = Neely_project

build: strmerge strdata
	@echo "Merging STRs from $(YSTR)..."
	@for dir in $(YSTR)/group_*; do \
		if [ -d "$$dir" ]; then \
			grp=$$(basename "$$dir" | cut -d'_' -f2); \
			if ls $$dir/*.csv >/dev/null 2>&1; then \
				./strmerge -o strdata$${grp}.txt $$dir/*.csv; \
			fi \
		fi \
	done
	@if [ -f strdata1.txt ] && [ -f strdata2.txt ]; then \
		cat strdata1.txt strdata2.txt > strdata12.txt; \
	fi
	@echo "Updating project files with new STR data..."
	@for file in $(PROJ)/Neely_s*.txt; do \
		if [ -f "$$file" ]; then \
			./strdata -o tmp_$$(basename $$file) $$file; \
			mv tmp_$$(basename $$file) $$file; \
		fi \
	done
	@rm -f strdata*.txt tmp_*.txt
	@echo "STR Build complete!"

update:
	@echo "Updating SNP ages..."
	python3 tools/update_tremer_ages.py $(PROJ)/Neely_s*.txt

em :
	@echo source emsdk/emsdk_env.sh

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

strmerge: tools/strmerge.c
	$(CC) $(CFLAGS) -o strmerge tools/strmerge.c

strdata: tools/strdata.c
	$(CC) $(CFLAGS) -o strdata tools/strdata.c

# Compile for WebAssembly
wasm: $(SRC) run
	@echo "Compiling C code to WebAssembly..."
	emcc $(SRC) -o tremer.js -O3 \
	  -s EXPORTED_FUNCTIONS=_main,_malloc,_free \
	  -s EXPORTED_RUNTIME_METHODS=callMain,FS \
	  -s FORCE_FILESYSTEM=1 \
	  -s ALLOW_MEMORY_GROWTH=1 \
	  -s TOTAL_STACK=67108864 \
	  -s INITIAL_MEMORY=268435456 \
	  -s INVOKE_RUN=0
	@echo "WebAssembly compilation complete. tremer.js and tremer.wasm generated."
	mv tremer.js tremer.wasm wasm

# Compile and then run the executable on all input files

Neely: $(TARGET)
	@if [ -z "$(Neely_INPUTS)" ]; then \
		echo "No input files found."; \
	else \
		echo "Executing $(TARGET) on all Neely files..."; \
		for file in $(Neely_INPUTS); do \
			./$(TARGET) $$file ; \
		done ; \
		cat tree[1-9].txt > trees.txt ; \
		rm -f tree[1-9].txt Neely_*_*.html *.json ; \
		mv trees.txt  Neely_*.* Neely_project ; \
		echo "Neely strict processing complete."; \
	fi

relaxed: $(TARGET)
	@if [ -z "$(Neely_INPUTS)" ]; then \
		echo "No input files found."; \
	else \
		echo "Executing $(TARGET) -relaxed on all Neely files..."; \
		for file in $(Neely_INPUTS); do \
			./$(TARGET) -relaxed $$file ; \
		done ; \
		cat tree[1-9].txt > trees_relaxed.txt ; \
		rm -f tree[1-9].txt Neely_*_*.html *.json ; \
		mv trees_relaxed.txt  Neely_project ; \
		echo "Neely relaxed processing complete."; \
	fi

# Clean up the compiled executable and the generated output files
clean:
	rm -f $(TARGET) strmerge strdata tree*.txt *_*.{svg,png,html} tremer.{js,wasm}
	rm -f *.json strdata*.txt
