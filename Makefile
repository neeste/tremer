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

build : 
	make -C Neely_build build 

update : 
	make -C Neely_build update

em :
	@echo source emsdk/emsdk_env.sh

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

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
	rm -f $(TARGET) tree*.txt *_*.{svg,png,html} tremer.{js,wasm}
	rm -f *.json

