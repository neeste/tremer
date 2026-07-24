# Compiler and flags
CC = gcc
CFLAGS = -Wall -O3

# Executable and source file names
TARGET = tremer_rewrite
SRC = src/main.c src/parse.c src/tree_logic.c src/render_svg.c src/render_text.c

# Automatically gather input files that start with "Neely_s" and end with ".txt"
Neely_INPUTS = $(wildcard Neely_project/Neely_s*.txt)

# Phony targets to prevent conflicts with files 
.PHONY: run Neely wasm clean update em

# Default target:

run : Neely 

YSTR = ../YSTR
PROJ = Neely_project

# Deployment configuration
DEPLOY_HOST = tremer_deploy@bonkachen.com
DEPLOY_PATH = public_html/tremer

GROUPS := $(patsubst $(YSTR)/group_%,%,$(wildcard $(YSTR)/group_*))
$(foreach grp,$(GROUPS),$(eval strdata_out/strdata$(grp).txt: $(wildcard $(YSTR)/group_$(grp)/*.csv)))

build: strdata_out $(Neely_INPUTS)
	@echo "STR Build complete!"

strdata_out:
	@mkdir -p strdata_out

strdata_out/strdata12.txt: strdata_out/strdata1.txt strdata_out/strdata2.txt | strdata_out
	@echo "Creating strdata12.txt..."
	@cat strdata_out/strdata1.txt strdata_out/strdata2.txt > strdata_out/strdata12.txt

strdata_out/strdata%.txt: strmerge | strdata_out
	@if ls $(YSTR)/group_$*/*.csv >/dev/null 2>&1; then \
		echo "Merging STRs for group $*..."; \
		./strmerge -o $@ $(YSTR)/group_$*/*.csv; \
	else \
		touch $@; \
	fi

$(PROJ)/Neely_s%.txt: strdata_out/strdata%.txt strdata
	@echo "Updating $@ with new STR data..."
	@./strdata -o tmp_$(notdir $@) $@
	@mv tmp_$(notdir $@) $@

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
	rm -rf $(TARGET) strmerge strdata tree*.txt *_*.{svg,png,html} tremer.{js,wasm} strdata_out *.json

# Deploy to Site5 web host
deploy:
	@echo "Deploying to Site5..."
	rsync -avz admin Neely_project wasm $(DEPLOY_HOST):$(DEPLOY_PATH)/
	@echo "Deployment complete!"
