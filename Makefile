# Default is x86_64
.PHONY: default
default: build

c_source_files := $(shell find src -name *.c)
c_object_files := $(patsubst src/%.c, build/%.o, $(c_source_files))
c_debug_object_files := $(patsubst src/%.c, build-debug/%.o, $(c_source_files))

C_INCLUDE_PATH := include

OUT_NAME := ecq
DEBUG_OUT_NAME := ecq-debug

LIBS := -lm -lgmp -L /usr/local/lib -l:libflint.so
OPTS := -march=native -O3 -g
DEBUG_OPTS := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined

$(c_object_files): build/%.o : src/%.c
	@mkdir -p $(dir $@) && \
	gcc $(OPTS) -c -Wall -I $(C_INCLUDE_PATH) $(patsubst build/%.o, src/%.c, $@) -o $@ $(LIBS)

.PHONY: build
build: $(c_object_files)
	@gcc $(OPTS) -Wall -o $(OUT_NAME) $(c_object_files) $(LIBS)

$(c_debug_object_files): build-debug/%.o : src/%.c
	@mkdir -p $(dir $@) && \
	gcc $(DEBUG_OPTS) -c -Wall -I $(C_INCLUDE_PATH) $(patsubst build-debug/%.o, src/%.c, $@) -o $@ $(LIBS)

.PHONY: debug
debug: $(c_debug_object_files)
	@gcc $(DEBUG_OPTS) -Wall -o $(DEBUG_OUT_NAME) $(c_debug_object_files) $(LIBS)

TEST_INPUTS := $(shell find tests/curves -type f)

.PHONY: test
test: build
	@for f in $(TEST_INPUTS); do \
		printf '%-24s ' "$$f"; \
		if [ ! -f tests/out/$$(basename $$f).out ]; then echo "NO out (run 'make record')"; continue; fi; \
		./$(OUT_NAME) < $$f 2>&1 | sed -n '1,/Starting main descent/p' > /tmp/ecq-$$(basename $$f).out; \
		if diff -q tests/out/$$(basename $$f).out /tmp/ecq-$$(basename $$f).out > /dev/null 2>&1; \
			then echo "OK"; \
			else echo "FAIL"; diff -u tests/out/$$(basename $$f).out /tmp/ecq-$$(basename $$f).out | head -30; fi; \
	done

.PHONY: record
record: build
	@mkdir -p tests/out
	@for f in $(TEST_INPUTS); do \
		echo "recording $$f"; \
		./$(OUT_NAME) < $$f 2>&1 | sed -n '1,/Starting main descent/p' > tests/out/$$(basename $$f).out; \
	done

clean:
	@rm -rfv ecq ecq-debug
	@rm -rfv build/ build-debug/
	@mkdir build
	@echo "CLEAN"
