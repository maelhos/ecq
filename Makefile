# Default is x86_64
.PHONY: default
default: build

c_source_files := $(shell find src -name *.c)
c_object_files := $(patsubst src/%.c, build/%.o, $(c_source_files))

C_INCLUDE_PATH := include

OUT_NAME := ecq

LIBS := -lm -lgmp -L /usr/local/lib -l:libflint.so
OPTS := -march=native -Ofast #-g

$(c_object_files): build/%.o : src/%.c
	@mkdir -p $(dir $@) && \
	gcc $(OPTS) -c -Wall -I $(C_INCLUDE_PATH) $(patsubst build/%.o, src/%.c, $@) -o $@ $(LIBS)

.PHONY: build
build: $(c_object_files)
	@gcc $(OPTS) -Wall -o $(OUT_NAME) $(c_object_files) $(LIBS)

clean:
	@rm -rfv build/
	@mkdir build
	@echo "CLEAN"