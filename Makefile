BUILD_DIR = build
INC_DIR = include
SRC_DIR = src
OUT_DIR = output

all: build_all

build_all: $(BUILD_DIR)
	cd $(BUILD_DIR); cmake $(CURDIR); make;

$(BUILD_DIR):
	mkdir -p $@; mkdir -p $@/tmp;

clean:
	rm -rf $(BUILD_DIR) $(OUT_DIR)/*

.PHONY: build_all
