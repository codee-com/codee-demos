TARGET_DIR = bin
CLASS = B
FILE ?= cg.c
TARGET ?= cg

default: run

clean:
	rm -fr $(TARGET_DIR)
	rm -f common/*.o CG/*.o

build: clean
	mkdir -p $(TARGET_DIR)
	$(MAKE) -C CG CLASS=$(CLASS) FILE=$(FILE) TARGET=$(TARGET)

run: build
	./$(TARGET_DIR)/cg.$(CLASS).x
