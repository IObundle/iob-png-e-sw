PNG_E_DIR:=../..

CC := gcc
CFLAGS := -g -pedantic --std=gnu99

CSRCS := png_encoder.c lodepng.c
TRGT := png_encoder
CHDRS :=lodepng.h

# Image parameters
RAW_NAME ?= x_pattern.raw
PNG_NAME ?= x_pattern.png
IMG_WIDTH ?= 255
IMG_HEIGHT ?= 255
IMG_COLOR_TYPE ?= 0
IMG_BIT_DEPTH ?= 8
WINDOW_SIZE ?= 1024
IDAT_BUFFER_SIZE ?= 2048


all: clean run

$(TRGT): $(CSRCS) $(CHDRS)
	$(CC) $(CFLAGS) $(DEFINE) $(CSRCS) -o $@ -lm -lgcc -lc

run: $(TRGT)
	./$(TRGT) $(RAW_NAME) $(PNG_NAME) $(IMG_WIDTH) $(IMG_HEIGHT) \
	$(IMG_COLOR_TYPE) $(IMG_BIT_DEPTH) $(WINDOW_SIZE) $(IDAT_BUFFER_SIZE)

clean-sw:
	@rm -f *# *~ $(TRGT)
	
clean: clean-sw
	@rm -f *.png

.PHONY: all run clean
