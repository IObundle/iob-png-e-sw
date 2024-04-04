# IOb-PNG-E-Software

Contains all the necessary software to encode a valid PNG image
without dependencies on any external libraries by using the LodePNG
encoder (more information about this open-source encoder can be found at the end of the [header](lodepng.h) file).

## Configuring the Input and Output files

To configure the names of the input and output files, the parameters on the [Makefile](Makefile) can be edited, or they can be passed through the
command line. The parameters are:
- RAW_NAME (default = x_pattern.raw): Sets the input raw image file name
- PNG_NAME (default = x_pattern.png): Sets the output .png file name
- IMG_WIDTH (default = 255): Sets the image width
- IMG_HEIGHT (default = 255): Sets the image height
- IMG_COLOR_TYPE (default = 0): Sets the image color type (0 = grayscale, 2 = RGB, 3 = indexed color, 4 = grayscale with alpha, 6 = RGBA)
- IMG_BIT_DEPTH (default = 8): Sets the image bit depth
- WINDOW_SIZE (default = 1024): Sets the window size for the LZ77 compression
- IDAT_BUFFER_SIZE (default = 2048): Sets the IDAT buffer size for the LZ77 compression

## Running and cleaning the encoder

To build and run the PNG encoder, run the command:
```
make run [RAW_NAME=<input raw file name>] [PNG_NAME=<output png file name>] [IMG_WIDTH=<image width>] [IMG_HEIGHT=<image height>] [IMG_COLOR_TYPE=<image color type>] [IMG_BIT_DEPTH=<image bit depth>] [WINDOW_SIZE=<window size>] [IDAT_BUFFER_SIZE=<IDAT buffer size>]
```
To use a fresh build and run the PNG encoder, run the command:
``` 
make all [RAW_NAME=<input raw file name>] [PNG_NAME=<output png file name>] [IMG_WIDTH=<image width>] [IMG_HEIGHT=<image height>] [IMG_COLOR_TYPE=<image color type>] [IMG_BIT_DEPTH=<image bit depth>] [WINDOW_SIZE=<window size>] [IDAT_BUFFER_SIZE=<IDAT buffer size>]
```
To clean the PNG decoder, run the command: 
```
make clean
```
