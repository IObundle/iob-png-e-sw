#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lodepng.h"

unsigned IDAT_buffer_size = 32;

int loadFile(unsigned char** buffer, char* filename) //designed for loading files from hard disk in an stream
{
  
   // opening the file in read mode
   FILE* fp = fopen(filename, "r");

   // checking if the file exist or not
   if (fp == NULL) {
      printf("File %s not found\n", filename);
      return -1;
   }

   // calculating the size of the file
   fseek(fp, 0L, SEEK_END);
   int file_size = ftell(fp);
   fseek(fp, 0, SEEK_SET);  /* same as rewind(f); */

   //allocate space for file if file pointer not initialized
   if((*buffer) == NULL) {
      (*buffer) = (unsigned char *) malloc(file_size);
      if ((*buffer) == NULL) {
         printf("Error: malloc failed");
         fclose(fp);
         return -1;
      }
   }

   //write file to memory
   size_t read_n = fread((*buffer), 1, file_size, fp);
   if (read_n != file_size) {
      printf("Error: could not read entire file\n");
      fclose(fp);
      return -1;
   }
  
   // closing the file
   fclose(fp);
  
   return file_size;
}

int main (int argc, char* argv[]) {
   unsigned error;

   //printf("\nWelcome, this program will generate files with the PNG image.\n");

   if(argc != 9) {
      printf("Not enough, or too many arguments\n");
      printf("Usage is: ./<prog_name> <raw_name> <png_name> <img_width> <img_height> <img_color_type> <img_bit_depth> <window_size> <IDAT_buffer_size>\n");
      return 1;
   }
   
   char* raw_name = argv[1];
   char* png_name = argv[2];
   unsigned img_width = (unsigned) atoi(argv[3]);
   unsigned img_height = (unsigned) atoi(argv[4]);
   LodePNGColorType img_color_type = (LodePNGColorType) atoi(argv[5]);
   unsigned img_bit_depth = (unsigned) atoi(argv[6]);
   unsigned windowsize = (unsigned) atoi(argv[7]);
   IDAT_buffer_size = (unsigned) atoi(argv[8]);
   unsigned char* raw_img = NULL;
   
   if(loadFile(&raw_img, raw_name) == -1) {
      printf("Couldn't load the raw image file\n");
      return -1;
   }
   
   error = lodepng_encode_file(png_name, raw_img, img_width, img_height, 
                                 img_color_type, img_bit_depth, windowsize);
   /*if there's an error, display it*/
   if(error) printf("error %u: %s\n", error, lodepng_error_text(error));

   //printf("Finished .png image encoding\n");

	return 0;
}
