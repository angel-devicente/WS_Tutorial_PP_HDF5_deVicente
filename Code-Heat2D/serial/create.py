#!/usr/bin/env python3

# Simple function to create an ASCII .dat file out of a .png file
# Used to get an ASCII representation of the IAC logo:

from PIL import Image

def create_dat(pngfile,datfile):
   img = Image.open(pngfile).convert('L')  # convert image to 8-bit grayscale
   WIDTH, HEIGHT = img.size

   data = list(img.getdata()) # convert image data to a list of integers
   # convert that to 2D list (list of lists of integers)
   data = [data[offset:offset+WIDTH] for offset in range(0, WIDTH*HEIGHT, WIDTH)]

   # At this point the image's pixels are all in memory and can be accessed
   # individually using data[row][col].

   f = open(datfile,'w')
   f.write("# {0} {1} \n".format(WIDTH,HEIGHT)) 
   # For example:
   for row in data:
        f.write(' '.join('{:3}'.format(255-value) for value in row))
        f.write("\n")

   f.close()

        
