/*
    SRTM2STL - Creates STL file from SRTM height data.
    Copyright (C) 2014  Thomas P. Sullivan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

// ************************************************************
// ************************************************************
// ** 
// ** SRTM to STL
// ** 
// ** Thomas P. Sullivan and Eric J. Fallon
// ** 
// ** Original Version 5/9/2014
// ** 
// ** Version 3.05 - June   17, 2014
// ** Version 3.06 - June   19, 2014
// ** Version 3.07 - June   21, 2014
// ** Version 3.08 - July   16, 2014
// ** Version 3.09 - July   18, 2014
// ** Version 3.10 - August  1, 2014
// ** Version 3.11 - August 11, 2014
// ** Version 3.12 - August 24, 2014
// **
// ** Converts SRTM .HGT files to ASCII or Binary STL files.
// ** 
// ************************************************************
// ************************************************************
/*
   Description:

   This program extracts data contained in SRTM 3 arc-second, 1 degree by 1 degree, 1200 x 1200 data point
   elevation files and constructs a raised-relief map as an STL file. The STL file can be converted to a
   format suitable for, and printed using, a 3D printer.  

   Source Data: 
   
   SRTM data can be found at http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/

   Using this program:
   
   Each 1 degree by 1 degree file contains 1201 by 1201 data points which represent 16 bit binary integer
   elevations in meters. After extraction from a ZIP format, each file is named using the format (no spaces):

   Latitude (N or S),  DD (degrees 00 to 89), Longitude (E or W), DDD (000 to 179)
   
   This name identifies the lower left-hand coordinate of the 1x1 grid. For example, Crater Lake is located
   about 42 56.5 N 122 6.38 W. To extract the data for Crater Lake you would need to use file N42W123.hgt.
   
   To extract the entire grid of data you would run this program like this:
   
   srtm2stl_3_12  N42W123.hgt craterlake.stl /S 0,0,1200,1200 /A 1.0 /B -1000 /F /T- /V
   
   Here is an explanation of the command line used above. 
   
   The data is N42W123.hgt which contains Crater Lake.
   
   The output file is named craterlake.stl. 
   
   The /S option allows us to extract whatever we want from within the grid. In the line above, we used the /S
   option to extract all of the grid of data. The format is, /S X1,Y1,X2,Y2. The /S option uses the Cartesian
   coordinate system with the origin at the lower left-hand part of the grid. So X1,Y1 is 0,0 which means start
   at the extreme lower left-hand coordinate and extract all of the data up to 1200,1200 which is the
   extreme upper right-hand coordinate of the grid. 
   
   The /A option is amplitude and allows you to scale the 3D surface by factor. In the example above we scaled
   at the default value of 1.0. 
   
   The /B option is bias and allows you to add or remove thickness to the 3D map. This is useful especially useful
   for sea level maps and mountainous terrain. It is also useful in shedding thickness when the scale factor is
   greater than 1.0. In the example above we are subtracting 1000 meters from the map thickness. 
   
   The /F option means that we want the program to fill in holes using the very simple (and somewhat limited) linear
   fill algorithm. This algorithm works well for small numbers of adjacent holes. It is not so good at fixing
   large holes in steep mountainous terrain. The source data files contain missing data for various reasons. You can
   read more about it in the report found here: 
   
   The /T- options means the STL output file will be in binary format. /T+ would be text (ASCII).
   
   The /V option is Verbose and forces extra information to be echoed to the command line.
   
   If the command line above is run and the output STL file is viewed using a STL file viewer (like MeshLab:
   http://meshlab.sourceforge.net) and oriented with north up and south down, crater lake will be in the extreme upper
   right-hand side of the map. To extract just the area around Crater lake you could compute the offsets within the grid
   by hand or use trial and error. You might arrive at the new /S option as /S 950,1050,1200,1200. This would do a reasonable
   job. You could also try the /C option instead. This would allow you to extract using minutes within the grid for Lat/Lon.
   The format for the coordinates is the same as /S i.e. X1,Y1,X2,Y2 but in this case, Latitude and longitude of the 
   lower left-hand and upper right-hand coordinates are expressed as minutes. No sign data is needed. In the case of
   Crater Lake, we could try /C 52,12,59.9,0.

*/


#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "geometry.h"
#include "stlwriter.h"

#define STL_HEADER_LENGTH 80
#define MAXINPUTLINESIZE        132
#define TRUE 1
#define FALSE 0

#define PIE ((double) 3.14159265358979)

#define _BIAS   ((long) 0)

#define ThreeMinuteFiles
#ifdef ThreeMinuteFiles

#define _ROWS 1201
#define _COLS 1201

#define _ROWSTART 0
#define _ROWEND 1200
#define _COLSTART 0
#define _COLEND 1200

#else

#define OneMinutesFiles

#define _ROWS 3601
#define _COLS 3601

#define _ROWSTART 0
#define _ROWEND 3600
#define _COLSTART 0
#define _COLEND 3600

#endif


// ************************
// SRTM2STL_ASCII
// 
// Convert SRTM data to STL
// ************************

int main(int argc, char *argv[])
{
   FILE *in = NULL;                        // input file
   FILE *out = NULL;                       // output file

   char SolidName[132] = {' '};            // Optional name for solid (text format STL files)
   int Verbose = 0;                        // Flag: Verbose
   int HoleErrors = FALSE;                 // Flag: True means there were hole errors
   int FixHoles = FALSE;                   // Flag: True to Fix any holes
   int Binary = TRUE;                      // TRUE for binary and FALSE for text
   int row,column;
   int facets = 0;                         // For counting the number of facets (triangles) created

   coordinate v1, v2, v3;                  // For point on triangle
   coordinate normal;                      // Used with function to get computed normal vector

   char *ptr;
   int ExtractionRow1 = _ROWSTART;
   int ExtractionColumn1 = _COLSTART;
   int ExtractionRow2 = _ROWEND;
   int ExtractionColumn2 = _COLEND;

   float AmpFactor = 1.0;                  // Amplification factor (multiply elevations by this)
   long Bias = 0;                          // Bias is the amount (positive or negative) to add or subtract from the map
   short grid[_ROWS][_COLS];               // Array to hold the elevations from 3 arc second files
   unsigned int i = 0;
   short stemp;
   int NumHoles = 0;
   int FixedHoles = 0;

   char LatDir,LonDir;
   double Latitude, Longitude;
        
   double xfact, yfact;
   double xSide, ySide;

   // For linear hole filling
   int _leftbank = -1;
   int _rightbank = -1;                            
   int _fill,_slope;

   int LonLowerLeft;
   int LatLowerLeft;
   int LonUpperRight;
   int LatUpperRight;

   float LatLowerLeftCoord;
   float LonLowerLeftCoord;
   float LatUpperRightCoord;
   float LonUpperRightCoord;

   // Default Solid Name for ASCII files
   strcpy(SolidName, "some_name\0");

   if (argc < 3)               // has command line enough words??
   {
      printf("\n");
      printf("  STRM2STL Copyright (C) 2014  Thomas P. Sullivan\n");
      printf("  This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.\n");
      printf("  This is free software, and you are welcome to redistribute it\n");
      printf("  under certain conditions; type `show c' for details.\n");

      printf ("\n  USAGE -  srtm2stl_ascii_2_xx input output /A|a x.x /B|b bias /F|f H|h+|- /S|s r,c,r,c /T|t+|- /V|v \n");
      printf ("\n  A|a - Amplification - Factor to multiply the heights by. This is");
      printf ("\n        applied before the bias is added/subtracted");
      printf ("\n  B|b - Bias adds (+) or removes (-) thickness. Must be an integer");
      printf ("\n  C|c - Coordinate subset - Latitude,Longitude,Latitude,Longitude (decimal minutes)");
      printf ("\n                            {  Lower Left    } {   Upper Right  }");
      printf ("\n  F|f - Fix holes (Linear Fill)");
      printf ("\n  N|n - Change the Solid Name is the STL ASCII files");
      printf ("\n  H|h - /H+ means print the array coordinates of hole (missing data). H- means don't (default)");
      printf ("\n  S|s - Subset - Start_Row,Start_Column,Stop_Row,Stop_Column");
      printf ("\n  T|t - T+ means output ASCII format and T- means output binary format");
      printf ("\n  V|v - Verbose\n\n\n");
      exit(1);
   }

   printf("\n");
   printf("  STRM2STL Copyright (C) 2014  Thomas P. Sullivan\n");
   printf("  This program comes with ABSOLUTELY NO WARRANTY; for details see GNU License file.\n");
   printf("  This is free software, and you are welcome to redistribute it\n");
   printf("  under certain conditions; for details see GNU License file.\n");
   printf("\n");

   in = fopen (argv[1], "rb");        // open input file
   if (in == NULL)
   {
      printf ("       Hey, I can't open the input file - %s \n", argv[1]);
      exit(2);
   }

   out = fopen (argv[2], "wb");          // open output file 
   if (out == NULL)
   {
      printf ("        Hey, I can't open the output file - %s \n", argv[2]);
      exit(3);
   }

   // Convert the filename to Lat/Lon and verify
   // Filename ares, for example, N42W074.hgt

   // 000000000011111111112
   // 012345678901234567890
   // N42W074.hgt

   LatDir = argv[1][0];    // North or South (N/S)
   LonDir = argv[1][3];    // East or West (E/W)
   argv[1][3] = 0;
   ptr = &argv[1][1];
   Latitude = atof(strtok(ptr, "\n"));
   argv[1][7] = 0;
   ptr = &argv[1][4];
   Longitude = atof(strtok(ptr, "\n"));

   // Extract SNEW
   if (LatDir == 'S')
   {
      Latitude = -Latitude;
   }

   if (LonDir == 'W')
   {
      Longitude = -Longitude;
   }

   // printf("\n Latitude 1: %lf\n", Latitude);
   // printf("\nLongitude 1: %lf\n", Longitude);

   xfact = distance(Latitude * PIE / 180.0, 
                    Longitude * PIE / 180.0, 
                    Latitude * PIE / 180.0, 
                    (Longitude - 1.0) * PIE / 180.0);
   
   yfact = distance(Latitude * PIE / 180.0, 
                    Longitude * PIE / 180.0, 
                    (Latitude + 1.0 ) *PIE / 180.0, 
                    Longitude * PIE / 180.0);
   
   // printf("\nyfact 1: %lf\n", yfact);

   xfact = xfact / _COLS * 1000.0;
   yfact = yfact / _ROWS * 1000.0;

   // printf("\nyfact 2: %lf\n",yfact);

   // Extract SNEW
   if (Latitude <=0 )
   {
      Latitude = -Latitude;
   }

   if (Longitude <=0 )
   {
      Longitude = -Longitude;
   }

   // Handle args
   i=3;
   while (i<argc)
   {
      switch (argv[i][1])
      {
         // Allows you to amplify the elevation by a factor
      case 'A':
      case 'a':
         ptr = argv[i+1];
         AmpFactor = atof(strtok(ptr, ","));
         i += 2;
         break;
   
         // Bias...force the altitude (thinness/thickness) up or down
         // This is especially nice when you amplify but don't want
         // the base to be too thick.
      case 'B':       
      case 'b':
         Bias = atoi(argv[i+1]);
         i += 2;
         break;
   
      case 'C':
      case 'c':
         ptr = argv[i+1];
               
         // Extract the coordinates
         LatLowerLeftCoord = abs(atof(strtok(ptr, ",")));
         LonLowerLeftCoord = abs(atof(strtok(NULL, ",")));
         LatUpperRightCoord = abs(atof(strtok(NULL, ",")));
         LonUpperRightCoord = abs(atof(strtok(NULL, ",")));
         i += 2;

         if(LatDir=='N')   //North
         {
            switch(LonDir)
            {
               case 'W':
                  // Compute indexes
                  ExtractionRow1 = (int)roundf(((float) 60.0 - (float) LatUpperRightCoord) * (float) 20.0) - 20;
                  ExtractionRow2 = (int)roundf(((float) 60.0 - (float) LatLowerLeftCoord) * (float) 20.0);
                  ExtractionColumn1 = _COLS - 1 - (int)roundf(((float) 60.0 - (float) LonUpperRightCoord) * (float) 20.0);
                  ExtractionColumn2 = _COLS - 1 -(int)roundf(((float) 60.0 - (float) LonLowerLeftCoord) * (float) 20.0) + 20;
                  break;
               case 'E':
                  // Compute indexes
                  ExtractionRow1 = (int)roundf(((float) 60.0 - (float) LatUpperRightCoord) * (float) 20.0) - 20;
                  ExtractionRow2 = (int)roundf(((float) 60.0 - (float) LatLowerLeftCoord) * (float) 20.0);
                  ExtractionColumn1 = _COLS - 1 - (int)roundf((float) LonUpperRightCoord * (float) 20.0) - 20;
                  ExtractionColumn2 = _COLS - 1 - (int)roundf((float) LonLowerLeftCoord * (float) 20.0);
                  break;
            }
         }
         else
         {
            switch(LonDir) //South
            {
               case 'W':
                  // Compute indexes
                  ExtractionRow1 = (int)roundf((float) LatUpperRightCoord * (float) 20.0);
                  ExtractionRow2 = (int)roundf((float) LatLowerLeftCoord * (float) 20.0) + 20;
                  ExtractionColumn1 = _COLS - 1 - (int)roundf(((float) 60.0 - (float) LonUpperRightCoord) * (float) 20.0);
                  ExtractionColumn2 = _COLS - 1 -(int)roundf(((float) 60.0 - (float) LonLowerLeftCoord) * (float) 20.0) + 20;
                  break;
               case 'E':
                  // Compute indexes
                  ExtractionRow1 = (int)roundf((float) LatUpperRightCoord * (float) 20.0);
                  ExtractionRow2 = (int)roundf((float) LatLowerLeftCoord * (float) 20.0) + 20;
                  ExtractionColumn1 = _COLS - 1 - (int)roundf((float) LonUpperRightCoord * (float) 20.0) - 20;
                  ExtractionColumn2 = _COLS - 1 - (int)roundf((float) LonLowerLeftCoord * (float) 20.0);
                  break;
            }
         }
         break;
   
         // Fix holes (linear fill)
      case 'F':
      case 'f':
         FixHoles = TRUE;
         i++;
         break;
   
         // Name of the file (STL ASCII)
      case 'N':
      case 'n':
         ptr=argv[i+1];
         strcpy(SolidName,strtok(ptr, " \n"));
         i += 2;
         break;
   
         // So we puke out all kinds of extra stuff
      case 'V':               // Verbose
      case 'v':
         Verbose = TRUE;
         i++;
         break;
   
         // So we puke out all kinds of extra stuff
      case 'H':              // Hole errors ON/OFF /H+ or /H-
      case 'h':
         switch (argv[i][2])
         {
         case '+':
            HoleErrors = TRUE;
            break;
         case '-':
            HoleErrors = FALSE;
            break;
         default:
            printf("\nNot a valid option: [%s]\n", argv[i]);
            break;
         }
         i++;
         break;
   
         // Subset corners. From lower left to upper right
      case 'S':
      case 's':
         ptr = argv[i+1];
                                   
         // Extract the coordinates (these variable names look wrong because the variables are doing double duty
         LonLowerLeft = atoi(strtok(ptr, ","));
         LatLowerLeft = atoi(strtok(NULL, ","));
         LonUpperRight = atoi(strtok(NULL, ","));
         LatUpperRight = atoi(strtok(NULL, ","));
         i += 2;
   
         ExtractionRow1 = _ROWS - LatUpperRight;
         ExtractionRow2 = _ROWS - LatLowerLeft - 1;                // Minus one
         ExtractionColumn1 = _COLS - LonUpperRight;
         ExtractionColumn2 = _COLS - LonLowerLeft - 1;             // Minus one
         break;
   
      case 'T':               // Output file: T+ = Text  or /T- = Binary
      case 't':
         switch (argv[i][2])
         {
         case '+':
            Binary = FALSE;
            break;
         case '-':
            Binary = TRUE;
            break;
         default:
            printf("\nNot a valid option: [%s]\n", argv[i]);
            break;
         }
         i++;
         break;
   
         // If we don't know what you entered...
      default:
         printf("\nNot a valid option: [%s]\n", argv[i]);
         i++;
         break;
      }
   }

   // Dump all the options entered (if verbose)
   if (Verbose)
   {
      //printf("\nExtractionRow1: %d", ExtractionRow1);
      //printf("\nExtractionRow2: %d", ExtractionRow2);
      //printf("\nExtractionColumn1: %d", ExtractionColumn1);
      //printf("\nExtractionColumn2: %d", ExtractionColumn2);
      printf("\nAmplification Factor: %.1f", AmpFactor);
      printf("\nBias: %ld", Bias);
      printf("\nLat: %.1f [%c]", Latitude, LatDir);
      printf("\nLon: %.1f [%c]", Longitude, LonDir);
      printf("\nLatitude  increment: %.2lf meters", yfact);
      printf("\nLongitude increment: %.2lf meters", xfact);
      // Compute Approximate Area
      xSide = (double) (ExtractionColumn2 - ExtractionColumn1) * xfact;
      ySide = (double) (ExtractionRow2 - ExtractionRow1) * yfact;
      printf("\nTotal Area (approximate): %.1lf Sq. Kilometers (%.1lf km X %.1lf km)", (xSide / 1000) * (ySide / 1000), (xSide / 1000), (ySide / 1000));
      printf("\n");
   }

   i = 0;

   // ****************************
   // Load the array with the grid
   // ****************************
   for (row=0; row < _ROWS; row++)
   {
      for (column=0; column< _COLS; column++)
      {
         grid[row][column] = ((unsigned short) fgetc(in)<<8) | (unsigned short) fgetc(in);

         // Is this a hole?
         if (grid[row][column] == -32768) // Is this a bad data point? Bad data is largest negative signed 16 bit value.
         {
            if (HoleErrors)
            {
               NumHoles++;
               printf("Bad data (i.e. hole) at [%d][%d]==[%d]\n", row, column, grid[row][column]);
            }
         }

         // Amplify the elevations by a factor
         if (AmpFactor != 1.0)
         {
            if (grid[row][column] != -32768)
            {
               grid[row][column] = (short) ((float) grid[row][column] * AmpFactor);
            }
         }

         if (Bias != 0)
         {
            if (grid[row][column] != -32768)
            {
               grid[row][column] += Bias;
            }
         }
      }
   }

   // ****************************
   // Flip the grid (Horizontally)
   // ****************************
   for (row=0; row < _ROWS; row++)
   {
      for (column=0; column < (_COLS / 2); column++)
      {
         stemp = grid[row][column];
         grid[row][column] = grid[row][_COLS - column];
         grid[row][_COLS - column] = stemp;
      }
   }


   // *********************************
   // If Verbose AND there are holes...
   // *********************************
   if (Verbose)
   {
      if (NumHoles)
      {
         printf("Number of holes: %d\n", NumHoles);
      }
   }

   // **********
   // Fix holes?
   // **********
   if (FixHoles && NumHoles)
   {
      FixedHoles = 0;
      // **************
      // Fill any holes (linear fill)
      // **************
      for (row=0; row < _ROWS; row++)
      {
         _leftbank = -1;
         _rightbank = -1;

         for (column=0; column < _COLS; column++)
         {
            if ((grid[row][column] == -32768) && (_leftbank == -1))
            {
               // Later handler the case when column = 0 or _COLS-1
               if (column != 0)
               {
                  _leftbank = column-1;
               }
            }
            else if ((grid[row][column] != -32768) && (_leftbank != -1) && (_rightbank == -1))
            {
               // Later handle the case when column = 0 or _COLS-1
               _rightbank = column;
               _slope = (int)(((float)grid[row][_rightbank] - (float)grid[row][_leftbank]) / ((float)_rightbank - (float)_leftbank));
               if (_slope == 0)
               {
                  _fill = 0;        // grid[row][_leftbank];
               }
               else
               {
                  _fill = _slope;
               }

               for (int tcol = _leftbank + 1, indx=1; tcol < _rightbank; tcol++, indx++)
               {
                  grid[row][tcol] = grid[row][_leftbank] + (indx * _fill);
                  FixedHoles++;
               }
               _leftbank = -1;
               _rightbank = -1;                                
            }
         }
      }
   }

   // ********************************
   // Report the number of holes fixed
   // ********************************
   if (Verbose)
   {
      if (FixedHoles)
      {
         printf("Holes fixed: %d\n",FixedHoles);
      }
   }

   // ******************
   // Write the STL file
   // ******************
   facets = 0;

   if (Verbose)
   {
      if (Binary)
      {
         printf("Writing binary file...\n");
      }
      else
      {
         printf("Writing file...\n");
      }
   }

   if (!Binary)
   {
      fprintf(out, "solid %s\n", SolidName);
   }
   else
   {
      // Binary file
      unsigned char header[STL_HEADER_LENGTH] = {0};
      fwrite(header, 1, STL_HEADER_LENGTH, out);
      fwrite(&facets, sizeof(facets), 1, out);
   }

   // *********************
   // Create Surface Facets
   // *********************
   for (row = ExtractionRow1; row < ExtractionRow2; row++)
   {
      for (column = ExtractionColumn1; column < ExtractionColumn2; column++)
      {
         // Two facets each time through the loop
         // Calculate the first facet
         v1.x = xfact * (float) column;                   // Vertex 
         v1.y = yfact * (float) row;
         v1.z = (float) grid[row][column];
        
         v2.x = xfact * (float) (column + 1);               // B
         v2.y = yfact * (float) row;
         v2.z = (float) grid[row][column + 1];

         v3.x = xfact * (float) column;                   // A
         v3.y = yfact * (float) (row+1);
         v3.z = (float) grid[row + 1][column];

         normal = calcNormal(v1.x, v3.x, v2.x);

         stlwrite(out, Binary, &v1, &v3, &v2, &normal);

         // Calculate the second facet (just the one new point)
         v1.x = xfact * (float) (column + 1);       // Vertex  (V2 is A, V3 is B)
         v1.y = yfact * (float) (row + 1);
         v1.z = (float) grid[row + 1][column + 1];

         normal = calcNormal(v1.x, v2.x, v3.x);

         stlwrite(out, Binary, &v1, &v2, &v3, &normal);

         facets += 2;
      }
   }       
        
   // ****************
   // Closing polygons
   // ****************

   // *********
   // Left side
   // *********
   for (row = ExtractionRow1; row < ExtractionRow2; row++)
   {
      column = ExtractionColumn1;

      // Two facets each time through the loop
      // Calculate the first facet
      v1.x = xfact * (float) column;           // B
      v1.y = yfact * (float) row;
      v1.z = (float) grid[row][column];
        
      v2.x = xfact * (float) column;           // Vertex
      v2.y = yfact * (float) row;
      v2.z = (float) 0.0;

      v3.x = xfact * (float) column;           // A
      v3.y = yfact * (float) (row + 1);
      v3.z = (float) grid[row + 1][column];

      normal = calcNormal(v2.x, v3.x, v1.x);

      stlwrite(out, Binary, &v2, &v3, &v1, &normal);

      // Calculate the second facet (just the one new point)
      v1.x = xfact * (float) column;   // B (V3 is the vertex, V2 is A)
      v1.y = yfact * (float) (row + 1);
      v1.z = (float) 0.0;
                                
      normal = calcNormal(v3.x, v2.x, v1.x);

      stlwrite(out, Binary, &v3, &v2, &v1, &normal);

      facets += 2;
   }       

   // **********
   // Right side
   // **********
   for (row = ExtractionRow1; row < ExtractionRow2; row++)
   {
      column = ExtractionColumn2;
 
      // Two facets each time through the loop
      // Calculate the first facet
      v1.x = xfact * (float) column;                   // A
      v1.y = yfact * (float) row;
      v1.z = (float) grid[row][column];
        
      v2.x = xfact * (float) column;                   // Vertex
      v2.y = yfact * (float) row;
      v2.z = (float) 0.0;

      v3.x = xfact * (float) column;                   // B
      v3.y = yfact * (float) (row + 1);
      v3.z = (float) grid[row + 1][column];

      normal = calcNormal(v2.x, v1.x, v3.x);

      stlwrite(out, Binary, &v2, &v1, &v3, &normal);

      // Calculate the second facet (just the one new point)
      v1.x = xfact * (float) column;           // A (V2 is B and the Vertex V3
      v1.y = yfact * (float) (row + 1);
      v1.z = (float) 0.0;
                                
      normal = calcNormal(v3.x, v1.x, v2.x);

      stlwrite(out, Binary, &v3, &v1, &v2, &normal);

      facets += 2;
   }       

   // ********
   // Top edge
   // ********
   for (column = ExtractionColumn1; column < ExtractionColumn2; column++)
   {
      row = ExtractionRow1;

      // Two facets each time through the loop
      // Calculate the first facet
      v1.x = xfact * (float) column;                 // A
      v1.y = yfact * (float) row;
      v1.z = (float) grid[row][column];
        
      v2.x = xfact * (float) (column);               // Vertex
      v2.y = yfact * (float) row;
      v2.z = (float) 0.0;

      v3.x = xfact * (float) (column + 1);             // B
      v3.y = yfact * (float) (row);
      v3.z = (float) grid[row][column + 1];

      // Write the first facet
      normal = calcNormal(v2.x, v1.x, v3.x);

      stlwrite(out, Binary, &v2, &v1, &v3, &normal);

      // Calculate the second facet (just the one new point)                   
      v1.x = xfact * (float) (column + 1);       // A (Vertex is V3 and V2 is B)
      v1.y = yfact * (float) (row);
      v1.z = (float) 0.0;
                                
      normal = calcNormal(v3.x, v1.x, v2.x);

      stlwrite(out, Binary, &v3, &v1, &v2, &normal);

      facets += 2;      
   }       

   // **********
   //Bottom edge
   // **********
   for (column = ExtractionColumn1; column < ExtractionColumn2; column++)
   {
      row = ExtractionRow2;

      // Two facets each time through the loop
      // Calculate the first facet
      v1.x = xfact * (float) column;                 // B
      v1.y = yfact * (float) row;
      v1.z = (float) grid[row][column];
        
      v2.x = xfact * (float) (column);               // Vertex
      v2.y = yfact * (float) row;
      v2.z = (float) 0.0;

      v3.x = xfact * (float) (column + 1);             // A
      v3.y = yfact * (float) (row);
      v3.z = (float) grid[row][column + 1];

      normal = calcNormal(v2.x, v3.x, v1.x);

      stlwrite(out, Binary, &v2, &v3, &v1, &normal);

      // Calculate the second facet (just the one new point)
      v1.x = xfact * (float) (column + 1);       // B (Vertex is V3 and A is V2) 
      v1.y = yfact * (float) (row);
      v1.z = (float) 0.0;
                                
      normal = calcNormal(v3.x, v2.x, v1.x);

      stlwrite(out, Binary, &v3, &v2, &v1, &normal);

      facets += 2;
   }       

   // *************
   // Create Bottom
   // *************
   for (row= ExtractionRow1; row < ExtractionRow2; row++)
   {
      for (column = ExtractionColumn1; column < ExtractionColumn2; column++)
      {
         //Two facets each time through the loop
         //Calculate the first facet
         v1.x = xfact * (float) column;                 //Vertex
         v1.y = yfact * (float) row;
         v1.z = (float) 0.0;
        
         v2.x = xfact * (float) (column+1);             //A
         v2.y = yfact * (float) row;
         v2.z = (float) 0.0;

         v3.x = xfact * (float) column;                 //B
         v3.y = yfact * (float) (row+1);
         v3.z = (float) 0.0;

         normal = calcNormal(v1.x, v2.x, v3.x);

         stlwrite(out, Binary, &v1, &v2, &v3, &normal);

         // Calculate the second facet (just the one new point)
         v1.x = xfact * (float) (column+1);               // Vertex (A is V3 and B is V2)
         v1.y = yfact * (float) (row+1);
         v1.z = (float) 0.0;

         normal = calcNormal(v1.x, v3.x, v2.x);

         stlwrite(out, Binary, &v1, &v3, &v2, &normal);

         facets += 2;
      }
   }       

   if (!Binary)
   {
      fprintf(out,"endsolid\n");
   }
   else
   {
      // now that we know the total number of facets, fill in the
      // appropriate field in the binary file
      fseek(out, STL_HEADER_LENGTH, SEEK_SET);
      fwrite(&facets, sizeof(facets), 1, out);
   }

   //Close files and free memory and anything else we have to clean up...do it here!
   if (in)
   {
      fclose(in);
   }
   if (out)
   {
      fclose(out);
   }
   if (Verbose)
   {
      printf("Facets: %d\n\n",facets);
   }

   return(0);

}
