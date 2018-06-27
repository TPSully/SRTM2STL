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
// ** stlwriter.c
// ** 
// ** Thomas P. Sullivan and Eric J. Fallon
// ** 
// ** STL File format utilities
// ** 
// ************************************************************
// ************************************************************

#include <stdio.h>
#include "stlwriter.h"

// Appends a triangle to the STL file in either binary or
// ASCII format
 void stlwrite(FILE *out,
              const int binaryOutput,
              const coordinate *v1, 
              const coordinate *v2, 
              const coordinate *v3,
              const coordinate *normal)
 {
    short dummy = 0;

    if (binaryOutput)
    {
       //Write binary
       fwrite(&normal->x, sizeof(float), (size_t) 1, out);
       fwrite(&normal->y, sizeof(float), (size_t) 1, out);
       fwrite(&normal->z, sizeof(float), (size_t) 1, out);
       fwrite(&v1->x, sizeof(float), (size_t) 1, out);
       fwrite(&v1->y, sizeof(float), (size_t) 1, out);
       fwrite(&v1->z, sizeof(float), (size_t) 1, out);
       fwrite(&v2->x, sizeof(float), (size_t) 1, out);
       fwrite(&v2->y, sizeof(float), (size_t) 1, out);
       fwrite(&v2->z, sizeof(float), (size_t) 1, out);
       fwrite(&v3->x, sizeof(float), (size_t) 1, out);
       fwrite(&v3->y, sizeof(float), (size_t) 1, out);
       fwrite(&v3->z, sizeof(float), (size_t) 1, out);
       fwrite(&dummy, sizeof(short), (size_t) 1, out);  //Atribute byte count (should be zero...at least, that's the recommendation)
    }
    else
    {
       fprintf(out,"   facet normal %f %f %f\n", normal->x, normal->y, normal->z);
       fprintf(out,"      outer loop\n");
       fprintf(out,"         vertex %f %f %f\n", v1->x, v1->y, v1->z);
       fprintf(out,"         vertex %f %f %f\n", v2->x, v2->y, v2->z);
       fprintf(out,"         vertex %f %f %f\n", v3->x, v3->y, v3->z);
       fprintf(out,"      endloop\n");
       fprintf(out,"   endfacet\n");
    }
 }
