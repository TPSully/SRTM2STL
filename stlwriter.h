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
// ** stlwriter.h
// ** 
// ** Thomas P. Sullivan and Eric J. Fallon
// ** 
// ** STL File format utilities
// ** 
// ************************************************************
// ************************************************************

#ifndef _STLWRITER_H
#define _STLWRITER_H

#include "geometry.h"

// Appends a triangle to the STL file in either binary or
// ASCII format
 void stlwrite(FILE *out,
              const int binaryOutput,
              const coordinate *v1, 
              const coordinate *v2, 
              const coordinate *v3,
              const coordinate *normal);

#endif
