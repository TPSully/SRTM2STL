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
// ** geometry.h
// ** 
// ** Thomas P. Sullivan and Eric J. Fallon
// ** 
// ** Misc geometry functions
// ** 
// ************************************************************
// ************************************************************


#ifndef _GEOMETRY_H
#define _GEOMETRY_H


// 3 dimensional Cartesian coordinate 
typedef struct {
  float x, y, z;
} coordinate;


// calculates normal vector
coordinate calcNormal(const float v, const float a, const float b);

// Computes Great Circle Distance
double distance(double lat1, double lon1, double lat2, double lon2);

#endif
