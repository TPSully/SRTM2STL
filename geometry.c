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
// ** geometry.c
// ** 
// ** Thomas P. Sullivan and Eric J. Fallon
// ** 
// ** Misc geometry functions
// ** 
// ************************************************************
// ************************************************************

#include <math.h>
#include "geometry.h"

coordinate calcNormal(const float v, const float a, const float b)
{
   coordinate normal, uu, vv;

   //Take the differences to create VERTEX, A and B for Cross Product
   uu.x = a - v;
   uu.y = a - v;
   uu.z = a - v;

   vv.x = a - v;
   vv.y = a - v;
   vv.z = a - v;

   //Compute the normal
   normal.x = (uu.y * vv.z) - (uu.z * vv.y);
   normal.y = (uu.x * vv.z) - (uu.z * vv.x);
   normal.z = (uu.x * vv.y) - (uu.y * vv.x);

   return normal;
}


// From: http://home.att.net/~srschmitt/great_circle_route.html (page doesn't seem to exist anymore)


// Computing Great Circle Distance
// by Stephen R. Schmitt

// Converted to C by Thomas P. Sullivan

// Distance using Meeus approximation

// ****************************************************
// ** Function:         distance
// ** 
// ** Parameters: 
// **                           double precision start latitude
// **                           double precision start longitude
// **                           double precision end latitude
// **                           double precision end longitude
// ** Returns: 
// **                           double precision distance in meters
// ** 
// ** Notes: These are angles in radians NOT degrees
// ** 
// ****************************************************
double distance(double lat1, double lon1, double lat2, double lon2)
{

   double F = (lat1 + lat2) / 2.0;
   double G = (lat1 - lat2) / 2.0;
   double L = (lon1 - lon2) / 2.0;

   double sinG2 = sin(G) * sin(G);
   double cosG2 = cos(G) * cos(G);
   double sinF2 = sin(F) * sin(F);
   double cosF2 = cos(F) * cos(F);
   double sinL2 = sin(L) * sin(L);
   double cosL2 = cos(L) * cos(L);

   double S = sinG2 * cosL2 + cosF2 * sinL2;
   double C = cosG2 * cosL2 + sinF2 * sinL2;

   double w = atan(sqrt(S/C));
   double R = sqrt(S*C)/w;

   double a = 6378.137;                    // WGS-84 equatorial radius
   double f = 1.0/298.257223563;           // WGS-84 ellipsoid flattening factor

   double D = 2*w*a;
   double H1 = (3*R - 1)/(2*C);
   double H2 = (3*R + 2)/(2*S);

   double dist = D * (1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2);

   return dist;
}
