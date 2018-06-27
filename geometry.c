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

coordinate calcNormal2(const float v, const float a, const float b)
{
   coordinate normal, uu, vv;

   //Take the differences to create VERTEX, A and B for Cross Product
   uu.x = a - v;
   uu.y = a - v;
   uu.z = a - v;

   vv.x = b - v;
   vv.y = b - v;
   vv.z = b - v;

   //Compute the normal
   normal.x = (uu.y * vv.z) - (uu.z * vv.y);
   normal.y = (uu.x * vv.z) - (uu.z * vv.x);
   normal.z = (uu.x * vv.y) - (uu.y * vv.x);

   return normal;
}


coordinate calcNormal22(coordinate *v, coordinate *a, coordinate *b)
{
   static coordinate normal, uu, vv;
   static coordinate t;

   //Take the differences to create VERTEX, A and B for Cross Product
   uu.x = a->x - v->x;
   uu.y = a->y - v->y;
   uu.z = a->z - v->z;

   vv.x = b->x - v->x;
   vv.y = b->y - v->y;
   vv.z = b->z - v->z;

   //Compute the normal
/* tps 5-30-2017
   normal.x = (uu.y * vv.z) - (uu.z * vv.y);
   normal.y = -1*((uu.x * vv.z) - (uu.z * vv.x));
   normal.z = (uu.x * vv.y) - (uu.y * vv.x);
*/
   normal.x = (uu.y * vv.z) - (uu.z * vv.y);
   normal.y = ((uu.z * vv.x) - (uu.x * vv.z));
   normal.z = (uu.x * vv.y) - (uu.y * vv.x);

   //For now, make the normal vector zero and let the program(s) fix it
//#define CALC_NORMAL
#ifdef CALC_NORMAL
   t.x = normal.x;
   t.y = normal.y;
   t.z = normal.z;
#else
   t.x = 0.0;
   t.y = 0.0;
   t.z = 0.0;
#endif
   return t;
}

//This method, straight off the internet (stack exchange), says:

//Let p1 = (x1,y1,z1)
//and p2 = (x2,y2,z2)
//and p3 = (x3,y3,z3)

//then the normal vector i s given by:

//nx  (y2y1)(z3z1)(y3y1)(z2z1)
//ny  (z2z1)(x3x1)(x2x1)(z3z1)
//nz  (x2x1)(y3y1)(x3x1)(y2y1)


coordinate calcNormal(coordinate *p1, coordinate *p2, coordinate *p3)
{
   //static coordinate normal, uu, vv;
   static coordinate t;

   t.x = ((p2->y - p1->y)*(p3->z - p1->z))-((p3->y - p1->y)*(p2->z - p1->z));
   t.y = ((p2->z - p1->z)*(p3->x - p1->x))-((p2->x - p1->x)*(p3->z - p1->z));
   t.z = ((p2->x - p1->x)*(p3->y - p1->y))-((p3->x - p1->x)*(p2->y - p1->y));

   return t;
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
