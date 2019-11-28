#pragma once

#include <iostream>
#include <cmath>
#include <time.h>

const double EPSILON = 0.000001;

struct ITRIANGLE
{
	int pt1, pt2, pt3;
};
struct IEDGE
{
	int pt1, pt2;
};
struct XYZ
{
	double x, y, z;
};

int XYZCompare(const void *v1, const void *v2);
int Triangulate(int nv, XYZ pxyz[], ITRIANGLE v[], int &ntri);
bool CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3,
	double &xc, double &yc, double &r);
