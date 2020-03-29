#include "vector2D.h"

#include <iostream>
#include <math.h>


//
//	vector2D() - default constructor
//
vector2D::vector2D(void)												
{
	x=y=0;
}

//
//	vector2D(float, float) - constructor via component
//
vector2D::vector2D(float xIn, float yIn)									
{
	x=xIn;
	y=yIn;
}

//
//	operator[]() - for array-style access
//
float &vector2D::operator[](int c)
{
	if (c == 0) return x;
	else return y;
}

//
//	operator*()	- for scalar multiplication
//
vector2D vector2D::operator*(float scalar)
{
	vector2D result(x*scalar,y*scalar);
	
	return result;
}
