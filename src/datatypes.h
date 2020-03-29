#pragma once

//
//	rect - for describing a rectangle
//
struct rect
{
	float x, y, w, h;									// Dimensions in world space
	int i0, i1, j0, j1;									// The indices of the the block of cells whose centres are within the boundaries of this rectangle
};

//
//	fluidSource - describes an entity that will intoduce liquid to the system
//
struct fluidSource
{	
	rect shape;											// The rectangular shape of the source
	float vx, vy;										// The vector velocity of liquid originating from this source
	float t0, t1;										// The start and end times of at which this source is active
	bool hasVel;										// Whether this source imposes a velocity value (as apposed to being passive, allowing the velocity to be determined by the fluid simulation)
};

//
//	particle - desribes a simple 2D particle for massless advection
//
struct particle
{
	float x, y;											// Particle position
	bool draw;											// Whether to draw this particle
};
