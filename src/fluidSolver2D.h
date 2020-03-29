#include "vector2D.h"
#include <afxwin.h>																	// MFC core and standard components
#include <afxext.h>																	// MFC extensions
#include <afxdisp.h>																// MFC Automation classes


#define NUM_THREADS 4																// Number of threads to use for parallel tasks
#define NULL_CELL 0																	// The various states that the cells can be in
#define AIR_CELL 1
#define SURFACE_FLUID_CELL 2
#define FULL_FLUID_CELL 3
#define SOLID_CELL 4

void fluidSolver_init(int, int, float, float, float, float, float);					// Receives all fluid properties and records them
unsigned char getState(int, int);													// Returns the state of the given cell
void setState(int, int, unsigned char);												// Sets the state of the given cell
bool cellIsFixed(int, int);															// Returns whether a cell is fixed (i.e. velocity values set by user and held constant)
void setFixed(int, int, bool);														// Sets whethe a cell is fixed

float getVel(int, int, int);														// Returns a velocity value specified by cell & face
void setVel(int, int, int, float);													// Sets a velocity valuespecified by cell & face
vector2D getVel(float, float);														// Returns interpolated velocity at arbitrary point in world space
float getVelH(int, int);															// Returns a horizontal velocity component by index
float getVelV(int, int);															// Returns a vertical velocity component by index

void fluidSolver_iterate(int infoX, int infoY, int com, int verbose, int);			// Moves the simulation forward by one given time step
void fluidSolver_cleanup();															// Frees dynamically allocated memory
		
extern float cellsize;																// The world space width of the square fluid cells;
extern int cellw, cellh;															// The number of fluid cells horizontally & vertically which make up the fluid domain;