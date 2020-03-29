#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "fluidSolver2D.h"

#define MAX_SPEED 2.0f																	// Limit speed to this value to avoid instability
#define SURFACE_FRICTION 0.125															// Determines tangent speed at fluid/solid interfaces (0: speed = nearest parallel velocity, 1: speed = 0.0)
#define B0 1.7																			// Relaxation coefficient
#define WIDTH(c)   ((c == 0)? cellw + 1 : cellw)										// Convenience macro to determine 2D array width based on velocity component
#define HEIGHT(c)   ((c == 0)? cellh : cellh + 1)										// Convenience macro to determine 2D array width based on velocity component

struct threadInfo																		// For providing per-thread information
{
	int i0, i1;
	int j0, j1;

	vector2D offset;
	int c;
};

int neighbour_i[4] = { -1, 0, 1, 0 };													// These arrays allow to iterate over a cell's neighbours by converting a 0-3 index to neighbour cell indices
int neighbour_j[4] = { 0, 1, 0, -1 };

float gravityVector[2] = {0.0, -9.8};													// Gravitational acceleration vector (usually points downwards)
float dt;																				// iteration timestep
float viscosity;																		// Fluid viscosity
float B;																				// Ratio between delta pressure & divergence

unsigned char *state;																	// Data structure for recording the state of each cell
bool *fixedCells;																		// Data structure for recording cells whose velocity values are fixed									
float *velocity[2];																		// 2 arrays for horizontal & vertical velocity components
float *velocity1[2];																	// Auxiliary arrays for velocity values, to allow iteration
float *pressure;																		// Array of pressure values for each cell 

CWinThread *threadObj1[NUM_THREADS];													// For managing threads
threadInfo threadInfosX[NUM_THREADS];
threadInfo threadInfosY[NUM_THREADS];

float *createFloatField(int, int);														// Creates a 2D float array with the dimensions specified, and initialises to 0
void generateThreadInfos(threadInfo *, int, int, float, float, int);					// Creates NUM_THREADS thread info objects, dividing the fluid domain into columns 
float interpolateVelocityComponentValue(float *, int, int, float, float, float, float);	// Returns the value of a velocity field at an arbitrary point, via interpolation
void applyBoundaryConditions();															// Applies all boundary conditions 
void makeSurfaceCellDivergenceFree(int, int);											// Appiles cell boundary conditions for a particular surface cell
UINT iterateVelocity(LPVOID pParam);													// Computes velocity field at t + 1, from field at time t. (uses Navier-Stokes eqn)
void conservationOfMass(float dt, int verbose);											// Applies conservation of mass by setting all fluid cells to be divergenve free


void fluidSolver_init(int cellw0, int cellh0, float cellsize0, float viscosity0, float d, float div, float dt0)
{
	// Record the main simualation settings in global variables:
	cellw = cellw0;
	cellh = cellh0;
	cellsize = cellsize0;
	dt = dt0;
	viscosity = viscosity0;
	
	// The fluid velocity field is stored in 2 arrays, for horizontal and vertical components.
	velocity[0] = createFloatField(cellw + 1, cellh);
	velocity[1] = createFloatField(cellw, cellh + 1);

	// Each iteration reads velocity values from the above arrays, and outputs the new velocitiy values into auxiliary versions:
	velocity1[0] = createFloatField(cellw + 1, cellh);
	velocity1[1] = createFloatField(cellw, cellh + 1);

	// Pressure values are stored at the centre of each fluid cell:
	pressure = createFloatField(cellw, cellh);

	// The state array records the state of each cell: air, surface fluid, full fluid or solid:
	state = new unsigned char[(cellw * cellh)];
	
	// Simple map recording which cells have their velocities fixed because they are fluid source cells:
	fixedCells = (bool *)(malloc(cellw * cellh *sizeof(bool)));
	for (int i = 0; i < (cellw * cellh); i++) fixedCells[i] = false;

	// Divide the velocity fields up into groups for multithreading:
	generateThreadInfos(threadInfosX, cellw + 1, cellh, 0.0, 0.5f * cellsize, 0);
	generateThreadInfos(threadInfosY, cellw, cellh + 1, 0.5f * cellsize, 0.0, 1);
}

//
//	getState() - returns the state of the speciied cell
//
unsigned char getState(int i, int j)
{
	// If queried on a non existant cell outside of the fluid domain, return a non-state:
	if((i < 0) || (i >= cellw) || (j < 0) || (j >= cellh)) return NULL_CELL;

	return state[j*cellw + i];
}

//
//	setState() - set the state of the specified cell
//
void setState(int i, int j, unsigned char val)
{
	if((i >= 0) && (i < cellw) && (j >= 0) && (j < cellh)) state[j*cellw + i] = val;
}

//
//	cellIsFixed() -	returns whether the given cell is fixed
//
bool cellIsFixed(int i, int j)
{
	if((i < 0) || (i >= cellw) || (j < 0) || (j >= cellh)) return false;
	return fixedCells[j * cellw + i];
}

//
//	setFixed() - sets the fixed status of the given cell
//
void setFixed(int i, int j, bool val)
{
	if((i >= 0) && (i < cellw) && (j >= 0) && (j < cellh)) fixedCells[j * cellw + i] = val;
}

//
//	getVel() - Returns a velocity value specified by cell & face
//
float getVel(int Xvox, int Yvox, int face)
{
	if(face == 0) return velocity[0][Yvox * (cellw + 1) + Xvox];
	if(face == 1) return velocity[1][(Yvox + 1) * cellw + Xvox];
	if(face == 2) return velocity[0][Yvox * (cellw + 1) + Xvox + 1];
	if(face == 3) return velocity[1][Yvox * cellw + Xvox];

	return 0;
}

//
//	setVel() -	Sets a velocity valuespecified by cell & face
//
void setVel(int Xvox, int Yvox, int face, float val)
{
	if(face == 0) velocity[0][Yvox * (cellw + 1) + Xvox] = val;
	if(face == 1) velocity[1][(Yvox + 1) * cellw + Xvox] = val;
	if(face == 2) velocity[0][Yvox * (cellw + 1) + Xvox + 1] = val;
	if(face == 3) velocity[1][Yvox * cellw + Xvox] = val;
}

//
//	getVel(float, float) - returns the vector velocity value at any point in the fluid domain. (Computed by linearly interpolating between the nearest recorded values)
//
vector2D getVel(float x, float y)
{
	float v00_i, v00_j;
	float px, py;
	float valX, valY;

	float xc = std::clamp(x, 0.0f, cellw*cellsize);
	float yc = std::clamp(y, 0.0f, cellh*cellsize);

	// Compute the horizontal and vertical velocity components at (x,y) by interpolating:
	float h = interpolateVelocityComponentValue(velocity[0], cellw + 1, cellh, 0.0, -0.5*cellsize, xc, yc);
	float v = interpolateVelocityComponentValue(velocity[1], cellw, cellh + 1, -0.5*cellsize, 0.0, xc, yc);

	return { h, v };
}

//
//	getVelH(int, int) - returns a horizontal velocity component value.  Accessed directly by specifiying the i,j indices of a particular recorded value
//
float getVelH(int i, int j) { return velocity[0][j * (cellw + 1) + i]; }

//
//	getVelV(int, int) - returns a vertical velocity component value.  Accessed directly by specifiying the i,j indices of a particular recorded value
//
float getVelV(int i, int j) { return velocity[1][j * cellw + i]; }

//
//	iterate() - Computes 1 iteration of fluid simulation over one timestep
//
void fluidSolver_iterate(int infoX, int infoY, int com, int verbose, int it)
{
	// Compute delta pressure to divergence ratio based on relaxation coefficient, cell size & timestep:
	B = B0 / (2.0 * dt * ((1.0 / pow(cellsize, 2)) + (1.0 / pow(cellsize, 2))));

	// Some velocity values at fluid/solid and fluid/air boundaries have fixed values:
	applyBoundaryConditions();
	
	// Compute the horizontal velocities in parallel:
	for (int t = 0; t < NUM_THREADS; t++) threadObj1[t] = AfxBeginThread(iterateVelocity, &threadInfosX[t]);
	for (int t = 0; t < NUM_THREADS; t++) WaitForSingleObject(threadObj1[t]->m_hThread, INFINITE);

	// Compute the vertical velocities in parallel:
	for (int t = 0; t < NUM_THREADS; t++) threadObj1[t] = AfxBeginThread(iterateVelocity, &threadInfosY[t]);
	for (int t = 0; t < NUM_THREADS; t++) WaitForSingleObject(threadObj1[t]->m_hThread, INFINITE);
	
	// Now that the result of this iteration has been computed (in the auxiliary arrays), swap the arrays so that the auxiliary arrays become the main arrays:
	for(int c = 0; c < 2; c++)
	{
		float *tmp = velocity[c];
		velocity[c] = velocity1[c];
		velocity1[c] = tmp;
	}

	// Apply conservation of mass by ensuring that all fluid cells are divergence free:
	conservationOfMass(dt, verbose);

	// Apply boundary conditions again, now that velocity values have been updated:
	applyBoundaryConditions();
}

//
//	fluidSolver_cleanup() - Frees dynamically allocated memory
//
void fluidSolver_cleanup()
{
	free(velocity[0]);
	free(velocity[1]);
	free(velocity1[0]);
	free(velocity1[1]);
	free(pressure);

	delete[] state;
	free(fixedCells);
}

//
//	createFloatField() -  Creates a 2D float array with the dimensions specified, and initialises to 0
//
float *createFloatField(int w, int h)
{
	float *ptr = (float *)(malloc(w * h * sizeof(float)));
	for (int i = 0; i < w * h; i++) ptr[i] = 0.0;

	return ptr;
}

//
//	generateThreadInfos() -  Creates NUM_THREADS thread info objects, dividing the fluid domain into columns 
//
void generateThreadInfos(threadInfo *threadInfos, int width, int height, float offsetX, float offsetY, int component)
{
	for (int t = 0; t < NUM_THREADS; t++)
	{
		threadInfos[t].i0 = (int)((float)(t) * (float)(width) / (float)(NUM_THREADS));
		threadInfos[t].i1 = (int)((float)(t + 1) * (float)(width) / (float)(NUM_THREADS)) - 1;

		threadInfos[t].j0 = 0;
		threadInfos[t].j1 = height - 1;

		threadInfos[t].offset.x = offsetX;
		threadInfos[t].offset.y = offsetY;
		threadInfos[t].c = component;
	}
}

//
//	interpolateVelocityComponentValue() - Returns the value of a velocity field at an arbitrary point, via interpolation
//
float interpolateVelocityComponentValue(float *field, int w, int h, float offsetX, float offsetY, float x, float y)
{
	// Determine the 4 data points closest to (x,y): 
	int i0 = std::clamp((int)((x + offsetX) / cellsize), 0, w - 1);
	int i1 = std::clamp((int)((x + offsetX) / cellsize) + 1, 0, w - 1);
	int j0 = std::clamp((int)((y + offsetY) / cellsize), 0, h - 1);
	int j1 = std::clamp((int)((y + offsetY) / cellsize) + 1, 0, h - 1);

	// Convert (x,y) to normalized parametric space between the 4 closest data points:
	float px = (x + offsetX) / cellsize - i0;
	float py = (y + offsetY) / cellsize - j0;

	// Interpolate & return:
	return((1.0 - px) * (1.0 - py) * field[j0 * w + i0] + (px) * (1.0 - py) * field[j0 * w + i1] + (1.0 - px) * (py)* field[j1 * w + i0] + (px) * (py)* field[j1 * w + i1]);
}

//
//	applyBoundaryConditions() - Applies all boundary conditions
//
void applyBoundaryConditions()
{
	// Loop over all cells and apply boundary conditions where necessary:
	for (int i = 0; i < cellw; i++)
	{
		for (int j = 0; j < cellh; j++)
		{
			// Set pressure to atmospheric in empty cells:
			if (getState(i, j) == AIR_CELL) pressure[j * cellw + i] = 0.0;

			// For each surface cell, propagate the velocity values from the fluid faces to the non fluid faces:
			if(getState(i, j) == SURFACE_FLUID_CELL)
			{
				makeSurfaceCellDivergenceFree(i, j);
				pressure[j * cellw + i] = 0.0;
			}

			// For solid cells, set the normal velocity to 0, the normal pressure gradient to 0, and the tangent velocity according to the user defined surface friction:
			if(getState(i, j) == SOLID_CELL)
			{
				for (int f = 0; f < 4; f++)
				{
					// If this face is on the boundary of the fluid domain, set its velocty to zero:
					if (getState(i + neighbour_i[f], j + neighbour_j[f]) == NULL_CELL)
					{
						setVel(i, j, f, 0.0);
						continue;
					}
					
					// If this face is adjacent to a non-solid cell, then its velocity is a surface normal velocity. Set it to zero, and set the pressure to be equal on both sides of the face:
					if(getState(i + neighbour_i[f], j + neighbour_j[f]) != SOLID_CELL)
					{
						setVel(i, j, f, 0.0);
						pressure[j * cellw + i] = pressure[(j + neighbour_j[f]) * cellw + i + neighbour_i[f]];

						continue;
					}

					// If this face's neighbouring cell is another solid cell, then this face represents a 'tangent velocity'.
					if(getState(i + neighbour_i[f], j + neighbour_j[f]) == SOLID_CELL)
					{
						// If the 2 neighbouring cells on the side perpendicular to this face are both non-solid cells 
						// Set this face's velocity such that the tangent velocity at the solid/non solid interface adheres to the user defined surface friction.
						int ai = i + neighbour_i[(f + 1) % 4];
						int aj = j + neighbour_j[(f + 1) % 4];
						int bi = i + neighbour_i[(f + 1) % 4] + neighbour_i[f];
						int bj = j + neighbour_j[(f + 1) % 4] + neighbour_j[f];

						if ((getState(ai, aj) != NULL_CELL) && (getState(ai, aj) != SOLID_CELL) && (getState(bi, bj) != NULL_CELL) && (getState(bi, bj) != SOLID_CELL))
						{
							setVel(i, j, f, (1.0 - 2.0*SURFACE_FRICTION) * getVel(ai, aj, f));
						}

						// In all other cases, set this face's velocity to 0.0:
						else setVel(i, j, f, 0.0);
					}
				}
			}
		}
	}
}

//
//	makeSurfaceCellDivergenceFree() - Appiles cell boundary conditions for a particular surface cell
//
void makeSurfaceCellDivergenceFree(int i, int j)
{
	// Compute the net influx of velocity into this cell:
	float influx = getVel(i, j, 0) - getVel(i, j, 1) - getVel(i, j, 2) + getVel(i, j, 3);

	// Determine how many of this surface cell's faces are shared with fluid cells (full or surface):
	int fluidFaces[4] = { 0, 0, 0, 0 };
	int numFluidFaces = 0;

	for (int f = 0; f < 4; f++)
	{
		if ((getState(i + neighbour_i[f], j + neighbour_j[f]) == FULL_FLUID_CELL) || (getState(i + neighbour_i[f], j + neighbour_j[f]) == SURFACE_FLUID_CELL))
		{
			fluidFaces[f] = 1;
			numFluidFaces++;
		}
	}

	// Determine how many of this surface cell's faces are shared with empty/air cells:
	int openFaces[4] = { 0, 0, 0, 0 };
	int numOpenFaces = 0;

	for (int f = 0; f < 4; f++)
	{
		if (((getState(i + neighbour_i[f], j + neighbour_j[f]) == NULL_CELL) || (getState(i + neighbour_i[f], j + neighbour_j[f]) == AIR_CELL)) && (cellIsFixed(i + neighbour_i[f], j + neighbour_j[f]) == 0))
		{
			openFaces[f] = 1;
			numOpenFaces++;
		}
	}

	// If there is only one open face, then pass the net influx velocity from the other 3 faces, back out through this face:
	if (numOpenFaces == 1)
	{
		for (int f = 0; f < 4; f++)
		{
			if (openFaces[f] == 1) setVel(i, j, f, getVel(i, j, f) + ((f % 3 == 0) ? -1.0 : 1.0) * influx);
		}
	}

	// If there are two open faces, the next step depends on whether the faces are adjacent or opposite:
	if (numOpenFaces == 2)
	{
		// If they are adjacent, then simply set each to the same value as the face opposite it:
		if ((openFaces[0] == openFaces[1]) || (openFaces[1] == openFaces[2]))
		{
			for (int f = 0; f < 4; f++)
			{
				if (openFaces[f] == 1) setVel(i, j, f, getVel(i, j, (f + 2) % 4));
			}
		}

		// If they are opposite, then they both share equally the net influx, and pass it out:
		else
		{
			for (int f = 0; f < 4; f++)
			{
				if (openFaces[f] == 1) setVel(i, j, f, getVel(i, j, f) + ((f % 3 == 0) ? -1.0 : 1.0) * 0.5 * influx);
			}
		}
	}

	// If there are 3 open faces, then the open face that is opposite the 1 non-open face will copy its velocity value, the other two faces values are untouched:
	if (numOpenFaces == 3)
	{
		for (int f = 0; f < 4; f++)
		{
			if (fluidFaces[f] == 1) setVel(i, j, (f + 2) % 4, getVel(i, j, f));
		}
	}
}

//
//	iterateVelocity() - Computes velocity field at t + 1, from field at time t. (uses Navier-Stokes eqn)
//
UINT iterateVelocity(LPVOID pParam)
{
	threadInfo *info = (threadInfo *)(pParam);										// Object that contains information for this thread

	vector2D offset = { info->offset[0], info->offset[1] };							// For converting from velocity data index space to world space
	int c = info->c;																// Which cartesian velocity component is being handled by this thread

	vector2D pos;																	// The world space position of each velocity data point
	float gravityTerm;																// The Navier-Stokes gravity component
	float viscousTerm;																// The Navier-Stokes viscous component
	float convectionTerm;															// The Navier-Stokes convection component
	float pressureTerm;																// The Navier-Stokes pressure component
	vector2D v_forward, v_backward, v_lat_backward, v_lat_forward;					// For neghbouring velocity values when computing derivatives

	// Use the Navier-Stokes equation to compute the new velocity at each point in this thread's group:
	for(int i = info->i0; i <= info->i1; i++)																
	{
		for(int j = info->j0; j <= info->j1; j++)
		{
			// Compute the world space position of this velocity component data point:
			pos = { offset.x + (float)(i)* cellsize, offset.y + (float)(j)* cellsize };

			// Gravity term: (Gravity only applies to velocity data points that are near a cell containing fluid:
			gravityTerm = 0.0;
			for(int j1 = j - 1; j1 <= (j + (1 - c)); j1++)
			{
				for (int i1 = i - 1; i1 <= (i + c); i1++)
				{
					if((getState(i1, j1) == FULL_FLUID_CELL)||(getState(i1, j1) == SURFACE_FLUID_CELL)) gravityTerm = gravityVector[c];
				}
			}

			// Viscous term:
			viscousTerm = viscosity * (1.0 / pow(cellsize, 2)) * (
					velocity[c][std::clamp(j + 1, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i, 0, WIDTH(c)-1)] - 2.0 * velocity[c][std::clamp(j, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i, 0, WIDTH(c)-1)] + velocity[c][std::clamp(j - 1, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i, 0, WIDTH(c)-1)]
				+	velocity[c][std::clamp(j, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i + 1, 0, WIDTH(c)-1)] - 2.0 * velocity[c][std::clamp(j, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i, 0, WIDTH(c)-1)] + velocity[c][std::clamp(j, 0, HEIGHT(c)-1)*WIDTH(c) + std::clamp(i - 1, 0, WIDTH(c)-1)]
				);

			// Convection term:
			v_backward = getVel(pos.x - 0.5*cellsize * (1 - c), pos.y - 0.5*cellsize * c);
			v_forward = getVel(pos.x + 0.5*cellsize * (1 - c), pos.y + 0.5*cellsize * c);
			v_lat_backward = getVel(pos.x - 0.5*cellsize * c, pos.y - 0.5*cellsize * (1 - c));
			v_lat_forward = getVel(pos.x + 0.5*cellsize * c, pos.y + 0.5*cellsize * (1 - c));
			convectionTerm = (1.0 / cellsize) * (pow(v_backward[c], 2) - pow(v_forward[c], 2) + v_lat_backward.x*v_lat_backward.y - v_lat_forward.x*v_lat_forward.y);

			// Pressure term:
			pressureTerm = (1.0 / cellsize) * (pressure[std::clamp(j - c, 0, cellw - 1) * cellw + std::clamp(i - (1 - c), 0, cellw - 1)] - pressure[std::clamp(j, 0, cellw - 1) * cellw + std::clamp(i, 0, cellw - 1)]);

			// Compute the new velocity value at this point:
			float newVel = velocity[c][j * WIDTH(c) + i] + (gravityTerm + viscousTerm + convectionTerm + pressureTerm) * dt;

			// Limit the velocity to a certain magnitude to avoid instability:
			newVel = std::clamp(newVel, -1.0f * MAX_SPEED, MAX_SPEED);

			// Do not adjust velocities that bound fixed cells or solid cells:
			if((cellIsFixed(i - (1 - c), j - c)) || (cellIsFixed(i, j))) newVel = velocity[c][j * WIDTH(c) + i];
			if((getState(i - (1 - c), j - c) == SOLID_CELL) || (getState(i, j) == SOLID_CELL)) newVel = velocity[c][j * WIDTH(c) + i];

			// Place the new velocity value into the other buffer:
			velocity1[c][j * WIDTH(c) + i] = newVel;
		}
	}

	return 0;
}

//
//	conservationOfMass() - Applies conservation of mass by setting all fluid cells to be divergenve free
//
void conservationOfMass(float dt, int verbose)
{
	int openFaces[4];													// For recording which of a cell's faces are free to have their velocity value changed
	int numOpenFaces;													// For counting the number of free faces
	float influx, dp, velocityChange;									// For computing each cell's divergence properties

	// Make the fluid divergence free by sequentially eliminating divergence in each cell, then repeating several times:
	for(int a = 0; a < 100; a++)
	{
		// Make every cell in the grid divergence free:
		for(int i = 0; i < cellw; i++)
		{
			for(int j = 0; j < cellh; j++)
			{
				// No need to apply conservation of mass to fixed cells & solid cells:
				if((cellIsFixed(i, j) == 1)||(getState(i, j) == SOLID_CELL)) continue;

				// Compute the total net influx velocity of this cell, the resulting change in pressure and change in outgoing velocity:
				float influx = getVel(i, j, 0) - getVel(i, j, 1) - getVel(i, j, 2) + getVel(i, j, 3);
				dp = B * influx;
				velocityChange = (dt / cellsize) * dp;

				// For air cells, eliminate divergence by changing velocities on faces bordering other air cells, but holding all other face velocities fixed:
				if(getState(i, j) == AIR_CELL) pressure[j * cellw + i] = 0.0;

				// Full cells eliminate divergence by passing their net influx velocity back out through all 4 faces, except faces that border solid cells and fixed cells:
				if(getState(i, j) == FULL_FLUID_CELL)
				{
					// Determmine how many of the four faces are free to have their velocity adjusted:
					numOpenFaces = 0;
					for (int f = 0; f < 4; f++)
					{
						int i1 = i + neighbour_i[f];
						int j1 = j + neighbour_j[f];

						// If the neghbouring cell not solid and not fixed, then this face is open:
						if((getState(i1, j1) != SOLID_CELL) && (cellIsFixed(i1, j1) == 0))
						{
							openFaces[f] = 1;
							numOpenFaces++;
						}
						else openFaces[f] = 0;
					}
					for (int f = 0; f < 4; f++)
					{
						if (openFaces[f] == 1) setVel(i, j, f, getVel(i, j, f) + (((f % 3) == 0) ? -1.0 : 1.0) * velocityChange * 4.0 / (float)(numOpenFaces));
					}

					// Update the pressure value for this cell:
					pressure[j * cellw + i] += dp;
				}

				// Surface cells acheive zero divergence by copying the velocities of full/surface fluid faces to the open faces opposite. Any net influx is passed back out through the remaining faces:
				if (getState(i, j) == SURFACE_FLUID_CELL)
				{
					makeSurfaceCellDivergenceFree(i, j);

					// Set pressure to atmospheric in surface cells:
					pressure[j * cellw + i] = 0.0;
				}
			}
		}
	}
}