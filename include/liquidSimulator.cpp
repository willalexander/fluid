#include "datatypes.h"
#include "fluidSolver2D.h"
#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <ostream>
#include <random>
#include <algorithm>
#include "videoOut.h"
#include "liquidSimulator.h"

#define PARTICLES_PER_CELL 160
#define OUT_VIDEO_WIDTH 800
#define DRAW_ITERATION_NUMBER false
#define DRAW_CELLS false
#define DRAW_SOLIDS true
#define DRAW_VELOCITIES false
#define DRAW_PARTICLES true
#define VELOCITY_DRAW_SCALE 0.001

// For dividing cells up for multi threading:
struct threadInfo
{
	int id;
	int start, end;
};

int neighbourIndex_i[4] = { -1, 0, 1, 0 };			// These arrays allow to iterate over a cell's neighbours by using a 0-3 index to access neighbour cell indices
int neighbourIndex_j[4] = { 0, 1, 0, -1 };

char yellow[3] = { 0, 255, 255 };					// For use in drawing velocity components

// Binary array for drawing digits 0-9 to the frame buffer
bool digits[10][15] = {		{ true, true, true, true, false, true, true, false, true, true, false, true, true, true, true },
							{ false, true, false, false, true, false, false, true, false, false, true, false, false, true, false },
							{ true, true, true, true, false, false, true, true, true, false, false, true, true, true, true },
							{ true, true, true, false, false, true, true, true, true, false, false, true, true, true, true },
							{ false, false, true, false, false, true, true, true, true, true, false, true, true, false, true },
							{ true, true, true, false, false, true, true, true, true, true, false, false, true, true, true },
							{ true, true, true, true, false, true, true, true, true, true, false, false, true, true, true },
							{ false, false, true, false, false, true, false, false, true, false, false, true, true, true, true },
							{ true, true, true, true, false, true, true, true, true, true, false, true, true, true, true },
							{ false, false, true, false, false, true, true, true, true, true, false, true, true, true, true } };

std::default_random_engine generator(0);			// For random numbers

float dimw, dimh;									// The width & height of the domain in metres
float cellsize;										// The width of the square fluid cells;
int cellw, cellh;									// The number of fluid cells horizontally & vertically which make up the fluid domain;
float duration;										// Simulation duration (in seconds)
float timestep;										// Duration of each iteration (seconds)
int outVideoHeight;									// The hight in pixels of the frame buffer
char *hexBuffer;									// The frame buffer for displaying the fluid:
int cellPixelWidth;									// cell size measured in pixels
float pixelsPerMetre;								// Conversion factor between world space and pixel space

std::vector<rect> *solidObjects;					// Container of solid object instances
std::vector<fluidSource> *fluidSources;				// Container of fluid source instances
std::vector<particle> **particles;					// Data structure for particles

threadInfo threadInfos[NUM_THREADS];				// Per-thread info for multithreading
CWinThread *threadObj[NUM_THREADS];

float randFloat();									// Generates a random number between 0 and 1
void setup();										// Applies information about fluid sources & solid objects to the fluid grid
void determineCellCoverage(rect *);					// Computes which cells are within the boundaries of a rectangle
void processFluidSources(float);					// Applies fluid source properties as velocity values on the simulation grid
void determineCellStates();							// Sets the state of each cell based upon positions of fluid particles
void iterateParticles();							// Advects all fluid particles using the velocity field
UINT advectParticles_thread(LPVOID);				// Per-thread function for advecting a subset of the particles

void visualiseFluid(int);							// For drawing the fluid & other data to the framebuffer
UINT visualiseFluid_thread(LPVOID);
void drawIterationNumber(int);
void drawVelocityLine(int, int, int, int, int);


//
//	initialise() - called externally. This is the start point and sets up the fundamental settings for the simulation
//
void initialise(float dimw0, float dimh0, float cellsize0, float duration0, float timestep0, std::vector<rect> *solidObjects0, std::vector<fluidSource> *fluidSources0)
{
	// Record fluid dimension info, and adjust to account for cell size:
	dimw = dimw0;
	dimh = dimh0;

	cellsize = cellsize0;

	cellw = (int)(ceil(dimw / cellsize));
	cellh = (int)(ceil(dimh / cellsize));

	dimw = cellw * cellsize;
	dimh = cellh * cellsize;

	// During simulation some operations can be computed in parallel, so split the cells into groups for multithreading:
	for (int t = 0; t < NUM_THREADS; t++)
	{
		threadInfos[t].id = t;
		threadInfos[t].start = (int)(t * (float)(cellw) / (float)(NUM_THREADS));
		threadInfos[t].end = (int)((t + 1) * (float)(cellw) / (float)(NUM_THREADS)-1);
	}

	// Determine various scale factors between world space and image space:
	pixelsPerMetre = OUT_VIDEO_WIDTH / dimw;
	cellPixelWidth = (int)(cellsize * pixelsPerMetre);
	outVideoHeight = (int)(OUT_VIDEO_WIDTH * cellh / cellw);

	// Record timing info. (There must be an integer number of time steps per frame (1/25th of a second)
	duration = duration0;
	timestep = 0.04 / round(0.04 / timestep0);

	// Print values for confirmation:
	std::cout << "Fluid." << std::endl;
	std::cout << "\tDimensions: " << dimw << " " << dimh << std::endl;
	std::cout << "\tCellsize: " << cellsize << std::endl;
	std::cout << "\tCells: " << cellw << " " << cellh << std::endl;
	std::cout << "\tDuration: " << duration << " seconds. Timestep: " << timestep << "\n" << std::endl;

	// Copy pointers to solid object & fluid source data:
	solidObjects = solidObjects0;
	fluidSources = fluidSources0;

	// Given all the input data, set the system up:
	setup();
}

//
//	simulate() - runs the simulation as a sequence of iterations
//
void simulate()
{
	int numIterations = (int)(duration / timestep);
	int iterationsPerFrame = (int)(0.04 / timestep);

	videoOut_initialize(&hexBuffer, OUT_VIDEO_WIDTH, outVideoHeight, duration);

	for (int it = 0; it < numIterations; it++)
	{
		std::cout << "Time: " << it * timestep << "s. (iteration = " << it << ") ";

		// Enforce velocities at fluid sources:
		processFluidSources(it * timestep);

		// Assign the correct state to each cell based on current fluid particle positions:
		determineCellStates();

		// Compute all velocities for this time step
		fluidSolver_iterate(0, 0, 1, false, it);

		// Advect the particles using the latest velocity field:
		iterateParticles();

		// Once the next 1/25th of a second has been reached, draw the fluid in its current state as a frame of video:
		if (it % iterationsPerFrame == 0)
		{
			std::cout << "(Output frame)";

			visualiseFluid(it);
			videoOut_addFrame();
		}
		std::cout << std::endl;
	}

	videoOut_finalize();
}

//
//	cleanup() - called externally after simulation. Frees dynamically allocated objects
//
void cleanup()
{
	// Delete each cell's vector of particles:
	for (int i = 0; i < cellw*cellh; i++) delete particles[i];

	// Delete the array of particle vectors:
	free(particles);

	fluidSolver_cleanup();
}


//
//	randFloat() - returns a random number between 0 and 1
//
float randFloat()
{
	std::uniform_real_distribution<float> distribution(0.0, 1.0);
	return distribution(generator);
}

//
//	stage 2 of initialisation. Takes the fluid source & solid objects and applie their values to the simulation grid
//
void setup()
{
	// Create the solver object:
	fluidSolver_init(cellw, cellh, cellsize, 0.001, 0.0, 0.0, timestep);

	// Use the solid object data to determine which cells are solid:
	for (int s = 0; s < solidObjects->size(); s++)
	{
		// Determine the range of fluid cells whose centres are contained within this solid object:
		determineCellCoverage(&(solidObjects->at(s)));

		// Set these cells to solid in the solver:
		for (int i = solidObjects->at(s).i0; i <= solidObjects->at(s).i1; i++)
		{
			for (int j = solidObjects->at(s).j0; j <= solidObjects->at(s).j1; j++) setState(i, j, SOLID_CELL);
		}
	}

	// Determine and record which fluid cells are within each fluid source instance:
	for (int s = 0; s < fluidSources->size(); s++) determineCellCoverage(&(fluidSources->at(s).shape));

	// Particles are stored in lists, one for each cell:
	particles = (std::vector<particle> **)(malloc(cellw * cellh * sizeof(std::vector<particle> *)));
	for (int i = 0; i < cellw*cellh; i++) particles[i] = new std::vector<particle>;
}

//
//	determineCellCoverage() - Computes which cells are within the boundaries of a rectangle
//
void determineCellCoverage(rect *shape)
{
	shape->i0 = std::clamp((int)((shape->x / cellsize) + 0.5), 0, cellw - 1);
	shape->i1 = std::clamp((int)((shape->x / cellsize) + 0.5 + (int)(shape->w / cellsize) - 1), 0, cellw - 1);
	shape->j0 = std::clamp((int)((shape->y / cellsize) + 0.5), 0, cellh - 1);
	shape->j1 = std::clamp((int)((shape->y / cellsize) + 0.5 + (int)(shape->h / cellsize) - 1), 0, cellh - 1);
}

//
//	processFluidSources() - Applies fluid source properties in the form of velocities on the simulation grid
//
void processFluidSources(float t)
{
	// For each fluid source object, make sure its cells are filled wth the correct density of particles, and that its fixed velocities are applied:
	vector2D particlePos;
	vector2D newPos;

	for (int s = 0; s < fluidSources->size(); s++)
	{
		fluidSource fs = fluidSources->at(s);

		// Fluid sources only operate within a set time period:
		if ((t < fs.t0) || (t > fs.t1))
		{
			// Make sure this fluid source does not fix any cell velocities outside of its set time period:
			for (int i = fs.shape.i0; i <= fs.shape.i1; i++)
			{
				for (int j = fs.shape.j0; j <= fs.shape.j1; j++) setFixed(i, j, false);
			}

			continue;
		}

		for (int i = fs.shape.i0; i <= fs.shape.i1; i++)
		{
			for (int j = fs.shape.j0; j <= fs.shape.j1; j++)
			{
				// Remove all particles from this source cell:
				particles[j * cellw + i]->clear();

				// Fill the cell with a fresh set of uniformly distributed particles:
				for (int p = 0; p < PARTICLES_PER_CELL; p++)
				{
					particle pNew = { ((float)(i)+randFloat()) * cellsize, ((float)(j)+randFloat()) * cellsize , false };

					particles[j * cellw + i]->push_back(pNew);
				}
				
				// Set the velocity values specified by the user (if applicable), and fix them:
				if (fs.hasVel)
				{
					setVel(i, j, 0, fs.vx);
					setVel(i, j, 2, fs.vx);
					setVel(i, j, 1, fs.vy);
					setVel(i, j, 3, fs.vy); 
					setFixed(i, j, true);
				}
			}
		}
	}
}

//
//	determineCellStates() - Sets the state of every non-solid cell to air, surface fluid or full fluid based on which cells contain particles:
//
void determineCellStates()
{
	// Set all cells containing particles to be fluid cells, and all others to be air cells:
	for (int i = 0; i < cellw; i++)
	{
		for (int j = 0; j < cellh; j++)
		{
			// Solid cells are unchanged:
			if(getState(i, j) == SOLID_CELL) continue;

			// Set cell to be an air cell by default:
			setState(i, j, AIR_CELL);

			// If cell has any particles in it, it is a fluid cell:
			if (particles[j * cellw + i]->size() > 0) setState(i, j, FULL_FLUID_CELL);
		}
	}

	// Loop over all full fluid cells. Any that is adjacent to at least one air cell, is actually a surface fluid cell:
	for (int i = 0; i < cellw; i++)
	{
		for (int j = 0; j < cellh; j++)
		{
			if (getState(i, j) != FULL_FLUID_CELL) continue;

			// If any of this cell's 4 neighbour cells are an air cell, then this cell is a surface fluid cell:
			for (int n = 0; n < 4; n++)
			{
				int ni = i + neighbourIndex_i[n];
				int nj = j + neighbourIndex_j[n];

				if(getState(ni, nj) == AIR_CELL)
				{
					setState(i, j, SURFACE_FLUID_CELL);
					break;
				}
			}
		}
	}
}

//
//	iterateParticles() - Advects all fluid particles using the velocity field
//
void iterateParticles()
{
	// Split the cells into columns for multithreading:
	threadInfo threadInfos[NUM_THREADS];
	CWinThread *threadObj[NUM_THREADS];
	for (int t = 0; t < NUM_THREADS; t++)
	{
		threadInfos[t].id = t;
		threadInfos[t].start = (int)(t * (float)(cellw) / (float)(NUM_THREADS));
		threadInfos[t].end = (int)((t + 1) * (float)(cellw) / (float)(NUM_THREADS)-1);
	}

	// Advect the particles using multiple threads:
	for (int t = 0; t < NUM_THREADS; t++) threadObj[t] = AfxBeginThread(advectParticles_thread, &threadInfos[t]);
	for (int t = 0; t < NUM_THREADS; t++) WaitForSingleObject(threadObj[t]->m_hThread, INFINITE);

	// Corrective step: Iterate over the particles again, removing any that have drifted out of the range, and correcting the membership of any that have moved cell:
	for (int i = 0; i < cellw; i++)
	{
		for (int j = 0; j < cellh; j++)
		{
			for (auto it = particles[j * cellw + i]->begin(); it != particles[j * cellw + i]->end();)
			{
				// If particle is now outside the simulation area, delete:
				if ((it->x < 0) || (it->x >= dimw) || (it->y < 0) || (it->y >= dimh) || (std::isnan(it->x)) || (std::isnan(it->y)))
				{
					it = particles[j * cellw + i]->erase(it);
					continue;
				}

				// If particle is now in a different cell, remove it from this list and add it the new cell's list
				int i1, j1;
				i1 = (int)(it->x / cellsize);
				j1 = (int)(it->y / cellsize);

				if ((i1 != i) || (j1 != j))
				{
					particle pNew = { it->x, it->y, it->draw };
					it = particles[j * cellw + i]->erase(it);

					// If the particle has drifted into a solid cell, just remove it:
					if(getState(i1, j1) == SOLID_CELL) continue;
					particles[j1 * cellw + i1]->push_back(pNew);
					continue;
				}
				it++;
			}
		}
	}
}

//
//	advectParticles_thread() - per-thread function for advecting subsets of the fluid particles:
//
UINT advectParticles_thread(LPVOID pParam)
{
	vector2D displacement;

	threadInfo *info = (threadInfo *)(pParam);

	// Loop over all of this thread's cells:
	for(int i = info->start; i <= info->end; i++)
	{
		for (int j = 0; j < cellh; j++)
		{
			// Advect each particle using the fluid velocity at its position:
			for (int p = 0; p < particles[j * cellw + i]->size(); p++)
			{
				particle *p0 = &(particles[j * cellw + i]->at(p));
				displacement = getVel(p0->x, p0->y) * timestep;
				p0->x += displacement.x;
				p0->y += displacement.y;
			}
		}
	}

	return 0;
}


//
//	visualiseFluid() - Draws the current fluid particles and all other data to the frame buffer
//
void visualiseFluid(int it)
{
	// Draw the fluid elements using multiple threads:
	for (int t = 0; t < NUM_THREADS; t++) threadObj[t] = AfxBeginThread(visualiseFluid_thread, &threadInfos[t]);
	for (int t = 0; t < NUM_THREADS; t++) WaitForSingleObject(threadObj[t]->m_hThread, INFINITE);

	// Display the interation number:
	if(DRAW_ITERATION_NUMBER) drawIterationNumber(it);
}

//
//	visualiseFluid_thread() - each thread draws the fluid data for a subsection of the domain:
//
UINT visualiseFluid_thread(LPVOID pParam)
{
	// Check which part of the domain this thread is responsible for drawing:
	threadInfo *info = (threadInfo *)(pParam);
	int segStart = (int)(OUT_VIDEO_WIDTH * (info->start / (float)(cellw)));
	int segEnd = (int)(OUT_VIDEO_WIDTH * ((info->end + 1) / (float)(cellw))) - 1;

	// Set background to white:
	for (int i = segStart; i <= segEnd; i++)
	{
		for (int j = 0; j < outVideoHeight; j++)
		{
			for (int c = 0; c < 4; c++) hexBuffer[(j * OUT_VIDEO_WIDTH + i) * 4 + c] = (char)(255);
		}
	}

	// Draw the cell grid lines & fill with colour based on state:
	if (DRAW_CELLS)
	{
		float r, g, b;

		for (int i = info->start; i <= info->end; i++)
		{
			for (int j = 0; j < cellh; j++)
			{
				if (getState(i, j) == AIR_CELL) continue;													// Air cells are already coloured white:
				if (getState(i, j) == SURFACE_FLUID_CELL) r = g = b = 0.125;								// Fluid surface cells in dark grey:											
				if (getState(i, j) == FULL_FLUID_CELL) r = g = b = 0.25;									// Full fluid cells in light grey

				for (int c = i * cellPixelWidth; c < (i + 1) * cellPixelWidth; c++)
				{
					for (int d = j * cellPixelWidth; d < (j + 1) * cellPixelWidth; d++)
					{
						hexBuffer[((outVideoHeight - 1 - d) * OUT_VIDEO_WIDTH + c) * 4 + 2] = (char)(255.0 * r);
						hexBuffer[((outVideoHeight - 1 - d) * OUT_VIDEO_WIDTH + c) * 4 + 1] = (char)(255.0 * g);
						hexBuffer[((outVideoHeight - 1 - d) * OUT_VIDEO_WIDTH + c) * 4 + 0] = (char)(255.0 * b);
					}
				}
			}
		}
		
		// Draw the grid lines in grey: 
		for (int i = info->start; i <= info->end; i++)
		{
			int li = i * cellPixelWidth;

			for (int j = 0; j < outVideoHeight; j++)
			{
				for(int c = 0; c < 3; c++)	hexBuffer[((outVideoHeight - 1 - j) * OUT_VIDEO_WIDTH + li) * 4 + c] = (char)(65);
			}
		}
		for (int j = 0; j < cellh; j++)
		{
			int lj = j * cellPixelWidth;

			for (int i = segStart; i <= segEnd; i++)
			{
				for (int c = 0; c < 3; c++) hexBuffer[((outVideoHeight - 1 - lj) * OUT_VIDEO_WIDTH + i) * 4 + c] = (char)(65);
			}
		}
	}

	// Draw the solid objects dark grey:
	if (DRAW_SOLIDS)
	{
		for (int i = info->start; i <= info->end; i++)
		{
			for (int j = 0; j < cellh; j++)
			{
				if (getState(i, j) != SOLID_CELL) continue;										
				for (int pi = i * cellPixelWidth; pi < (i + 1) * cellPixelWidth; pi++)
				{
					for (int pj = j * cellPixelWidth; pj < (j + 1) * cellPixelWidth; pj++)
					{
						for(int c = 0; c < 3; c++) hexBuffer[((outVideoHeight - 1 - pj) * OUT_VIDEO_WIDTH + pi) * 4 + c] = (char)(65);
					}
				}
			}
		}
	}

	// Draw the velocity componenents in yellow:
	if(DRAW_VELOCITIES)
	{
		int Pi, Pj;
		float val;
		int lineLength;

		// Draw the horizontal velocity components:
		for (int i = info->start; i <= (info->end + 1); i++)
		{
			// Only the last thread needs to draw one extra column of velocities on the extreme right:
			if((info->id < (NUM_THREADS - 1))&&(i == (info->end + 1))) continue;

			for (int j = 0; j < cellh; j++)
			{
				val = getVelH(i, j);
				Pj = (int)(((float)(j)+0.5) * (float)(cellPixelWidth));

				// start point of the velocity line:
				if(val > 0.0) Pi = i * cellPixelWidth;
				else Pi = (int)(i * cellPixelWidth + VELOCITY_DRAW_SCALE * val * pixelsPerMetre);

				lineLength = (int)(fabs(VELOCITY_DRAW_SCALE * val * pixelsPerMetre));

				drawVelocityLine(Pi, Pj, 1, 0, lineLength);
			}
		}

		// Draw the vertical velocity components:
		for (int i = info->start; i <= info->end; i++)
		{
			for (int j = 0; j <= cellh; j++)
			{
				val = getVelV(i, j);
				Pi = (int)(((float)(i)+0.5) * cellPixelWidth);

				// start point of the velocity line:
				if (val > 0.0) Pj = j * cellPixelWidth;
				else Pj = (int)(j * cellPixelWidth + VELOCITY_DRAW_SCALE * val * pixelsPerMetre);

				lineLength = (int)(fabs(VELOCITY_DRAW_SCALE * val * pixelsPerMetre));

				drawVelocityLine(Pi, Pj, 0, 1, lineLength);
			}
		}
	}
	
	// Draw the particles:
	if (DRAW_PARTICLES)
	{
		vector2D particlePos;
		int pp_i, pp_j;

		for (int i = info->start; i <= info->end; i++)
		{
			for (int j = 0; j < cellh; j++)
			{
				for (auto it = particles[j * cellw + i]->begin(); it != particles[j * cellw + i]->end(); it++)
				{
					pp_i = (int)(OUT_VIDEO_WIDTH * it->x / dimw);
					pp_j = (int)(outVideoHeight * it->y / dimh);

					for (int pp_i1 = pp_i - 1; pp_i1 <= (pp_i + 1); pp_i1++)
					{
						for (int pp_j1 = pp_j - 1; pp_j1 <= (pp_j + 1); pp_j1++)
						{
							if ((pp_i1 < 0) || (pp_i1 >= OUT_VIDEO_WIDTH) || (pp_j1 < 0) || (pp_j1 >= outVideoHeight)) continue;

							hexBuffer[((outVideoHeight - 1 - pp_j1) * OUT_VIDEO_WIDTH + pp_i1) * 4 + 2] = (char)(0);
							hexBuffer[((outVideoHeight - 1 - pp_j1) * OUT_VIDEO_WIDTH + pp_i1) * 4 + 1] = (char)(0);
							hexBuffer[((outVideoHeight - 1 - pp_j1) * OUT_VIDEO_WIDTH + pp_i1) * 4 + 0] = (char)(255);
						}
					}
				}
			}
		}
	}
	
	return 0;
}

//
// drawVelocityLine() - given a start point (P), direction vector (D) and length, all in pixel space, draws a line in yellow:
//
void drawVelocityLine(int Pi, int Pj, int Di, int Dj, int lineLength)
{
	for(int a = 0; a < lineLength; a++)
	{
		int i = Pi + Di * a;
		int j = Pj + Dj * a;

		if ((i < 0) || (i >= OUT_VIDEO_WIDTH)||(j < 0)||(j >= outVideoHeight)) continue;

		for (int c = 0; c < 3; c++) hexBuffer[((outVideoHeight - 1 - j) * OUT_VIDEO_WIDTH + i) * 4 + c] = yellow[c];
	}
}

//
//	drawIterationNumber() - draws the current iteration count in the top left corner (4 digits only)
//
void drawIterationNumber(int it)
{
	int numDigits;
	// Determine the number of digits from the order of magnitude:
	for (int oom = 1; oom < 10; oom++)
	{
		if (it < oom * 10) { numDigits = oom; break; }
	}

	int digit[10];

	for (int oom = 0; oom < numDigits; oom++)
	{
		float it0 = (float)(it);
		for (int a = 0; a < oom; a++) it0 -= digit[a] * pow(10.0, numDigits - 1 - a);
		digit[oom] = (int)(it0 / pow(10.0, numDigits - 1 - oom));
	}

	// Draw the frame number:
	for (int d = 0; d < numDigits; d++)
	{
		int a = 2 + d * 4;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				for (int i1 = i * 4; i1 < (i + 1) * 4; i1++)
				{
					for (int j1 = j * 4; j1 <= (j + 1) * 4; j1++)
					{
						if (digits[digit[d]][j * 3 + i]) hexBuffer[((50 - j1) * OUT_VIDEO_WIDTH + i1 + a * 4) * 4 + 1] = (char)(255);
					}
				}
			}
		}
	}
}

