#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <string>
#include <algorithm>

#include "datatypes.h"
#include "liquidSimulator.h"
#include "fluidIO.h"

int main(int argc, char *argv[])
{
	// Process arguments:
	if((argc < 2)||(argc > 3))
	{
		std::cout << "Error. Wrong number of parameters. Useage: \n\n>fluidSim.exe <fluidFile.fluid> or >fluidSim.exe <fluidFile.fluid> -debug" << std::endl;
		return 0;
	}

	std::string filename(argv[1]);
	std::string arg2 = "";

	if(argc == 3) arg2 = argv[2];
	

	// Open & parse fluid file:
	std::ifstream file(filename);
	std::string fileword;

	float dimw, dimh;														// Width and height of rectangular fluid system in 2D space
	float cellsize;															// The width/height of the square cells
	float duration, timestep;												// The duration of the sim and the length of each time step
	std::vector<rect> solidObjects;											// Container for all solid objects defined by user
	std::vector<fluidSource> fluidSources;									// Container for all fluid sources defined by user

	// Load all the information from the fluid file:
	loadFluidProperties(file, &dimw, &dimh, &cellsize, &duration, &timestep, &solidObjects, &fluidSources);
	
	// Create a liquid simulator and pass all this data to it:
	initialise(dimw, dimh, cellsize, duration, timestep, &solidObjects, &fluidSources);

	// Check whether the user specified debug mode:
	if (argc == 3)
	{
		std::string arg2(argv[2]);
		if (arg2 == "-debug")
		{
			std::cout << "DEBUG MODE" << std::endl;
		}
	}

	// Run the simulation:
	simulate();

	cleanup();

	return 0;
}
