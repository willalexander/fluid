
#include "fluidIO.h"

//
//	loadFluidProperties() - Parses a '.fluid' file and passes the information to the caller by reference
//
void loadFluidProperties(std::ifstream &file, float *dimw, float *dimh, float *cellsize, float *duration, float *timestep, std::vector<rect> *solidObjects, std::vector<fluidSource> *fluidSources)
{
	std::string fileword;

	while (file >> fileword)
	{
		// Basic fluid dimensions. Width and height in 2D space, and number of discrete cells in width are given by user. Num cells in height is determined from this: 
		if (fileword == "fluid")
		{
			file >> fileword; *dimw = atof(fileword.c_str());
			file >> fileword; *dimh = atof(fileword.c_str());
			file >> fileword; *cellsize = atof(fileword.c_str());

			continue;
		}

		// Basic timing information: Duration of simulation, and length of each time step (in seconds)
		if (fileword == "time")
		{
			file >> fileword; *duration = atof(fileword.c_str());
			file >> fileword; *timestep = atof(fileword.c_str());
			continue;
		}

		// Load any solid rectangular obstacles that the user wants in the fluid:
		if (fileword == "obj")
		{
			rect newSolidObject;
			file >> fileword; newSolidObject.x = atof(fileword.c_str());
			file >> fileword; newSolidObject.y = atof(fileword.c_str());
			file >> fileword; newSolidObject.w = atof(fileword.c_str());
			file >> fileword; newSolidObject.h = atof(fileword.c_str());
			solidObjects->push_back(newSolidObject);

			std::cout << "Solid object. " << newSolidObject.x << " " << newSolidObject.y << " " << newSolidObject.w << " " << newSolidObject.h << " " << std::endl;

			continue;
		}

		// Load any fluid sources. (These are required for any fluid to be introduced to the system. They are rectangular and have a velocity associated with them:
		if (fileword == "src")
		{
			fluidSource newFluidSource;
			file >> fileword; newFluidSource.shape.x = atof(fileword.c_str());
			file >> fileword; newFluidSource.shape.y = atof(fileword.c_str());
			file >> fileword; newFluidSource.shape.w = atof(fileword.c_str());
			file >> fileword; newFluidSource.shape.h = atof(fileword.c_str());
			file >> fileword; newFluidSource.t0 = atof(fileword.c_str());
			file >> fileword; newFluidSource.t1 = atof(fileword.c_str());

			// Applying an initial velocity by the fluid source is optional, so test to see whether a velocity value is provided:
			fileword = "";
			file >> fileword; newFluidSource.vx = atof(fileword.c_str());
			if (fileword == "") newFluidSource.hasVel = false;
			else
			{
				file >> fileword; newFluidSource.vy = atof(fileword.c_str());
				newFluidSource.hasVel = true;
			}

			fluidSources->push_back(newFluidSource);

			std::cout << "Fluid source. " << std::endl;
			std::cout << "\tOrigin: " << newFluidSource.shape.x << ", " << newFluidSource.shape.y << std::endl;
			std::cout << "\tDimensions: " << newFluidSource.shape.w << ", " << newFluidSource.shape.h << std::endl;
			std::cout << "\tTime period: " << newFluidSource.t0 << ", " << newFluidSource.t1 << std::endl;
			if (newFluidSource.hasVel) std::cout << "\tVelocity: " << newFluidSource.vx << ", " << newFluidSource.vy << std::endl;
			std::cout << "\n" << std::endl;

			continue;
		}
	}
}