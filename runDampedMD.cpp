#include <vector>
#include <cmath>
#include <numeric>
#include "MD.h"


// Run dampedMD()
// Uses InteractionManager::force_threshold and InteractionManager::central_force_mag


void MD::runDampedMD()
{

	int step = 1;

	do
	{

		// Run MD step
		dampedMD();

		// Check to see if Verlet list should be remade
		makeVerletList();

		// Write Files
		if(step % writestep==0)
		{
			// Compute extra quantities
			computeTemp();
			sim->computeRg();
			interman->Etot = interman->KE + interman->PE;
			adjust_dt();

			// Write files for NVE steps only
			writeFiles(step);

		}

		// increment step
		step++;
	}
	// when relaxing to remove overlaps:
	while(interman->Fmag_max > interman->Fthreshold || step < 50000);




	// Write final configuration to file
	computeTemp();
	sim->computeRg();
	interman->Etot = interman->KE + interman->PE;

	writeFiles(step+1, 1, 1);

}


