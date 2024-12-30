#include "MD.h"


// Run NVE_MD()


void MD::runNVE(int numsteps, bool write)
{

	for(int step = 1; step < 1 + numsteps; ++step)
	{
		// Run MD step
		NVE_MD();

		// Check if Verlet should be remade
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
			if(write)
			{
				writeFiles(step);
			}
		}
	}


}


