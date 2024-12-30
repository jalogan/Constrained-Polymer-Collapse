#include "InteractionManager.h"
#include "Vector3D.h"


// Include any external forces, e.g., central force, on each sphere


void InteractionManager::computeExternalInteractions()
{

	// Include central force on all spheres
	// The central force acts differently depending on the size of the sphere
	if(central_force_mag > 0.0)
	{
		double scale;
		Vector3D unitba;
		Vector3D centralF;

		for(auto& sph_ptr : sim->spheres)
		{	
			auto* sph = sph_ptr.get();

			// unit vector from particle to origin (assumed to be approximate COM)
			unitba = -sph->position.normalized();

			scale = central_force_mag * pow(sph->diameter / sim->max_diameter, 2.25);
			centralF = unitba*scale;

			sph->force += centralF;	

		}
	}


}


