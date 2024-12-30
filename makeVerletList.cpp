#include <cmath>
#include <numeric>
#include "MD.h"
#include "Vector3D.h"



// This is used to keep track of nonbonded_pairs in Simulation.
// All other bonds are permanent, e.g., backbone_pairs, bond_angles, bond_dihedrals.
// Parameter bool initialize is used to make the initial Verlet list and bypass the displacement check


void MD::makeVerletList(bool initialize)
{

	// Check the displacements of all spheres
	// If the sum of the two max displacements is greater than the verlet_skin, 
	// remake the Verlet list
	

	double max_disp=0.0, max_disp_2=0.0;
	double disp_mag;
	for(const auto& sph : sim->spheres)
	{
		const auto& displacement = sph->displacement;
		disp_mag = displacement.norm();
		disp_mag > max_disp_2 ? disp_mag > max_disp ? max_disp = disp_mag : max_disp_2 = disp_mag : 0;
	}


	if(max_disp + max_disp_2 >= verlet_skin - 1.0 || initialize)
	{
		// Clear current list nonbonded_pairs
		sim->nonbonded_pairs.clear();

		// Zero the displacement vectors
		for(const auto& sph : sim->spheres)
		{
			sph->displacement.zeroVector();
		}


		// Remake nonbonded_pairs in Simulation.h
		// Iterate through all pairs of spheres that are not permanently bonded
		// These are all the bonds that may be temporarily interacting
		std::unordered_map<int, std::shared_ptr<Sphere> > sphereMap;
		int sphaID, sphbID;
		bool spha_found, sphb_found;
		std::shared_ptr<Sphere> spha_ptr;
		std::shared_ptr<Sphere> sphb_ptr;
		//double distsq;
		double dist;


		for(int bond_ind=0; bond_ind<sim->num_non_perm_bond_pairs; ++bond_ind)
		{
			const auto& bond = sim->all_non_perm_bond_pairs[bond_ind];

			sphaID = bond->sphereIDs[0];
			sphbID = bond->sphereIDs[1];

			sphereMap.clear();

			// If atomIDs match the index in the spheres vector, it's easy
			if(sim->Natoms>=sphaID && sim->Natoms>=sphbID && sim->spheres[sphaID]->atomID==sphaID && sim->spheres[sphbID]->atomID==sphbID)
			{
				auto* spha = sim->spheres[sphaID].get();
				auto* sphb = sim->spheres[sphbID].get();


    			// Get distance between spheres
				dist = spha->distij(sphb, 0);
			
				// bond length for two particles that are not bonded is the sum of their radii
				if(dist/bond->bond_length <= verlet_skin)
				{
					sim->nonbonded_pairs.push_back(bond_ind);
				}

			}	

			// otherwise we need to find the correct sphere object...
			else
			{
				spha_found = 0;
				sphb_found = 0;
				for(const auto& sph : sim->spheres)
				{					
					if(sph->atomID==sphaID)
					{
						spha_ptr = sph;
						spha_found = 1;
					}
					else if(sph->atomID==sphbID)
					{
						sphb_ptr = sph;
						sphb_found = 1;
					}
					if(spha_found && sphb_found)
					{
						auto* spha = spha_ptr.get();
						auto* sphb = sphb_ptr.get();
						sphereMap[sphaID] = spha_ptr;
						sphereMap[sphbID] = sphb_ptr;

						// Get distance between spheres
						dist = spha->distij(sphb, 0);						

						// bond length for two particles that are not bonded is the sum of their radii
						if(dist/bond->bond_length <= verlet_skin)
						{
							sim->nonbonded_pairs.push_back(bond_ind);
						}
						break;
					}
				}

			}

		}		

		// Update the verlet skin
		// verlet_skin_mult is defined in MD.h and can be changed
		// if it is >1 then the verlet list is updated less often.
		// if it is <1 it is safer and the verlet list is updated before the spheres move 
		// as much as they moved the previous amounte between verlet list updates

		//verlet_skin = (max_disp + max_disp_2) * verlet_skin_mult;

		// smoother transitions

		//verlet_skin = 0.9*verlet_skin + 0.1*(max_disp + max_disp_2)*verlet_skin_mult;
		//std::cout<<"verlet_skin: "<<verlet_skin<<"\n";
	}

}




