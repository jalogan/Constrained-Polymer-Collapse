#ifndef MD_H
#define MD_H


#include <vector>
#include "Simulation.h"
#include "InteractionManager.h"



class MD
{
	public:

		Simulation* sim;
		InteractionManager* interman;
		std::string IN;
		std::string input_filename;
		std::string OUT;
		double deltat;
		double damping;
		int writestep;		

		double temperature;

		static constexpr double verlet_skin = 1.2;
		static constexpr double verlet_skin_mult = 1.0;
		static constexpr double mindt = 1e-5;
		static constexpr double maxdt = 1.0;


		MD(Simulation* simulation, InteractionManager* interaction_manager, std::string& inpath, std::string& infile, std::string& outpath, double dt, double damp_param, int writestep_, double initial_temp) : sim{simulation}, interman{interaction_manager}, IN{inpath}, input_filename{infile}, OUT{outpath}, deltat{dt}, damping{damp_param}, writestep{writestep_}, temperature{initial_temp} {}


		void initializeTemp();
		void makeVerletList(bool initialize=0);
		void runNVE(int numsteps, bool write=1);
		void NVE_MD();
		void runNVT(double desired_temp, int numsteps, bool write=1);
		void NVT_MD(double desired_temp);
		void runDampedMD();
		void dampedMD();
		void adjust_dt();
		void writeFiles(int step, bool overwrite_files=1, bool write_vis_file=0);
		void loadConfig();



		void computeKE()
		{
			interman->KE = 0.0;

			for(const auto& sph : sim->spheres)
			{
				const auto& vel = sph->velocity;
				interman->KE += 0.5*sph->mass*(vel.dot(vel));
			}
		}


		void computeTemp()
		{
			computeKE();
			temperature = (2.0*interman->KE) / (sim->dim * sim->Natoms);
		}


};






#endif

