#include <cstdio>
#include <iostream>
#include <string>
#include "run_sim.h"
#include "session.h"
#include "utils.h"
using namespace std;
// The MC code
// The paramiters for the MC run ( Composition, Simulation cell size, ect...) are all set in INPUT.txt
// These settings are read in and stored in the SimCell object. The SimCell object is then reffrenced to set up the actual MC run.
SimCell sim_cell;
int main(void) {
	Session ses("INPUT"); // settings for MC run from input file
	if (ses.use_states == true) { ses.add_spin_states("SPIN_STATES"); }
	ses.fill_rule_list(); // Initalize the "Rules" that will be used to properly apply the ECI to each site in the simulation cell. Rules are defined as follows, [species (0,1,2 for Ni,Mn,In)]:[NA if monomer,dist of dymer, AB,AC,BC if trimer]:no-spin(0) or spin(1): value
	ses.find_unique_dists(); // Find the list of discances used in mc_rules
	cout << "Unique Dist List filled\n";
	// Create the simulation cell object. Arguments (POSCAR_file, dist_list, shape, species numbs, cutoff (currently unused), sim_type (also unused), phase_init (aust/mart), spin_init (AFM/FM/RAND), species_init (Ordered/Random), bool use_poscar) 
	sim_cell.initSimCell("POSCAR", ses); // Create and initalize the simulation cell
	Sim sim(ses, sim_cell);
	sim.start();
	return 0;
}
