#include <cstdio>
#include <iostream>
#include <string>
#include "rule.h"
#include "monte_carlo.h"
#include "session.h"
#include "mc_utils.h"
using namespace std;
// The MC code
// The paramiters for the MC run ( Composition, Simulation cell size, ect...) are all set in INPUT.txt
// These settings are read in and stored in the SimCell object. The SimCell object is then reffrenced to set up the actual MC run.
SimCell sim_cell;
int main(void) {
	Session ses("INPUT"); // settings for MC run from input file
	if (ses.use_states == true) { ses.get_spin_states("SPIN_STATES"); }
	vector<Rule> mc_rules; // Initalize the "Rules" that will be used to properly apply the ECI to each site in the simulation cell
	fillRuleList(mc_rules, "Rule_file", 0); // Populate mc_rules. Rules are defined as follows, [species (0,1,2 for Ni,Mn,In)]:[NA if monomer,dist of dymer, AB,AC,BC if trimer]:no-spin(0) or spin(1): value
	vector<float> dist_list{ 0.0 }; // Initalize list of discances used in mc_rules
	fillDistList(dist_list, mc_rules); // Fill the dist_list 
	cout << "Dist List filled\n";
	// Create the simulation cell object. Arguments (POSCAR_file, dist_list, shape, species numbs, cutoff (currently unused), sim_type (also unused), phase_init (aust/mart), spin_init (AFM/FM/RAND), species_init (Ordered/Random), bool use_poscar) 
	sim_cell.initSimCell("POSCAR", dist_list, ses); // Create and initalize the simulation cell
	MCrun mcrun(ses, sim_cell);
	mcrun.start(mc_rules);
	return 0;