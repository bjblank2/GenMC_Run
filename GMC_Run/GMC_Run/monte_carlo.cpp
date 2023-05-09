#include "monte_carlo.h"

extern SimCell sim_cell;
extern Session session;

MCrun::MCrun(void) {}

MCrun::MCrun(Session& _session, SimCell& _sim_cell) {
	sim_cell = _sim_cell;
	session = _session;
}

void MCrun::start() {
	cout <<"using algo " << session.algo << "\n";
	if (session.algo == 0) { return; }
	else if (session.algo == 3) { runMetropolis3(); }
	else if (session.algo == 4) { runMetropolis4(); }
	else if (session.algo == 5) { runMetropolis5(); }
	else if (session.algo == 6) { runMetropolis6(); }
	else if (session.algo == 7) { runMetropolis7(); }
	else if (session.algo == 8) { runMetropolis8(); }
	else if (session.algo == 9) { runMetropolis9(); }
	else if (session.algo == 10) { runMetropolis10(); }
	else if (session.algo == -1) { outputSimCell(); }
	else if (session.algo == -2) { debug_run(); }
}
// Evaluate the energy controbution from a singal atom, includeing spin chem rules.
long double MCrun::evalSiteEnergyAll(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	// Initalize variables
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	int dist_2_ind;
	float neighbor_dist2;
	float neighbor_dist3;
	key = "_" + to_string(site_species) + ",0,"; // make key for one atom term
	rule_itr = rule_map_chem.find(key);
	enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule for one atom term
	// for 2 atom terms
	// loop through all neighbors in the neighbor list (list only includes atoms that have distances listed in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		if (neighbor_site1 == site) { continue; }
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // make key string
		rule_itr = rule_map_chem.find(key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // if the key was in the chem rule map, apply the rule, otherwise do nothing
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if the key was in the spin rule map, apply the rule, otherwise do nothing
		enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
		// for 3 atom terms
		// loop through all neighbors in the neighbor list
		//for (int j = 0; j < neighbor_index_list[site].size(); j++) {
		//	neighbor_site2 = neighbor_index_list[site][j];
		//	if (neighbor_site1 == neighbor_site2 or neighbor_site2 == site) { continue; }
		//	neighbor_species2 = atom_species[neighbor_site2];
		//	neighbor_spin2 = atom_spin[neighbor_site2];
		//	vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
		//	if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
		//		dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
		//		neighbor_dist2 = neighbor_dist_list[site][j];
		//		neighbor_dist3 = neighbor_dist_list[neighbor_site1][dist_2_ind];
		//		// make key for three atom term
		//		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
		//		key +="," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
		//		rule_itr = rule_map_chem.find(key);
		//		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / 6 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
		//		rule_itr = rule_map_spin.find(key);
		//		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 * 1 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
		//	}
		//}
	}
	return enrg;
}
// Evaluate the energy contribution from a singal atom using only spin rules. (Assuming the spin is all that changed)
long double MCrun::evalSiteEnergySpin(int site, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	// Initalize variables and constants
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;

	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	float neighbor_dist2;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int dist_2_ind;
	float neighbor_dist3;
	enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
	// for 2 atom terms
	// Loop through all neighbors in the neighbor list (neighbor list only includes distances included in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		// gather info to build key
		neighbor_site1 = neighbor_index_list[site][i];
		if (neighbor_site1 == site) { continue; }
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // create key string
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if key is found in rule list, apply the rule, else do nothing.
		//for (int j = 0; j < neighbor_index_list[site].size(); j++) {
		//	neighbor_site2 = neighbor_index_list[site][j];
		//	if (neighbor_site1 == neighbor_site2 or neighbor_site2 == site) { continue; }
		//	neighbor_species2 = atom_species[neighbor_site2];
		//	neighbor_spin2 = atom_spin[neighbor_site2];
		//	vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
		//	if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
		//		dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
		//		neighbor_dist2 = neighbor_dist_list[site][j];
		//		neighbor_dist3 = neighbor_dist_list[neighbor_site1][dist_2_ind];
		//		// make key for three atom term
		//		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
		//		key += "," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
		//		rule_itr = rule_map_spin.find(key);
		//		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 * 1 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
		//	}
		//}
	}
	return enrg;
}
// Evaluate the energy controbution from a singal atom, includeing spin chem rules.
long double MCrun::evalSiteEnergyAll3(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	// Initalize variables
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	int dist_2_ind;
	float neighbor_dist2;
	float neighbor_dist3;
	enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
	key = "_" + to_string(site_species) + ",0,"; // make key for one atom term
	rule_itr = rule_map_chem.find(key);
	enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule for one atom term
	// for 2 atom terms
	// loop through all neighbors in the neighbor list (list only includes atoms that have distances listed in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		if (neighbor_site1 == site) { continue; }
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // make key string
		rule_itr = rule_map_chem.find(key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // if the key was in the chem rule map, apply the rule, otherwise do nothing
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if the key was in the spin rule map, apply the rule, otherwise do nothing
		// for 3 atom terms
		// loop through all neighbors in the neighbor list
		for (int j = 0; j < neighbor_index_list[site].size(); j++) {
			neighbor_site2 = neighbor_index_list[site][j];
			if (neighbor_site1 == neighbor_site2 or neighbor_site2 == site) { continue; }
			neighbor_species2 = atom_species[neighbor_site2];
			neighbor_spin2 = atom_spin[neighbor_site2];
			vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
			if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
				dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
				neighbor_dist2 = neighbor_dist_list[site][j];
				neighbor_dist3 = neighbor_dist_list[neighbor_site1][dist_2_ind];
				// make key for three atom term
				key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
				key += "," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
				rule_itr = rule_map_chem.find(key);
				enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
				rule_itr = rule_map_spin.find(key);
				enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
			}
		}
	}
	return enrg;
}
// Evaluate the energy contribution from a singal atom using only spin rules. (Assuming the spin is all that changed)
long double MCrun::evalSiteEnergySpin3(int site, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	// Initalize variables and constants
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;

	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	float neighbor_dist2;

	int neighbor_site3;
	int neighbor_species3;
	float neighbor_dist3;

	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int dist_2_ind;
	enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
	// for 2 atom terms
	// Loop through all neighbors in the neighbor list (neighbor list only includes distances included in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		// gather info to build key
		neighbor_site1 = neighbor_index_list[site][i];
		if (neighbor_site1 == site) { continue; }
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // create key string
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if key is found in rule list, apply the rule, else do nothing.
		for (int j = 0; j < neighbor_index_list[site].size(); j++) {
			neighbor_site2 = neighbor_index_list[site][j];
			if (neighbor_site1 == neighbor_site2 or neighbor_site2 == site) { continue; }
			neighbor_species2 = atom_species[neighbor_site2];
			neighbor_spin2 = atom_spin[neighbor_site2];
			vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
			if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
				dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
				neighbor_dist2 = neighbor_dist_list[site][j];
				neighbor_dist3 = neighbor_dist_list[neighbor_site1][dist_2_ind];
				// make key for three atom term
				key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
				key += "," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
				rule_itr = rule_map_spin.find(key);
				if (rule_itr != rule_map_spin.end()) {
					enrg += rule_itr->second * site_spin * neighbor_spin1;
					for (int k = 0; k < neighbor_index_list[neighbor_site1].size(); k++) {
						neighbor_site3 = neighbor_index_list[neighbor_site1][k];
						if (neighbor_site3 == site or neighbor_site1 == neighbor_site2) { continue; }
						neighbor_species3 = atom_species[neighbor_site3];
						vector<int>::iterator dist_3_itr = find(neighbor_index_list[site].begin(), neighbor_index_list[site].end(), neighbor_site3);
						if (dist_3_itr != neighbor_index_list[site].end()) {
							neighbor_dist2 = neighbor_dist_list[neighbor_site1][k];
							neighbor_dist3 = neighbor_dist_list[site][distance(neighbor_index_list[site].begin(), dist_3_itr)];
							key = "_" + to_string(neighbor_species1) + "," + to_string(site_species) + "," + to_string(neighbor_species3);
							key += "," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
							//cout << key << "\n";
							rule_itr = rule_map_spin.find(key);
							enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if key is found in rule list, apply the rule, else do nothing.
						}
					}
				}
				//enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 * 1 : 0.0; // apply rule (if possible). 	
			}
		}
	}
	return enrg;
}
// Evaluate the total energy of the simulation cell in it's current state
long double MCrun::evalLattice3(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	long double enrg = 0;
	// loop thruogh simulation cell and add everything up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergyAll3(site, rule_map_spin, rule_map_chem, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;
}
// Evaluate the spin contribution to the total energy of the simulation cell
long double MCrun::evalLatticeSpin3(map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	long double enrg = 0;
	// loop through everything and add it up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergySpin3(site, rule_map_spin, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;// / atom_species.size() * 16 - 66.8295069760671;
}

float MCrun::evalLatticeFast(vector<int> spins, vector<int> species, vector<vector<vector<int>>> JcountChem, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatChem, vector<vector<float>> JMatSpin) {
	return evalLattiveChemFast(JcountChem, JMatChem) + evalLatticeSpinFast(spins, species, JcountSpin, JMatSpin);
}

float MCrun::evalLatticeSpinFast(vector<int> spins, vector<int> species, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy = 0;
	for (int i = 0; i < JcountSpin.size(); i++) {
		site_species = species[i];
		energy += evalSiteEnergySpinFast(spins[i], spins, JcountSpin[i], JMatSpin);
	}
	return energy;
}

float MCrun::evalLattiveChemFast(vector<vector<vector<int>>> JcountChem, vector<vector<float>> JMatChem) {
	float energy = 0;
	for (int i = 0; i < JcountChem.size(); i++) {
		energy += evalSiteEnergyChemFast(JcountChem[i], JMatChem);
	}
	return energy;
}

float MCrun::evalSiteEnergyAllFast(int spin, vector<int> spins, vector<vector<int>> JcountChem, vector<vector<int>> JcountSpin, vector<vector<float>> JMatChem, vector<vector<float>> JMatSpin) {
	float energy = 0;
	for (int i = 0; i < JcountChem.size(); i++) { energy += JcountChem[i].size() * JMatChem[i][0]; }
	energy += evalSiteEnergySpinFast(spin, spins, JcountSpin, JMatSpin);
}

float MCrun::evalSiteEnergySpinFast(int spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy = 0;
	float energy_inc;
	for (int i = 0; i < JcountSpin.size(); i++) {
		energy_inc = 0;
		for (int j = 0; j < JcountSpin[i].size(); j++) {
			energy_inc += spins[JcountSpin[i][j]];
		}
		energy += energy_inc * JMatSpin[i][0];
	}
	//energy -= session.mag_ext * session.moments[site_species] * uB;
	return 0.5 * energy * spin;
}

float MCrun::evalSiteEnergyChemFast(vector<vector<int>> JcountChem, vector<vector<float>> JMatChem) {
	long double enrg = 0;
	for (int i = 0; i < JcountChem.size(); i++) { enrg += JcountChem[i].size() * JMatChem[i][0]; }
	return 0.5 * enrg;
}

float MCrun::evalDelEnergySpinFast(int old_spin, int new_spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy = 0;
	float energy_inc;
	int field_itr;
	bool field_term = false;
	for (int i = 0; i < JcountSpin.size(); i++) {
		if (JMatSpin[i][1] != 1) {
			energy_inc = 0;
			for (int j = 0; j < JcountSpin[i].size(); j++) {
				energy_inc += spins[JcountSpin[i][j]];
			}
			energy += energy_inc * JMatSpin[i][0];
		}
		else {
			field_term = true;
			field_itr = i;
		}
	}
	//energy -= session.mag_ext * session.moments[site_species] * uB;
	if (field_term == false) { return (0.5 * (energy * new_spin - energy * old_spin)); }
	else if (field_term == true) { return (0.5 * (energy * new_spin - energy * old_spin)) + (JMatSpin[field_itr][0] * pow(new_spin, 2) - JMatSpin[field_itr][0] * pow(old_spin, 2)); }

}

float MCrun::evalPottsLatticeSpinFast(vector<int> spins, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy = 0;
	for (int i = 0; i < JcountSpin.size(); i++) {
		energy += evalPottsSiteEnergySpinFast(spins[i], spins, JcountSpin[i], JMatSpin);
	}
	return energy;
}

float MCrun::evalPottsSiteEnergySpinFast(int spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy = 0;
	float energy_inc;
	for (int i = 0; i < JcountSpin.size(); i++) {
		energy_inc = 0;
		for (int j = 0; j < JcountSpin[i].size(); j++) {
			if (spins[JcountSpin[i][j]] == spin) { energy_inc += 1; }
		}
		energy += energy_inc * float(JMatSpin[i][0]);
	}
	return energy;
}

float MCrun::evalPottsDelEnergySpinFast(int old_spin, int new_spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin) {
	float energy_old = 0;
	float energy_new = 0;
	float energy_inc_old;
	float energy_inc_new;
	bool field_term = false;
	int field_itr;
	for (int i = 0; i < JcountSpin.size(); i++) {
		energy_inc_old = 0;
		energy_inc_new = 0;
		if (JcountSpin[i][1] != 1) {
			for (int j = 0; j < JcountSpin[i].size(); j++) {
				if (spins[JcountSpin[i][j]] == old_spin) { energy_inc_old += 1; }
				if (spins[JcountSpin[i][j]] == new_spin) { energy_inc_new += 1; }
			}
			energy_old += energy_inc_old * float(JMatSpin[i][0]);
			energy_new += energy_inc_new * float(JMatSpin[i][0]);
		}
		else {
			field_term = true;
			field_itr = i;
		}
	}
	if (field_term == false) { return (energy_new - energy_old); }
	else if (field_term == true) { return (energy_new - energy_old) + (JMatSpin[field_itr][0] * pow(new_spin, 2) - JMatSpin[field_itr][0] * pow(old_spin, 2)); }
}

float MCrun::calc_struct(int site, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_species = atom_species[site];
	int count = 0;
	if (site_species == 2) {
		for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
			if (neighbor_dist_list[site][i] == 0.5 or neighbor_dist_list[site][i] == 0.75) {
				if (atom_species[neighbor_index_list[site][i]] == 1) {
					count += 1;
				}
			}
			/*if (neighbor_dist_list[site][i] == 1.0 or neighbor_dist_list[site][i] == 1.25) {
				if
			}*/
		}
	}
	return count / 6.0;
}

float MCrun::calcWarrenCrowley(int site, float numb_Mn, float numb_In, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_species = atom_species[site];
	int count = 0;
	if (site_species == 1) {
		for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
			if (neighbor_dist_list[site][i] == 0.5 or neighbor_dist_list[site][i] == 0.75) {
				if (atom_species[neighbor_index_list[site][i]] == 2) {
					count += 1;
				}
			}
		}
	}
	return 1 - (count / 6.0) / (numb_In / (numb_In + numb_Mn));
}
// calculate the contribution to the spin-product order paramiter from a single site
float MCrun::calcMag2(int site, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_spin = atom_spin[site];
	float mag_product = 0;
	// loop through all neighbors in neighbor list
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (abs(neighbor_dist_list[site][i] - 0.707106769) < 0.001) {
			//if (neighbor_dist_list[site][i] == .5  or neighbor_dist_list[site][i] == 0.75) {
			mag_product += site_spin * atom_spin[neighbor_index_list[site][i]]; // if neighbor is directly above/below or directly adjacent... do the product
		}
	}
	return mag_product; // devide by 6 for normalization
}
// calculate the contribution to the spin-product order paramiter from a single site overload
float MCrun::calcMag2Diff(int site, int old_spin, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_spin = atom_spin[site];
	float mag_product = 0;
	float mag_product_old = 0;
	// loop through all neighbors in neighbor list
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (abs(neighbor_dist_list[site][i] - 0.707106769) < 0.001) {
			//if (neighbor_dist_list[site][i] == .5 or neighbor_dist_list[site][i] == 0.75) {
			mag_product += 2 * (site_spin - old_spin) * atom_spin[neighbor_index_list[site][i]]; // if neighbor is directly above/below or directly adjacent... do the product
		}
	}
	return (mag_product); // devide by 6 for normalization
}
// claculate the stag mag order paramiter contribution from a single atom
float MCrun::stagMag(int site, int spin, SimCell& sim_cell) {
	int mag = 0;
	int check = fmod(sim_cell.atom_list[site].pos[0] * sim_cell.atom_list[site].pos[1] * 2, 2); // check to see what sublattice the site is on
	if (check == 0) { mag += spin; }
	else { mag += 0.0; }
	return mag;
}

float MCrun::stagMagDiff(int site, int spin, int old_spin, SimCell& sim_cell) {
	int mag = 0;
	int check = fmod(sim_cell.atom_list[site].pos[0] * sim_cell.atom_list[site].pos[1] * 2, 2); // check to see what sublattice the site is on
	if (check == 0) { mag += (spin - old_spin); }
	else { mag += 0.0; }
	return mag;
}
// Curently unused!!
float MCrun::calcSpecies(int site, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_species = 0;
	int neighbor_species = 0;
	if (atom_species[site] == 1) { site_species = 1; }
	else if (atom_species[site] == 2) { site_species = -1; }
	float species_product = 0;
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (neighbor_dist_list[site][i] == .5 or neighbor_dist_list[site][i] == 0.75) {
			if (atom_species[neighbor_index_list[site][i]] == 1) { neighbor_species = 1; }
			else if (atom_species[neighbor_index_list[site][i]] == 2) { neighbor_species = -1; }
			species_product += site_species * neighbor_species;
		}
	}
	return species_product / 6;
}
// Curently unused!!
long double MCrun::delSiteEnergySpin(int site, int old_spin, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	//	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
	// for 2 atom terms
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ",";
		rule_itr = rule_map_spin.find(key);
		if (site_spin == 0) {
			enrg -= (rule_itr != rule_map_spin.end()) ? rule_itr->second * old_spin * neighbor_spin1 : 0.0;
		}
		else {
			enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * 2 * site_spin * neighbor_spin1 : 0.0;
		}
	}
	return enrg;
}
// Evaluate the total energy of the simulation cell in it's current state
long double MCrun::evalLattice(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	long double enrg = 0;
	// loop thruogh simulation cell and add everything up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;
}
// Evaluate the spin contribution to the total energy of the simulation cell
long double MCrun::evalLatticeSpin(map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	long double enrg = 0;
	// loop through everything and add it up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergySpin(site, rule_map_spin, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;// / atom_species.size() * 16 - 66.8295069760671;
}
// Run MC with ONLY spin flips
void MCrun::runMetropolis3() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	cout << atom_species_list.size() << ", ";
	cout << atom_spin_list.size() << "\n";
	// for (int i = 0; i < atom_species_list.size(); i++) { cout << atom_species_list[i] << ", "; }
	// cout << neighbor_index_list[0].size() << "\n";
	// for (int i = 0; i < atom_spin_list.size(); i++) { cout << atom_species_list[i] << ", "; }
	// cout << neighbor_dist_list[0].size() << "\n";

	// Start MC stuff
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	float spin_avg_Mn = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	vector<vector<float>> spin_states = session.spin_states;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// make sure STATES are given in input file. If not, abort.
	if (spin_states.empty()) {
		cout << "STATES not set\n";
		return;
	}
	if (spin_states.size() == 1 and spin_states[0][0] == 0) {
		cout << "something is wrong with the SPIN_STATES fiel\n";
		return;
	}
	// create seperate output file to avoid race condition
	cout << "Creating outputs\n";
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin SRO MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: ";
	for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
	Output << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_5\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	Output << "SRO Target: " << sro_target << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	cout << "finding initial SRO\n";
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		cout << "Setting SRO\n";
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";


	// setup part two: populate CountMat, JMat for faster MC run
	vector<vector<vector<int>>> CountMatChem;
	vector<vector<vector<int>>> CountMatSpin;
	vector<vector<int>> JcountsSpin;
	vector<vector<int>> JcountsChem;
	vector<int> JcSpin;
	vector<int> JcChem;
	vector<vector<float>> JMatChem;
	vector<vector<float>> JMatSpin;
	vector<float> JMChem;
	vector<float> JMSpin;
	vector<int> species_check;
	vector<float> dist_check;
	vector<int> spin_atoms;
	cout << "Making JMatSpin and Chem\n";
	for (int i = 0; i < session.rule_list.size(); i++) {
		if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() < 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() >= 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 1) {
			if (session.rule_list[i].GetLength() == 1) { spin_atoms.push_back(session.rule_list[i].GetDeco()[0]); }
			else {
				spin_atoms.push_back(session.rule_list[i].GetDeco()[0]);
				spin_atoms.push_back(session.rule_list[i].GetDeco()[1]);
			}
			JMSpin.push_back(session.rule_list[i].GetEnrgCont());
			JMSpin.push_back(session.rule_list[i].GetLength());
			JMatSpin.push_back(JMSpin);
			JMSpin.clear();
		}
	}
	cout << "Made JMatSpin and Chem\n";
	// make a list of atoms with a spin
	sort(spin_atoms.begin(), spin_atoms.end());
	spin_atoms.erase(unique(spin_atoms.begin(), spin_atoms.end()), spin_atoms.end());

	// cout << numb_atoms << "\n";
	for (int i = 0; i < numb_atoms; i++) {
		for (int j = 0; j < session.rule_list.size(); j++) {
			if (session.rule_list[j].GetLength() == 1) {
				if (session.rule_list[j].GetDeco()[0] == atom_species_list[i]) {
					if (session.rule_list[j].GetType() == 0) { JcChem.push_back(i); }
					else if (session.rule_list[j].GetType() == 1) { JcSpin.push_back(i); }
				}
			}
			else if (session.rule_list[j].GetLength() == 2) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					species_check.push_back(atom_species_list[i]);
					species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
					dist_check.push_back(neighbor_dist_list[i][k]);
					if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); } //////////////////// Could this if / else be a problem (is the type check inside the function messing up?) 
					else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); } ///////////////
					//else { cout << session.rule_list[j].GetType() << ", " << species_check[0] << species_check[1] << ", " << neighbor_dist_list[i][k]; }
					species_check.clear();
					dist_check.clear();
				}
			}
			if (session.rule_list[j].GetLength() == 3) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					for (int a = 0; a < neighbor_index_list[i].size(); a++) {
						species_check.push_back(atom_species_list[i]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][a]]);
						dist_check.push_back(neighbor_dist_list[i][k]);
						dist_check.push_back(neighbor_dist_list[i][a]);
						vector<int>::iterator dist_3_itr = find(neighbor_index_list[neighbor_index_list[i][k]].begin(), neighbor_index_list[neighbor_index_list[i][k]].end(), neighbor_index_list[i][a]);
						if (dist_3_itr != neighbor_index_list[neighbor_index_list[i][k]].end()) {
							int dist_3_ind = distance(neighbor_index_list[neighbor_index_list[i][k]].begin(), dist_3_itr);
							dist_check.push_back(neighbor_dist_list[neighbor_index_list[i][k]][dist_3_ind]);
							if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); }
							else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); }
						}
						species_check.clear();
						dist_check.clear();
					}
				}
			}
			if (session.rule_list[j].GetType() == 0) { JcountsChem.push_back(JcChem); }
			else if (session.rule_list[j].GetType() == 1) { JcountsSpin.push_back(JcSpin); }
			if (JcChem.size() != 0) { JcChem.clear(); }
			if (JcSpin.size() != 0) { JcSpin.clear(); }
		}
		CountMatChem.push_back(JcountsChem);
		CountMatSpin.push_back(JcountsSpin);
		if (JcountsChem.size() != 0) { JcountsChem.clear(); }
		if (JcountsSpin.size() != 0) { JcountsSpin.clear(); }
	}
	for (int i = 0; i < numb_atoms; i++) { for (int j_count = 0; j_count < CountMatChem[0].size(); j_count++) { cout << CountMatChem[i][j_count].size() << ", "; } }
	cout << "\n";
	for (int i = 0; i < numb_atoms; i++) { for (int j_count = 0; j_count < CountMatSpin[0].size(); j_count++) { cout << CountMatSpin[i][j_count].size() << ", "; } }
	cout << "\n";
	cout << "\n ok made it to here: " << CountMatChem.size() << ", " << CountMatSpin.size() << "\n";
	vector<vector<vector<int>>> PCountMatChem;
	vector<vector<vector<int>>> PCountMatSpin;
	for (int i = 0; i < CountMatChem.size(); i++) {
		for (int j = 0; j < CountMatChem[i].size(); j++) {
			JcChem = CountMatChem[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatChem[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatChem[neighbor_index_list[i][k]][j][a]) { JcChem.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsChem.push_back(JcChem);
			if (JcChem.size() > 0) { JcChem.clear(); }
		}
		PCountMatChem.push_back(JcountsChem);
		if (JcountsChem.size() > 0) { JcountsChem.clear(); }
	}

	for (int i = 0; i < CountMatSpin.size(); i++) {
		for (int j = 0; j < CountMatSpin[i].size(); j++) {
			JcSpin = CountMatSpin[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatSpin[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatSpin[neighbor_index_list[i][k]][j][a]) { JcSpin.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsSpin.push_back(JcSpin);
			if (JcSpin.size() > 0) { JcSpin.clear(); }
		}
		PCountMatSpin.push_back(JcountsSpin);
		if (JcountsSpin.size() > 0) { JcountsSpin.clear(); }
	}
	cout << "P Counts:\n";
	for (int i = 0; i < numb_atoms; i++) { for (int j_count = 0; j_count < PCountMatChem[0].size(); j_count++) { cout << PCountMatChem[i][j_count].size() << ", "; } }
	cout << "\n";
	for (int i = 0; i < numb_atoms; i++) { for (int j_count = 0; j_count < PCountMatSpin[0].size(); j_count++) { cout << PCountMatSpin[i][j_count].size() << ", "; } }
	cout << "\n";



	// Begin REAL MC
	//long double initial_enrg = evalLatticeFast(atom_spin_list, CountMatChem, CountMatSpin, JMatChem, JMatSpin);
	long double initial_enrg = evalLattiveChemFast(CountMatChem, JMatChem) + evalPottsLatticeSpinFast(atom_spin_list, CountMatSpin, JMatSpin);
	long double initial_spin_cont = evalPottsLatticeSpinFast(atom_spin_list, CountMatSpin, JMatSpin);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, magMn, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	long double var_mag = 0.0;
	long double var_e = 0.0;
	long double avg_mag = 0.0;
	long double avg_mag1 = 0.0;
	long double avg_mag2 = 0.0;
	long double avg_mag3 = 0.0;
	long double avg_mag4 = 0.0;
	long double e_flip = 0;
	float mag = 0;
	float mag1 = 0;
	float mag2 = 0;
	float mag3 = 0;
	RunningStat rs_C;
	RunningStat rs_X;
	vector<vector<int>> q_counts(spin_states.size());
	vector<float> ::iterator state_itr_old;
	vector<float> ::iterator state_itr_new;
	cout << "setting up q_counts\n";
	for (int i = 0; i < spin_states.size(); i++) {
		for (int j = 0; j < spin_states[i].size(); j++) {
			q_counts[i].push_back(0.0);
		}
	}
	for (int site = 0; site < numb_atoms; site++) {
		for (int i = 0; i < spin_states[atom_species_list[site]].size(); i++) {
			if (spin_states[atom_species_list[site]][i] == atom_spin_list[site] and find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
				q_counts[atom_species_list[site]][i] += 1;
			}
		}
	}
	for (int i = 0; i < q_counts.size(); i++) {
		cout << "(";
		for (int j = 0; j < q_counts[i].size(); j++) {
			cout << q_counts[i][j] << " ";
		}
		cout << ") \n";
	}
	cout << "done with q_counts\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		avg_mag = 0.0;
		avg_mag1 = 0.0;
		avg_mag2 = 0.0;
		avg_mag3 = 0.0;
		avg_mag4 = 0.0;

		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				if (find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
					// Flip Spin
					old_spin = atom_spin_list[site];
					spin_same = true;
					while (spin_same == true) {
						spin_rand = unif(rng);
						for (int it_spin_state = 0; it_spin_state < spin_states[atom_species_list[site]].size(); it_spin_state++) {
							if (spin_rand > float(it_spin_state) * 1.0 / float(spin_states[atom_species_list[site]].size())) { new_spin = spin_states[atom_species_list[site]][it_spin_state]; }
						}
						//if (spin_rand > 0.5) { new_spin = -1; }
						//else { new_spin = 1; }
						if (new_spin != old_spin or spin_states[atom_species_list[site]].size() <= 1) { spin_same = false; }
					}
					atom_spin_list[site] = new_spin;
					e_flip = evalPottsDelEnergySpinFast(old_spin, new_spin, atom_spin_list, PCountMatSpin[site], JMatSpin);
					if (e_flip < 0) {
						flip_count += 1;
						state_itr_old = find(spin_states[atom_species_list[site]].begin(), spin_states[atom_species_list[site]].end(), old_spin);
						state_itr_new = find(spin_states[atom_species_list[site]].begin(), spin_states[atom_species_list[site]].end(), new_spin);
						q_counts[atom_species_list[site]][distance(spin_states[atom_species_list[site]].begin(), state_itr_old)] -= 1;
						q_counts[atom_species_list[site]][distance(spin_states[atom_species_list[site]].begin(), state_itr_new)] += 1;
					}
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * temp) * (e_flip));
						if (keep_rand < keep_prob) {
							flip_count2 += 1;
							state_itr_old = find(spin_states[atom_species_list[site]].begin(), spin_states[atom_species_list[site]].end(), old_spin);
							state_itr_new = find(spin_states[atom_species_list[site]].begin(), spin_states[atom_species_list[site]].end(), new_spin);
							q_counts[atom_species_list[site]][distance(spin_states[atom_species_list[site]].begin(), state_itr_old)] -= 1;
							q_counts[atom_species_list[site]][distance(spin_states[atom_species_list[site]].begin(), state_itr_new)] += 1;
						}
						else { atom_spin_list[site] = old_spin; e_flip = 0; }
					}
					initial_enrg += e_flip;
					if (pass >= passes * .2) {
						e_avg += initial_enrg;// / numb_atoms;  //(pow(numb_atoms, 2) * 0.8 * passes);
						rs_C.Push(initial_enrg);// / numb_atoms);
						mag = 0;
						mag1 = 0;
						mag2 = 0;
						mag3 = 0;
						for (int i = 0; i < spin_states.size(); i++) {
							if (spin_states[i].size() == 1 or find(spin_atoms.begin(), spin_atoms.end(), i) == spin_atoms.end()) { continue; }
							int q_max = *max_element(q_counts[i].begin(), q_counts[i].end());
							mag += (float(spin_states[i].size()) * float(q_max) - float(session.atom_numbs[i])) / (float(spin_states[i].size()) - 1.0);
							if (i == 0) { mag1 += (float(spin_states[i].size()) * float(q_max) - float(session.atom_numbs[i])) / (float(spin_states[i].size()) - 1.0); }
							else if (i == 1) { mag2 += (float(spin_states[i].size()) * float(q_max) - float(session.atom_numbs[i])) / (float(spin_states[i].size()) - 1.0); }
							else if (i == 2) { mag3 += (float(spin_states[i].size()) * float(q_max) - float(session.atom_numbs[i])) / (float(spin_states[i].size()) - 1.0); }
						}
						mag /= numb_atoms;
						mag1 /= numb_atoms;
						mag2 /= numb_atoms;
						mag3 /= numb_atoms;
						avg_mag += mag;
						avg_mag1 += mag1;
						avg_mag2 += mag2;
						avg_mag3 += mag3;
						avg_mag4 += (mag2 + mag3);
						rs_X.Push(mag);// / numb_atoms);
					}
				}
			}
		}
		e_avg /= double(numb_atoms * numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_mag = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_mag / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< avg_mag / (numb_atoms * passes * .8) << ", "
			<< avg_mag1 / (numb_atoms * passes * .8) << ", "
			<< avg_mag2 / (numb_atoms * passes * .8) << ", "
			<< avg_mag3 / (numb_atoms * passes * .8) << ", "
			<< avg_mag4 / (numb_atoms * passes * .8) << ", "
			//<< "|" <<q_counts[1][0] << ", " << q_counts[1][1] << ", " << q_counts[1][2] << ", " << q_counts[1][3] << ", " << q_counts[1][4] << "|,"
			<< var_e << ", "
			<< var_mag << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	cout << e_avg << ", " << evalPottsLatticeSpinFast(atom_spin_list, CountMatSpin, JMatSpin) << "\n";
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}

void MCrun::runMetropolis4() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double spin_avg_Mn = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	vector<vector<float>> spin_states = session.spin_states;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin SRO MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: ";
	for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
	Output << "\n";	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_5\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	Output << "SRO Target: " << sro_target << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		cout << "Setting SRO\n";
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";


	// setup part two: populate CountMat, JMat for faster MC run
	vector<vector<vector<int>>> CountMatChem;
	vector<vector<vector<int>>> CountMatSpin;
	vector<vector<int>> JcountsSpin;
	vector<vector<int>> JcountsChem;
	vector<int> JcSpin;
	vector<int> JcChem;
	vector<vector<float>> JMatChem;
	vector<vector<float>> JMatSpin;
	vector<float> JMChem;
	vector<float> JMSpin;
	vector<int> species_check;
	vector<float> dist_check;
	vector<int> spin_atoms;
	cout << "Making JMatSpin and Chem\n";
	for (int i = 0; i < session.rule_list.size(); i++) {
		if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() < 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() >= 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 1) {
			if (session.rule_list[i].GetLength() == 1) { spin_atoms.push_back(session.rule_list[i].GetDeco()[0]); }
			else {
				spin_atoms.push_back(session.rule_list[i].GetDeco()[0]);
				spin_atoms.push_back(session.rule_list[i].GetDeco()[1]);
			}
			JMSpin.push_back(session.rule_list[i].GetEnrgCont());
			JMSpin.push_back(session.rule_list[i].GetLength());
			JMatSpin.push_back(JMSpin);
			JMSpin.clear();
		}
	}
	cout << "Made JMatSpin and Chem\n";
	sort(spin_atoms.begin(), spin_atoms.end());
	spin_atoms.erase(unique(spin_atoms.begin(), spin_atoms.end()), spin_atoms.end());
	cout << numb_atoms << "\n";
	for (int i = 0; i < numb_atoms; i++) {
		for (int j = 0; j < session.rule_list.size(); j++) {
			if (session.rule_list[j].GetLength() == 1) {
				if (session.rule_list[j].GetDeco()[0] == atom_species_list[i]) {
					if (session.rule_list[j].GetType() == 0) { JcChem.push_back(i); }
					else if (session.rule_list[j].GetType() == 1) { JcSpin.push_back(i); }
				}
			}
			else if (session.rule_list[j].GetLength() == 2) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					species_check.push_back(atom_species_list[i]);
					species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
					dist_check.push_back(neighbor_dist_list[i][k]);
					if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); }
					else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); }
					species_check.clear();
					dist_check.clear();
				}
			}
			if (session.rule_list[j].GetLength() == 3) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					for (int a = 0; a < neighbor_index_list[i].size(); a++) {
						species_check.push_back(atom_species_list[i]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][a]]);
						dist_check.push_back(neighbor_dist_list[i][k]);
						dist_check.push_back(neighbor_dist_list[i][a]);
						vector<int>::iterator dist_3_itr = find(neighbor_index_list[neighbor_index_list[i][k]].begin(), neighbor_index_list[neighbor_index_list[i][k]].end(), neighbor_index_list[i][a]);
						if (dist_3_itr != neighbor_index_list[neighbor_index_list[i][k]].end()) {
							int dist_3_ind = distance(neighbor_index_list[neighbor_index_list[i][k]].begin(), dist_3_itr);
							dist_check.push_back(neighbor_dist_list[neighbor_index_list[i][k]][dist_3_ind]);
							if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); }
							else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); }
						}
						species_check.clear();
						dist_check.clear();
					}
				}
			}
			if (session.rule_list[j].GetType() == 0) { JcountsChem.push_back(JcChem); }
			else if (session.rule_list[j].GetType() == 1) { JcountsSpin.push_back(JcSpin); }
			if (JcChem.size() != 0) { JcChem.clear(); }
			if (JcSpin.size() != 0) { JcSpin.clear(); }
		}
		CountMatChem.push_back(JcountsChem);
		CountMatSpin.push_back(JcountsSpin);
		if (JcountsChem.size() != 0) { JcountsChem.clear(); }
		if (JcountsSpin.size() != 0) { JcountsSpin.clear(); }
	}


	vector<vector<vector<int>>> PCountMatChem;
	vector<vector<vector<int>>> PCountMatSpin;
	for (int i = 0; i < CountMatChem.size(); i++) {
		for (int j = 0; j < CountMatChem[i].size(); j++) {
			JcChem = CountMatChem[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatChem[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatChem[neighbor_index_list[i][k]][j][a]) { JcChem.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsChem.push_back(JcChem);
			if (JcChem.size() > 0) { JcChem.clear(); }
		}
		PCountMatChem.push_back(JcountsChem);
		if (JcountsChem.size() > 0) { JcountsChem.clear(); }
	}

	for (int i = 0; i < CountMatSpin.size(); i++) {
		for (int j = 0; j < CountMatSpin[i].size(); j++) {
			JcSpin = CountMatSpin[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatSpin[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatSpin[neighbor_index_list[i][k]][j][a]) { JcSpin.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsSpin.push_back(JcSpin);
			if (JcSpin.size() > 0) { JcSpin.clear(); }
		}
		PCountMatSpin.push_back(JcountsSpin);
		if (JcountsSpin.size() > 0) { JcountsSpin.clear(); }
	}



	cout << "Starting Real MC\n";
	// Begin REAL MC
	long double initial_enrg = evalLatticeFast(atom_spin_list, atom_species_list, CountMatChem, CountMatSpin, JMatChem, JMatSpin);
	long double initial_spin_cont = evalLatticeSpinFast(atom_spin_list, atom_species_list, CountMatSpin, JMatSpin);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, magMn, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	long double initial_spin = 0.0;
	long double initial_spin_Mn = 0.0;
	long double var_spin = 0.0;
	long double var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	for (int site = 0; site < numb_atoms; site++) {
		if (find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
			initial_spin += atom_spin_list[site];
			if (atom_species_list[site] == 1) {
				initial_spin_Mn += atom_spin_list[site];
			}
		}
	}

	cout << "spin " << initial_spin / numb_atoms << " spinMn " << initial_spin_Mn / sim_cell.species_numbs[1] << "\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				site_species = atom_species_list[site];
				if (find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
					// Flip Spin
					old_spin = atom_spin_list[site];
					spin_same = true;
					while (spin_same == true) {
						spin_rand = unif(rng);
						for (int it_spin_state = 0; it_spin_state < spin_states[atom_species_list[site]].size(); it_spin_state++) {
							if (spin_rand > float(it_spin_state) * 1.0 / float(spin_states[atom_species_list[site]].size())) { new_spin = spin_states[atom_species_list[site]][it_spin_state]; }
						}
						if (new_spin != old_spin) { spin_same = false; }
					}
					atom_spin_list[site] = new_spin;
					e_flip = evalDelEnergySpinFast(old_spin, new_spin, atom_spin_list, PCountMatSpin[site], JMatSpin);
					spin_flip = new_spin - old_spin;
					if (e_flip < 0) { flip_count += 1; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * temp) * (e_flip));
						if (keep_rand < keep_prob) { flip_count2 += 1; }
						else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; }
					}
					initial_enrg += e_flip;
					initial_spin += spin_flip;
					if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
					if (pass >= passes * .2) {
						e_avg += initial_enrg; // / (pow(numb_atoms, 2) * 0.8 * passes);
						rs_C.Push(initial_enrg);// / numb_atoms);
						spin_avg += initial_spin / (pow(numb_atoms, 2) * 0.8 * passes);
						if (atom_species_list[site] == 1) { spin_avg_Mn += initial_spin_Mn / (pow(numb_atoms, 2) * 0.8 * passes); }//(numb_atoms * sim_cell.species_numbs[1] * passes * 0.8); }
						rs_X.Push(initial_spin);// / numb_atoms);
					}
				}
			}
		}
		e_avg /= double(numb_atoms * numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< initial_enrg << ", "
			<< spin_avg << ", "// / (pow(numb_atoms, 2) * passes * .8) << ", "
			<< spin_avg_Mn << ", "// / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}

void MCrun::runMetropolis5() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, long double> rule_map_spin;
	map <string, long double> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < session.rule_list.size(); i++) {
		string rule_key = "_";
		if (session.rule_list[i].GetLength() == 1) {
			rule_key = "_" + to_string(session.rule_list[i].GetDeco()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, session.rule_list[i].GetEnrgCont()));
		}
		else if (session.rule_list[i].GetLength() == 2) {
			vector<int> species = session.rule_list[i].GetDeco();
			float dist = session.rule_list[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			if (species[0] != species[1]) {
				rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
				if (session.rule_list[i].GetType() == 0) {
					rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
				}
				else if (session.rule_list[i].GetType() == 1) {
					rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
				}
			}
		}
		else if (session.rule_list[i].GetLength() == 3) {
			vector<int> trip = session.rule_list[i].GetDeco();
			vector<float> dists = session.rule_list[i].GetDists();
			if (session.rule_list[i].GetType() == 0) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
				//cout << rule_key << "\n";
			}
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors


	int spec_count1 = 0;
	for (int spec = 0; spec < sim_cell.numb_atoms; spec++) { spec_count1 += sim_cell.atom_list[spec].getSpecies(); }
	cout << spec_count1 << "\n";

	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	int spec_count = 0;
	for (int spec = 0; spec < atom_species_list.size(); spec++) { spec_count += atom_species_list[spec]; }
	cout << spec_count << "\n";



	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	long double e_site = 0;
	long double e_site_new = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double squ_e_avg = 0.0;
	long double squ_spin_avg = 0.0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);


	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin actual MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_5\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	long double initial_enrg = evalLattice3(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin3(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << ", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";

	//Begin REAL MC
	initial_enrg = evalLattice3(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin3(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, mag2, mag3, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	float initial_spin = 0.0;
	float initial_spin_Mn = 0.0;
	float initial_spin_2_Mn = 0.0;
	float initial_stag_mag = 0;
	long double initail_squ_spin = 0.0;
	long double initial_squ_enrg = 0.0;
	long double var_spin = 0.0;
	long double var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	for (int site = 0; site < numb_atoms; site++) {
		initial_spin += atom_spin_list[site];
		initial_spin_2_Mn += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
		initial_stag_mag += stagMag(site, atom_spin_list[site], sim_cell);
		if (atom_species_list[site] == 1) {
			initial_spin_Mn += atom_spin_list[site];
		}
	}
	initail_squ_spin = pow(initial_spin, 2);
	initial_squ_enrg = pow(initial_enrg, 2);
	cout << "spin " << initial_spin / numb_atoms << " spin2 " << initial_spin_Mn / numb_atoms << " spin 2 prod " << initial_spin_2_Mn / numb_atoms << "\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		//squ_e_avg = 0.0;
		spin_avg = 0.0;
		//squ_spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		spin_avg_Mn2 = 0.0;
		stag_avg = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin3(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin3(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				spin_flip = new_spin - old_spin;
				spin_flipMn2 = calcMag2Diff(site, old_spin, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				stag_flip = stagMagDiff(site, new_spin, old_spin, sim_cell);
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; spin_flipMn2 = 0; stag_flip = 0; }
				}
				initial_enrg += e_flip;
				initial_spin += spin_flip;
				initial_spin_2_Mn += spin_flipMn2;
				initial_stag_mag += stag_flip;
				if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
				if (pass >= passes * .2) {
					e_avg += initial_enrg;// / numb_atoms;  //(pow(numb_atoms, 2) * 0.8 * passes);
					rs_C.Push(initial_enrg);// / numb_atoms);
					spin_avg += initial_spin / (numb_atoms * 0.8 * passes);// (pow(numb_atoms, 2) * 0.8 * passes);
					spin_avg_Mn2 += initial_spin_2_Mn;
					stag_avg += initial_stag_mag;
					if (atom_species_list[site] == 1) {
						spin_avg_Mn += initial_spin_Mn;
						rs_X.Push(initial_spin);// / numb_atoms);
						//Output << temp << "," << initial_enrg << "," << initial_spin << "\n";

					}
				}
			}
		}
		e_avg /= double(numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< abs(spin_avg) << ", "
			<< spin_avg_Mn2 / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< stag_avg / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}
// Uses seperate thermostat for Spin and Chem
void MCrun::runMetropolis6() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, long double> rule_map_spin;
	map <string, long double> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < session.rule_list.size(); i++) {
		string rule_key = "_";
		if (session.rule_list[i].GetLength() == 1) {
			rule_key = "_" + to_string(session.rule_list[i].GetDeco()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, session.rule_list[i].GetEnrgCont()));
		}
		else if (session.rule_list[i].GetLength() == 2) {
			vector<int> species = session.rule_list[i].GetDeco();
			float dist = session.rule_list[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
		else if (session.rule_list[i].GetLength() == 3) {
			vector<int> trip = session.rule_list[i].GetDeco();
			vector<float> dists = session.rule_list[i].GetDists();
			if (session.rule_list[i].GetType() == 0) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
	}

	// make atom_list more acessable (and a map) for species and spin and neighbors
	int spec_count1 = 0;
	for (int spec = 0; spec < sim_cell.numb_atoms; spec++) { spec_count1 += sim_cell.atom_list[spec].getSpecies(); }
	cout << spec_count1 << "\n";

	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	int spec_count = 0;
	for (int spec = 0; spec < atom_species_list.size(); spec++) { spec_count += atom_species_list[spec]; }
	cout << spec_count << "\n";



	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	long double e_site = 0;
	long double e_site_new = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double squ_e_avg = 0.0;
	long double squ_spin_avg = 0.0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// Initalice and calc SRO_init //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_7\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	long double initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << ", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";

	//Begin REAL MC
	initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, mag2, mag3, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	float initial_spin = 0.0;
	float initial_spin_2 = 0.0;
	float initial_stag_mag = 0;
	long double var_spin = 0.0;
	long double var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	for (int site = 0; site < numb_atoms; site++) {
		if (atom_species_list[site] == 1) {
			initial_spin += atom_spin_list[site];
			initial_spin_2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
			initial_stag_mag += stagMag(site, atom_spin_list[site], sim_cell);
		}
	}
	cout << " Spin Mn " << initial_spin / numb_atoms << " spin 2 prod " << initial_spin_2 / numb_atoms << "\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		spin_avg_Mn2 = 0.0;
		stag_avg = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				if (atom_species_list[site] == 1) {
					// Flip Spin
					old_spin = atom_spin_list[site];
					spin_same = true;
					e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					while (spin_same == true) {
						spin_rand = unif(rng);
						if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
						else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
						else { new_spin = 1; }
						if (new_spin != old_spin) { spin_same = false; }
					}
					atom_spin_list[site] = new_spin;
					e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					e_flip = e_site_new - e_site;
					spin_flip = new_spin - old_spin;
					spin_flipMn2 = calcMag2Diff(site, old_spin, atom_spin_list, neighbor_index_list, neighbor_dist_list);
					stag_flip = stagMagDiff(site, new_spin, old_spin, sim_cell);
					if (e_flip < 0) { flip_count += 1; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * temp) * (e_flip));
						if (keep_rand < keep_prob) { flip_count2 += 1; }
						else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; spin_flipMn2 = 0; stag_flip = 0; }
					}
					initial_enrg += e_flip;
					initial_spin += spin_flip;
					initial_spin_2 += spin_flipMn2;
					initial_stag_mag += stag_flip;
					if (pass >= passes * .2) {
						e_avg += initial_enrg;// / numb_atoms;  //(pow(numb_atoms, 2) * 0.8 * passes);
						rs_C.Push(initial_enrg);// / numb_atoms);
						spin_avg += initial_spin / (numb_atoms * 0.8 * passes);// (pow(numb_atoms, 2) * 0.8 * passes);
						spin_avg_Mn2 += initial_spin_2;
						stag_avg += initial_stag_mag;
						if (atom_species_list[site] == 1) {
							spin_avg_Mn += initial_spin;
							rs_X.Push(initial_spin);// / numb_atoms);
							//Output << temp << "," << initial_enrg << "," << initial_spin << "\n";

						}
					}
				}
			}
		}
		e_avg /= double(numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< abs(spin_avg) << ", "
			<< spin_avg_Mn2 / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< stag_avg / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}

void MCrun::runMetropolis7() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, long double> rule_map_spin;
	map <string, long double> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < session.rule_list.size(); i++) {
		string rule_key = "_";
		if (session.rule_list[i].GetLength() == 1) {
			rule_key = "_" + to_string(session.rule_list[i].GetDeco()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, session.rule_list[i].GetEnrgCont()));
		}
		else if (session.rule_list[i].GetLength() == 2) {
			vector<int> species = session.rule_list[i].GetDeco();
			float dist = session.rule_list[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
		else if (session.rule_list[i].GetLength() == 3) {
			vector<int> trip = session.rule_list[i].GetDeco();
			vector<float> dists = session.rule_list[i].GetDists();
			if (session.rule_list[i].GetType() == 0) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors


	int spec_count1 = 0;
	for (int spec = 0; spec < sim_cell.numb_atoms; spec++) { spec_count1 += sim_cell.atom_list[spec].getSpecies(); }
	cout << spec_count1 << "\n";

	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	int spec_count = 0;
	for (int spec = 0; spec < atom_species_list.size(); spec++) { spec_count += atom_species_list[spec]; }
	cout << spec_count << "\n";



	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	long double e_site = 0;
	long double e_site_new = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double squ_e_avg = 0.0;
	long double squ_spin_avg = 0.0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);


	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin actual MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_7\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	long double initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << ", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";

	//Begin REAL MC
	initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, mag2, mag3, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	float initial_spin = 0.0;
	float initial_spin_Mn = 0.0;
	float initial_spin_2_Mn = 0.0;
	float initial_stag_mag = 0;
	long double initail_squ_spin = 0.0;
	long double initial_squ_enrg = 0.0;
	long double var_spin = 0.0;
	long double var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	for (int site = 0; site < numb_atoms; site++) {
		initial_spin += atom_spin_list[site];
		initial_spin_2_Mn += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
		initial_stag_mag += stagMag(site, atom_spin_list[site], sim_cell);
		if (atom_species_list[site] == 1) {
			initial_spin_Mn += atom_spin_list[site];
		}
	}
	initail_squ_spin = pow(initial_spin, 2);
	initial_squ_enrg = pow(initial_enrg, 2);
	cout << "spin " << initial_spin / numb_atoms << " spin2 " << initial_spin_Mn / numb_atoms << " spin 2 prod " << initial_spin_2_Mn / numb_atoms << "\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		//squ_e_avg = 0.0;
		spin_avg = 0.0;
		//squ_spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		spin_avg_Mn2 = 0.0;
		stag_avg = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				spin_flip = new_spin - old_spin;
				spin_flipMn2 = calcMag2Diff(site, old_spin, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				stag_flip = stagMagDiff(site, new_spin, old_spin, sim_cell);
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; spin_flipMn2 = 0; stag_flip = 0; }
				}
				initial_enrg += e_flip;
				initial_spin += spin_flip;
				initial_spin_2_Mn += spin_flipMn2;
				initial_stag_mag += stag_flip;
				if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
				if (pass >= passes * .2) {
					e_avg += initial_enrg;// / numb_atoms;  //(pow(numb_atoms, 2) * 0.8 * passes);
					rs_C.Push(initial_enrg);// / numb_atoms);
					spin_avg += initial_spin / (numb_atoms * 0.8 * passes);// (pow(numb_atoms, 2) * 0.8 * passes);
					spin_avg_Mn2 += initial_spin_2_Mn;
					stag_avg += initial_stag_mag;
					if (atom_species_list[site] == 1) {
						spin_avg_Mn += initial_spin_Mn;
						rs_X.Push(initial_spin);// / numb_atoms);
						//Output << temp << "," << initial_enrg << "," << initial_spin << "\n";

					}
				}
			}
		}
		e_avg /= double(numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< abs(spin_avg) << ", "
			<< spin_avg_Mn2 / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< stag_avg / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}

void MCrun::runMetropolis8() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, long double> rule_map_spin;
	map <string, long double> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < session.rule_list.size(); i++) {
		string rule_key = "_";
		if (session.rule_list[i].GetLength() == 1) {
			rule_key = "_" + to_string(session.rule_list[i].GetDeco()[0]) + ",0,";
			rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
		}
		else if (session.rule_list[i].GetLength() == 2) {
			vector<int> species = session.rule_list[i].GetDeco();
			float dist = session.rule_list[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
		else if (session.rule_list[i].GetLength() == 3) {
			vector<int> trip = session.rule_list[i].GetDeco();
			vector<float> dists = session.rule_list[i].GetDists();
			if (session.rule_list[i].GetType() == 0) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	ofstream Output;
	Output.open("OUTPUT");
	// Obtain correct SRO //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_7\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << ", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	for (int i = 0; i < eq_passes; i++) {
		for (int site = 0; site < numb_atoms; site++) {
			// Flip Species
			if (atom_species_list[site] != 0) {
				int rand_index = site;
				while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
				int old_species_site = atom_species_list[site];
				int old_species_rand = atom_species_list[rand_index];
				if (old_species_site != old_species_rand) {
					sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
					atom_species_list[site] = old_species_rand;
					atom_species_list[rand_index] = old_species_site;
					sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
					sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
					if (sro_flip < 0) { sro_initial += sro_site_flip; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
						if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
						else {
							atom_species_list[site] = old_species_site;
							atom_species_list[rand_index] = old_species_rand;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << sro_final << "\n";
	Output << "SRO: " << sro_final << "\n";






	//Begin REAL MC
	initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	float initial_spin = 0.0;
	float initial_spin_Mn = 0.0;
	float initial_spin_2_Mn = 0.0;
	float initial_stag_mag = 0;
	for (int site = 0; site < numb_atoms; site++) {
		initial_spin += atom_spin_list[site];
		initial_spin_2_Mn += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
		initial_stag_mag += stagMag(site, atom_spin_list[site], sim_cell);
		if (atom_species_list[site] == 1) {
			initial_spin_Mn += atom_spin_list[site];
		}
	}
	cout << "spin " << initial_spin / numb_atoms << " spin2 " << initial_spin_Mn / numb_atoms << " spin 2 prod " << initial_spin_2_Mn / numb_atoms << "\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		stag_avg = 0;
		flip_count = 0;
		flip_count2 = 0;
		// randomize spins
		for (int spin_site = 0; spin_site <= atom_spin_list.size(); spin_site++) {
			spin_rand = unif(rng);
			if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
			else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
			else { new_spin = 1; }
			atom_spin_list[spin_site] = new_spin;
		}
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				spin_flip = new_spin - old_spin;
				spin_flipMn2 = calcMag2Diff(site, old_spin, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				stag_flip = stagMagDiff(site, new_spin, old_spin, sim_cell);
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; spin_flipMn2 = 0; stag_flip = 0; }
				}
				initial_enrg += e_flip;
				initial_spin += spin_flip;
				initial_spin_2_Mn += spin_flipMn2;
				initial_stag_mag += stag_flip;
				if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
				if (pass >= passes * .2) {
					e_avg += initial_enrg;
					spin_avg += initial_spin;
					spin_avg_Mn2 += initial_spin_2_Mn;
					stag_avg += initial_stag_mag;
					if (atom_species_list[site] == 1) {
						spin_avg_Mn += initial_spin_Mn;
					}
				}
			}
		}
		cout << temp << ": " << e_avg / (pow(numb_atoms, 2) * passes * .8) - 0 << ": " << numb_atoms << ": " << initial_spin_2_Mn / numb_atoms << ": " << e_avg << "\n";
		Output << temp;
		Output << " , ";
		Output << e_avg / (pow(numb_atoms, 2) * passes * .8) - 0;
		Output << " , ";
		Output << abs(spin_avg) / (pow(numb_atoms, 2) * passes * .8);
		Output << " , ";
		Output << abs(spin_avg_Mn) / (pow(sim_cell.species_numbs[1], 2) * passes * .8);
		Output << " , ";
		Output << spin_avg_Mn2 / (numb_atoms * sim_cell.species_numbs[1] * passes * .8);
		Output << " , ";
		Output << stag_avg / (numb_atoms * sim_cell.species_numbs[1] * passes * .8);
		Output << " , ";
		Output << flip_count;
		Output << " , ";
		Output << flip_count2;
		Output << "\n";
	}
	writeSuperCell("CONTCAR", atom_species_list, atom_spin_list, sim_cell);

	Output.close();


}

void MCrun::runMetropolis9() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, long double> rule_map_spin;
	map <string, long double> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < session.rule_list.size(); i++) {
		string rule_key = "_";
		if (session.rule_list[i].GetLength() == 1) {
			rule_key = "_" + to_string(session.rule_list[i].GetDeco()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, session.rule_list[i].GetEnrgCont()));
		}
		else if (session.rule_list[i].GetLength() == 2) {
			vector<int> species = session.rule_list[i].GetDeco();
			float dist = session.rule_list[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (session.rule_list[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
		else if (session.rule_list[i].GetLength() == 3) {
			vector<int> trip = session.rule_list[i].GetDeco();
			vector<float> dists = session.rule_list[i].GetDists();
			if (session.rule_list[i].GetType() == 0) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));

				rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
				rule_map_chem.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
			else if (session.rule_list[i].GetType() == 1) {
				rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
				rule_map_spin.insert(pair<string, long double>(rule_key, session.rule_list[i].GetEnrgCont()));
			}
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	long double e_site = 0;
	long double e_site_new = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double squ_e_avg = 0.0;
	long double squ_spin_avg = 0.0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);


	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin actual MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_9\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	long double initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << ", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));

	for (int h = 0; h < passes; h++) {
		for (int site = 0; site < numb_atoms; site++) {
			if (atom_species_list[site] != 0) {
				int rand_index = site;
				while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
				int old_species_site = atom_species_list[site];
				int old_species_rand = atom_species_list[rand_index];
				if (old_species_site != old_species_rand) {
					atom_species_list[site] = old_species_rand;
					atom_species_list[rand_index] = old_species_site;
					keep_rand = unif(rng);
					if (keep_rand < 0.5) { sro_initial += sro_site_flip; }
					else {
						atom_species_list[site] = old_species_site;
						atom_species_list[rand_index] = old_species_rand;
					}
				}
			}
		}

		sro_initial = 0;
		sro_final = 0;
		for (int i = 0; i < numb_atoms; i++) {
			sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
		}
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < numb_atoms; i++) {
			sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
		}
		initial_enrg = evalLattice_tri(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
		Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2])))
			<< ", ENRG: " << initial_enrg / numb_atoms + 0 << "\n";
	}
	Output.close();
	writeSuperCell("CONTCAR", atom_species_list, atom_spin_list, sim_cell);
}

void MCrun::runMetropolis10() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	long double e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	long double spin_avg_Mn = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	vector<vector<float>> spin_states = session.spin_states;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// create seperate output file to avoid race condition
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin SRO MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: ";
	for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
	Output << "\n";	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_5\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	Output << "SRO Target: " << sro_target << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	if (eq_passes > 0) {
		cout << "Setting SRO\n";
		for (int i = 0; i < eq_passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";


	// setup part two: populate CountMat, JMat for faster MC run
	vector<vector<vector<int>>> CountMatChem;
	vector<vector<vector<int>>> CountMatSpin;
	vector<vector<int>> JcountsSpin;
	vector<vector<int>> JcountsChem;
	vector<int> JcSpin;
	vector<int> JcChem;
	vector<vector<float>> JMatChem;
	vector<vector<float>> JMatSpin;
	vector<float> JMChem;
	vector<float> JMSpin;
	vector<int> species_check;
	vector<float> dist_check;
	vector<int> spin_atoms;
	cout << "Making JMatSpin and Chem\n";
	for (int i = 0; i < session.rule_list.size(); i++) {
		if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() < 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 0 and session.rule_list[i].GetLength() >= 3) {
			JMChem.push_back(session.rule_list[i].GetEnrgCont());
			JMChem.push_back(session.rule_list[i].GetLength());
			JMatChem.push_back(JMChem);
			JMChem.clear();
		}
		else if (session.rule_list[i].GetType() == 1) {
			if (session.rule_list[i].GetLength() == 1) { spin_atoms.push_back(session.rule_list[i].GetDeco()[0]); }
			else {
				spin_atoms.push_back(session.rule_list[i].GetDeco()[0]);
				spin_atoms.push_back(session.rule_list[i].GetDeco()[1]);
			}
			JMSpin.push_back(session.rule_list[i].GetEnrgCont());
			JMSpin.push_back(session.rule_list[i].GetLength());
			JMatSpin.push_back(JMSpin);
			JMSpin.clear();
		}
	}
	cout << "Made JMatSpin and Chem\n";
	sort(spin_atoms.begin(), spin_atoms.end());
	spin_atoms.erase(unique(spin_atoms.begin(), spin_atoms.end()), spin_atoms.end());
	cout << numb_atoms << "\n";
	for (int i = 0; i < numb_atoms; i++) {
		for (int j = 0; j < session.rule_list.size(); j++) {
			if (session.rule_list[j].GetLength() == 1) {
				if (session.rule_list[j].GetDeco()[0] == atom_species_list[i]) {
					if (session.rule_list[j].GetType() == 0) { JcChem.push_back(i); }
					else if (session.rule_list[j].GetType() == 1) { JcSpin.push_back(i); }
				}
			}
			else if (session.rule_list[j].GetLength() == 2) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					species_check.push_back(atom_species_list[i]);
					species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
					dist_check.push_back(neighbor_dist_list[i][k]);
					if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); }
					else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); }
					species_check.clear();
					dist_check.clear();
				}
			}
			if (session.rule_list[j].GetLength() == 3) {
				for (int k = 0; k < neighbor_index_list[i].size(); k++) {
					for (int a = 0; a < neighbor_index_list[i].size(); a++) {
						species_check.push_back(atom_species_list[i]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][k]]);
						species_check.push_back(atom_species_list[neighbor_index_list[i][a]]);
						dist_check.push_back(neighbor_dist_list[i][k]);
						dist_check.push_back(neighbor_dist_list[i][a]);
						vector<int>::iterator dist_3_itr = find(neighbor_index_list[neighbor_index_list[i][k]].begin(), neighbor_index_list[neighbor_index_list[i][k]].end(), neighbor_index_list[i][a]);
						if (dist_3_itr != neighbor_index_list[neighbor_index_list[i][k]].end()) {
							int dist_3_ind = distance(neighbor_index_list[neighbor_index_list[i][k]].begin(), dist_3_itr);
							dist_check.push_back(neighbor_dist_list[neighbor_index_list[i][k]][dist_3_ind]);
							if (session.rule_list[j].IsRuleChem(species_check, dist_check)) { JcChem.push_back(neighbor_index_list[i][k]); }
							else if (session.rule_list[j].IsRuleSpin(species_check, dist_check)) { JcSpin.push_back(neighbor_index_list[i][k]); }
						}
						species_check.clear();
						dist_check.clear();
					}
				}
			}
			if (session.rule_list[j].GetType() == 0) { JcountsChem.push_back(JcChem); }
			else if (session.rule_list[j].GetType() == 1) { JcountsSpin.push_back(JcSpin); }
			if (JcChem.size() != 0) { JcChem.clear(); }
			if (JcSpin.size() != 0) { JcSpin.clear(); }
		}
		CountMatChem.push_back(JcountsChem);
		CountMatSpin.push_back(JcountsSpin);
		if (JcountsChem.size() != 0) { JcountsChem.clear(); }
		if (JcountsSpin.size() != 0) { JcountsSpin.clear(); }
	}

	cout << "debugging... " << JMatChem.size() << " " << JMatChem[8].size() << " " << JMatChem[8][0] << " " << JMatChem[8][1] << '\n';
	vector<vector<vector<int>>> PCountMatChem;
	vector<vector<vector<int>>> PCountMatSpin;
	for (int i = 0; i < CountMatChem.size(); i++) {
		for (int j = 0; j < CountMatChem[i].size(); j++) {
			JcChem = CountMatChem[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatChem[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatChem[neighbor_index_list[i][k]][j][a]) { JcChem.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsChem.push_back(JcChem);
			if (JcChem.size() > 0) { JcChem.clear(); }
		}
		PCountMatChem.push_back(JcountsChem);
		if (JcountsChem.size() > 0) { JcountsChem.clear(); }
	}

	for (int i = 0; i < CountMatSpin.size(); i++) {
		for (int j = 0; j < CountMatSpin[i].size(); j++) {
			JcSpin = CountMatSpin[i][j];
			for (int k = 0; k < neighbor_index_list[i].size(); k++) {
				for (int a = 0; a < CountMatSpin[neighbor_index_list[i][k]][j].size(); a++) {
					if (i == CountMatSpin[neighbor_index_list[i][k]][j][a]) { JcSpin.push_back(neighbor_index_list[i][k]); }
				}
			}
			JcountsSpin.push_back(JcSpin);
			if (JcSpin.size() > 0) { JcSpin.clear(); }
		}
		PCountMatSpin.push_back(JcountsSpin);
		if (JcountsSpin.size() > 0) { JcountsSpin.clear(); }
	}



	cout << "Starting Real MC\n";
	// Begin REAL MC
	long double initial_enrg = evalLatticeFast(atom_spin_list, atom_species_list, CountMatChem, CountMatSpin, JMatChem, JMatSpin);
	long double initial_spin_cont = evalLatticeSpinFast(atom_spin_list, atom_species_list, CountMatSpin, JMatSpin);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, magMn, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	long double initial_spin = 0.0;
	long double initial_spin_Mn = 0.0;
	long double var_spin = 0.0;
	long double var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	for (int site = 0; site < numb_atoms; site++) {
		if (find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
			initial_spin += atom_spin_list[site];
			if (atom_species_list[site] == 1) {
				initial_spin_Mn += atom_spin_list[site];
			}
		}
	}

	cout << "spin " << initial_spin / numb_atoms << " spinMn " << initial_spin_Mn / sim_cell.species_numbs[1] << "\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		spin_avg_Mn = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				site_species = atom_species_list[site];
				if (find(spin_atoms.begin(), spin_atoms.end(), atom_species_list[site]) != spin_atoms.end()) {
					// Flip Spin
					old_spin = atom_spin_list[site];
					spin_same = true;
					while (spin_same == true) {
						spin_rand = unif(rng);
						for (int it_spin_state = 0; it_spin_state < spin_states[atom_species_list[site]].size(); it_spin_state++) {
							if (spin_rand > float(it_spin_state) * 1.0 / float(spin_states[atom_species_list[site]].size())) { new_spin = spin_states[atom_species_list[site]][it_spin_state]; }
						}
						if (new_spin != old_spin) { spin_same = false; }
					}
					atom_spin_list[site] = new_spin;
					e_flip = evalDelEnergySpinFast(old_spin, new_spin, atom_spin_list, PCountMatSpin[site], JMatSpin);
					spin_flip = new_spin - old_spin;
					if (e_flip < 0) { flip_count += 1; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * temp) * (e_flip));
						if (keep_rand < keep_prob) { flip_count2 += 1; }
						else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; }
					}
					initial_enrg += e_flip;
					initial_spin += spin_flip;
					if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
					if (pass >= passes * .2) {
						e_avg += initial_enrg; // / (pow(numb_atoms, 2) * 0.8 * passes);
						rs_C.Push(initial_enrg);// / numb_atoms);
						spin_avg += initial_spin / (pow(numb_atoms, 2) * 0.8 * passes);
						if (atom_species_list[site] == 1) { spin_avg_Mn += initial_spin_Mn / (pow(numb_atoms, 2) * 0.8 * passes); }//(numb_atoms * sim_cell.species_numbs[1] * passes * 0.8); }
						rs_X.Push(initial_spin);// / numb_atoms);
					}
				}
			}
		}
		e_avg /= double(numb_atoms * numb_atoms * 0.8 * passes);
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< initial_enrg << ", "
			<< spin_avg << ", "// / (pow(numb_atoms, 2) * passes * .8) << ", "
			<< spin_avg_Mn << ", "// / (numb_atoms * sim_cell.species_numbs[1] * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}


long double MCrun::evalSiteEnergyAll_tri(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	// Initalize variables
	string key;
	long double enrg = 0;
	map<string, long double>::iterator rule_itr;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int neighbor_site1;
	int neighbor_species1;
	//int neighbor_spin1;
	float neighbor_dist1;
	int neighbor_site2;
	int neighbor_species2;
	//int neighbor_spin2;
	int dist_2_ind;
	float neighbor_dist2;
	float neighbor_dist3;
	enrg -= session.mag_ext * session.moments[site_species] * uB * site_spin;
	key = "_" + to_string(site_species) + ",0,"; // make key for one atom term
	rule_itr = rule_map_chem.find(key);
	enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule for one atom term
	// for 2 atom terms
	// loop through all neighbors in the neighbor list (list only includes atoms that have distances listed in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		if (neighbor_site1 == site) { continue; }
		neighbor_species1 = atom_species[neighbor_site1];
		//neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // make key string
		rule_itr = rule_map_chem.find(key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // if the key was in the chem rule map, apply the rule, otherwise do nothing
		//rule_itr = rule_map_spin.find(key);
		//enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * float(site_spin * neighbor_spin1) : 0.0; // if the key was in the spin rule map, apply the rule, otherwise do nothing
		//// for 3 atom terms
		//// loop through all neighbors in the neighbor list
		for (int j = 0; j < neighbor_index_list[site].size(); j++) {
			neighbor_site2 = neighbor_index_list[site][j];
			if (neighbor_site1 == neighbor_site2 or neighbor_site2 == site) { continue; }
			neighbor_species2 = atom_species[neighbor_site2];
			//neighbor_spin2 = atom_spin[neighbor_site2];
			vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
			if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
				dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
				neighbor_dist2 = neighbor_dist_list[site][j];
				neighbor_dist3 = neighbor_dist_list[neighbor_site1][dist_2_ind];
				// make key for three atom term
				key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
				key += "," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
				rule_itr = rule_map_chem.find(key);
				enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / 6 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
				//rule_itr = rule_map_spin.find(key);
				//enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 * 1 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
			}
		}
	}
	return enrg;
}
// Evaluate the total energy of the simulation cell in it's current state
long double MCrun::evalLattice_tri(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	long double enrg = 0;
	// loop thruogh simulation cell and add everything up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergyAll_tri(site, rule_map_spin, rule_map_chem, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;
}
// Write the current state of the simulation cell to a txt file. Mostly for debugging.
void MCrun::writeSuperCell(string out_path, vector<int>& atom_species, vector<int>& atom_spins, SimCell& sim_cell) {
	ofstream OUT_file;
	OUT_file.open("CONTCAR");
	if (OUT_file.is_open()) {
		OUT_file << "Ni2Mn(1-X)In(X)\n" << 5.76 * sim_cell.sup_cell[0] << "\n";
		OUT_file << sim_cell.unit_LC[0] << " 0 0\n";
		OUT_file << "0 " << sim_cell.unit_LC[1] << " 0\n" << "0 0 ";
		OUT_file << sim_cell.unit_LC[2] << "\n";
		OUT_file << "Ni Mn In\n";
		for (int i = 0; i < sim_cell.species_numbs.size(); i++) { OUT_file << sim_cell.species_numbs[i] << " "; }
		OUT_file << "\n";
		for (int i = 0; i < sim_cell.numb_atoms; i++) {
			OUT_file << sim_cell.atom_list[i].pos[0] / (sim_cell.sup_cell[0] * sim_cell.unit_LC[0]) << " ";
			OUT_file << sim_cell.atom_list[i].pos[1] / (sim_cell.sup_cell[1] * sim_cell.unit_LC[1]) << " ";
			OUT_file << sim_cell.atom_list[i].pos[2] / (sim_cell.sup_cell[2] * sim_cell.unit_LC[2]) << " ";
			OUT_file << " # " << atom_species[i] << " " << atom_spins[i] << "\n";
		}
	}
	OUT_file.close();
}

void MCrun::outputSimCell() {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	cout << atom_species_list.size() << ", ";
	cout << atom_spin_list.size() << "\n";
	// for (int i = 0; i < atom_species_list.size(); i++) { cout << atom_species_list[i] << ", "; }
	// cout << neighbor_index_list[0].size() << "\n";
	// for (int i = 0; i < atom_spin_list.size(); i++) { cout << atom_species_list[i] << ", "; }
	// cout << neighbor_dist_list[0].size() << "\n";

	// Start MC stuff
	long double new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	long double e_avg = 0.0;
	long double spin_avg = 0.0;
	float spin_avg_Mn = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	long double Cmag = 0.0;
	long double Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float thermostat = session.sro_temp;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	vector<vector<float>> spin_states = session.spin_states;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// make sure STATES are given in input file. If not, abort.
	if (spin_states.empty()) {
		cout << "STATES not set\n";
		return;
	}
	if (spin_states.size() == 1 and spin_states[0][0] == 0) {
		cout << "something is wrong with the SPIN_STATES fiel\n";
		return;
	}
	// create seperate output file to avoid race condition
	cout << "Creating outputs\n";
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);

	// begin SRO MC run //
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: ";
	for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
	Output << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_5\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	Output << "SRO Target: " << sro_target << "\n";
	//cout << "SRO Target: " << sro_target << "\n";
	//float sro_final = 0;
	//float sro_initial = 0;
	//float sro_site_new = 0;
	//float sro_site_old = 0;
	//float sro_site_flip = 0;
	//float sro_flip = 0;
	//cout << "finding initial SRO\n";
	//for (int i = 0; i < numb_atoms; i++) {
	//	sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	//}
	//if (eq_passes > 0) {
	//	cout << "Setting SRO\n";
	//	for (int i = 0; i < eq_passes; i++) {
	//		for (int site = 0; site < numb_atoms; site++) {
	//			// Flip Species
	//			if (atom_species_list[site] != 0) {
	//				int rand_index = site;
	//				while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
	//				int old_species_site = atom_species_list[site];
	//				int old_species_rand = atom_species_list[rand_index];
	//				if (old_species_site != old_species_rand) {
	//					sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
	//					sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
	//					atom_species_list[site] = old_species_rand;
	//					atom_species_list[rand_index] = old_species_site;
	//					sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
	//					sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
	//					sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
	//					sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
	//					if (sro_flip < 0) { sro_initial += sro_site_flip; }
	//					else {
	//						keep_rand = unif(rng);
	//						keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
	//						if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
	//						else {
	//							atom_species_list[site] = old_species_site;
	//							atom_species_list[rand_index] = old_species_rand;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//for (int i = 0; i < numb_atoms; i++) {
	//	sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	//}
	//cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	//Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";

	writeSuperCell("CONTCAR" + to_string(outFile_count), atom_species_list, atom_spin_list, sim_cell);

	Output.close();
}

void MCrun::debug_run() {
	// create seperate output file to avoid race condition
	cout << "Starting debug...\n";
	string file_name = "OUTPUT";
	bool file_exists = true;
	int outFile_count = 0;
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outFile_count += 1;
			file_name = "OUTPUT" + to_string(outFile_count);
		}
		else {
			file_exists = false;
			close(fd);
		}
	}
	const char* c_file = file_name.c_str();
	ofstream Output;
	Output.open(c_file);
	cout << "attempting to print output...\n";
	// begin SRO MC run //
	Output << "Rule Info: \n";
	Output << "Rule file:" << session.rules_file << "\n";
	Output << "Algo: " << session.algo <<"\n";
	Output << "Number: " << session.rule_list.size() << "\n";
	Output << "Energy cont Rule[0]: " << session.rule_list[0].GetEnrgCont() << "\n";
	Output << "Deco Rule[0]: ";
	for (int i = 0; i < session.rule_list[6].deco.size(); i++) { Output << session.rule_list[6].deco[i] << ","; }
	Output << "\n" << sim_cell.atom_list[0].allowed_species.size() << "atoms allowed\n";
	Output << "\n";
	Output << "Motif Rule[0]: ";
	Output << "\n";
	Output.close();
}