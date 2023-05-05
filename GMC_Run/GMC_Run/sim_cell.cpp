#include "sim_cell.h"

// This is the sim_cell class. It holds the state of the entire simulation cell. In the MC implementation it's main components are 
// coppied to a seperate vector to minimize the size or complexity of arguments being passed to certian functions. Still, it does 
// serve as a usefull too for initalizing the simulation cell and as a reffernce.
SimCell::SimCell(void) {
	sim_type = "EMPTY";
}

// the sim_cell object
//SimCell::SimCell(string POSCAR_file, int _sup_cell[3], vector<int>& _species_numbs, float _cutoff, string _sim_type, string _phase_init, string spin_init, string species_init) {
//	sim_type = _sim_type;
//	phase_init = _phase_init;
//	cutoff = _cutoff;
//	bool use_poscar = true;
//	for (int i = 0; i < _species_numbs.size(); i++) {
//		species_types.push_back(i);
//		species_numbs.push_back(_species_numbs[i]);
//	}
//	fillUnitCell(POSCAR_file,  use_poscar); // make a "unit cell" from poscar file. This is used in turn to make the simulation cell
//	for (int i = 0; i < 3; i++) { // set the dimensions of the simulation cell
//		sup_cell[i] = _sup_cell[i];
//		cell_dim[i] = unit_LC[i] * sup_cell[i];
//	}
//	vector<vector<float>> _pos_list;
//	vector<int> _species_list;
//	make_supercell(_pos_list, _species_list, phase_init); // create the simulation cell from the unit cell.
//	vector<float> dist_list{ 0 };
//	fillAtomList(_pos_list, _species_list, dist_list, phase_init, spin_init, species_init); // populate the simulation cell with "atom" objects initalized to the desired spin, species and phase settings
//}

SimCell::SimCell(SimCell& sc_copy) {
	cutoff = sc_copy.cutoff;
	sim_type = sc_copy.sim_type;
	phase_init = sc_copy.phase_init;
	species_types = sc_copy.species_types;
	species_numbs = sc_copy.species_numbs;
	poscar_comp = sc_copy.poscar_comp;
	copy(sc_copy.sup_cell, sc_copy.sup_cell + 3, sup_cell);
	copy(sc_copy.cell_dim, sc_copy.cell_dim + 3, cell_dim);
	numb_atoms = sc_copy.numb_atoms;
	copy(sc_copy.numb_cells, sc_copy.numb_cells + 3, numb_cells);
	copy(sc_copy.unit_LC, sc_copy.unit_LC + 3, unit_LC);
	X_num = sc_copy.X_num;
	atom_list = sc_copy.atom_list;
	unit_cell = sc_copy.unit_cell;
}

void SimCell::_copy(SimCell& sc_copy) {
	cutoff = sc_copy.cutoff;
	sim_type = sc_copy.sim_type;
	phase_init = sc_copy.phase_init;
	species_types = sc_copy.species_types;
	species_numbs = sc_copy.species_numbs;
	poscar_comp = sc_copy.poscar_comp;
	copy(sc_copy.sup_cell, sc_copy.sup_cell + 3, sup_cell);
	copy(sc_copy.cell_dim, sc_copy.cell_dim + 3, cell_dim);
	numb_atoms = sc_copy.numb_atoms;
	copy(sc_copy.numb_cells, sc_copy.numb_cells + 3, numb_cells);
	copy(sc_copy.unit_LC, sc_copy.unit_LC + 3, unit_LC);
	X_num = sc_copy.X_num;
	atom_list = sc_copy.atom_list;
	unit_cell = sc_copy.unit_cell;
}
// fill the unit cell using the poscar file 
void SimCell::fillUnitCell(string POSCAR_file,vector<string>& _species, bool use_poscar) {
	ifstream POS_list;
	string pos_line;
	vector<string> comp_line;
	vector<string> pos_lines;
	vector<string> LCs;
	vector<string> pos_list_s;
	vector<string> pos_species_s;
	vector<vector<float>> pos_list_f;
	vector<vector<string>> pos_species_f;
	vector<int> species_list;
	vector<float> pos{ 0,0,0 };
	// read the poscar file
	POS_list.open(POSCAR_file, ifstream::in);
	if (POS_list.is_open()) {
		while (getline(POS_list, pos_line))
		{
			pos_lines.push_back(pos_line);
		}
		POS_list.close();
	}
	else cout << "Unable to open POSCAR file\n";
	// get the lattice constants              NEEDS TO BE CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for (int i = 2; i < 5; i++) {
		pos_line = pos_lines[i];
		LCs = split(pos_line);// , " ");
		unit_LC[i - 2] = stof(LCs[i - 2]);
	}
	// get the positions of each atom
	for (int i = 8; i < pos_lines.size(); i++) {
		pos_line = pos_lines[i];
		pos_list_s = split(pos_line);// , " ");
		for (int j = 0; j < 3; j++) {
			pos[j] = stof(pos_list_s[j]);
		}
		for (int j = 3; j < pos_list_s.size(); j++) {
			pos_species_s.push_back(pos_list_s[j]);
		}
		if (pos_species_s.size() == 0) {
			pos_species_f.push_back(_species);
		}
		else {
			pos_species_f.push_back(pos_species_s);
		}
		pos_list_f.push_back(pos);
	}
	// assign the correct atomic species to each atom position
	X_num = stoi(split(pos_lines[6])[0]);
	comp_line = split(pos_lines[6]);
	for (int i = 0; i < comp_line.size(); i++) {
		poscar_comp.push_back(stoi(comp_line[i]));
	}
	cout << use_poscar << "\n";
	if (use_poscar == false) {
		for (int i = 0; i < pos_list_f.size(); i++) {
			if (i < X_num) {
				species_list.push_back(0);
				cout << "what??!? \n";
			}
			else {
				species_list.push_back(1);
				cout << "How?!?!?! \n";
			}
		}
	}
	else {
		int comp_count = 0;
		int comp_running_sum = poscar_comp[comp_count];
		for (int i = 0; i < pos_list_f.size(); i++) {
			if (i < comp_running_sum) {
				species_list.push_back(comp_count);
			}
			else {
				bool inc_cc = true;
				while (inc_cc == true)
				{
					comp_count += 1;
					if (poscar_comp[comp_count] > 0) { inc_cc = false; }
				}
				comp_running_sum += poscar_comp[comp_count];
				species_list.push_back(comp_count);
			}
		}
	}
	// create "atom" object for each atom in the unit cell
	for (int i = 0; i < pos_list_f.size(); i++) {
		unit_cell.push_back(Atom(i, species_list[i], 0, 0, pos_list_f[i]));
	}
}

// make a super cell (simulation cell) from the unit cell 
void SimCell::make_supercell(vector<vector<float>>& _pos_list, vector<int>& _species_list, string phase_init) {
	int x = sup_cell[0];
	int y = sup_cell[1];
	int z = sup_cell[2];
	int current_cell[3];
	vector<float> new_atom_pos{ 0,0,0 };
	const int unit_length = unit_cell.size();
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				current_cell[0] = i;
				current_cell[1] = j;
				current_cell[2] = k;
				for (int m = 0; m < unit_length; m++) {
					for (int n = 0; n < 3; n++) {
						new_atom_pos[n] = round((unit_cell[m].pos[n] + current_cell[n]) * unit_LC[n] * 100000) / 100000;
					}
					_pos_list.push_back(new_atom_pos);
					_species_list.push_back(unit_cell[m].getSpecies());
				}
			}
		}
	}
	//for (int i = 0; i < poscar_comp.size(); i++) {
	//	if (poscar_comp[i] * x * y * z != species_numbs[i]) {
	//		int del_num = poscar_comp[i] * x * y * z - species_numbs[i];
	//		std::mt19937_64 rng_int;
	//		uint64_t timeSeed_int = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	//		std::seed_seq ss_int{ uint32_t(timeSeed_int & 0xffffffff), uint32_t(timeSeed_int >> 32) };
	//		rng_int.seed(ss_int);
	//		std::uniform_int_distribution<int> unif_int(0, _species_list.size() - 1);
	//		int flip, check, index_rand;
	//		for (int i = 0; i < abs(del_num); i++) {
	//			if (del_num > 0) { flip = 1; check = 2; }
	//			else { flip = 2; check = 1; }
	//			bool check_site = false;
	//			while (check_site == false)
	//			{
	//				index_rand = unif_int(rng_int);
	//				if (_species_list[index_rand] == check) { check_site = true; }
	//			}
	//			_species_list[index_rand] = flip;
	//		}
	//	}
	//}
	if (poscar_comp[2] * x * y * z != species_numbs[2]) {
		cout << "Changing POSCAR comp to match input file\n";
		int del_Znum = poscar_comp[2] * x * y * z - species_numbs[2];
		std::mt19937_64 rng_int;
		uint64_t timeSeed_int = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::seed_seq ss_int{ uint32_t(timeSeed_int & 0xffffffff), uint32_t(timeSeed_int >> 32) };
		rng_int.seed(ss_int);
		std::uniform_int_distribution<int> unif_int(0, _species_list.size() - 1);
		int flip, check, index_rand;
		for (int i = 0; i < abs(del_Znum); i++) {
			if (del_Znum > 0) { flip = 1; check = 2; }
			else { flip = 2; check = 1; }
			bool check_site = false;
			while (check_site == false)
			{
				index_rand = unif_int(rng_int);
				if (_species_list[index_rand] == check) { check_site = true; }
			}
			_species_list[index_rand] = flip;
		}
	}
}

// fill the neighbor matrix using the dist_list (list of relevent neighbors from the mc_rules)
void SimCell::setNeighborDists(vector<float>& dist_list) {
	// declare variables
	float pos1[3];
	float pos2[3];
	float distXYZ[3];
	float dist;
	int added = 0;
	// loop through each atom object
	for (int i = 0; i < atom_list.size(); i++) {
		pos1[0] = atom_list[i].pos[0];
		pos1[1] = atom_list[i].pos[1];
		pos1[2] = atom_list[i].pos[2];
		// loop through all neighbors
		for (int j = 0; j < atom_list.size(); j++) {
			pos2[0] = atom_list[j].pos[0];
			pos2[1] = atom_list[j].pos[1];
			pos2[2] = atom_list[j].pos[2];
			// find the xyz distance between the atoms
			distXYZ[0] = pos2[0] - pos1[0];
			distXYZ[1] = pos2[1] - pos1[1];
			distXYZ[2] = pos2[2] - pos1[2];
			// apply BC
			for (int k = 0; k < 3; k++) {
				if (distXYZ[k] > .5 * cell_dim[k]) {
					distXYZ[k] = cell_dim[k] - distXYZ[k];
				}
				else if (distXYZ[k] < -.5 * cell_dim[k]) {
					distXYZ[k] = cell_dim[k] + distXYZ[k];
				}
			}
			// find straight line distance in BC
			dist = sqrt(pow(distXYZ[0], 2) + pow(distXYZ[1], 2) + pow(distXYZ[2], 2));
			added = 0;
			// add to neighbor distance list and neighbor index list if the distance is relevent
			for (int k = 0; k < dist_list.size(); k++) {
				if (fabs(dist - dist_list[k]) <= .001) {
					atom_list[i].neighbor_dists.push_back(dist);
					atom_list[i].neighbors.push_back(j);
					added = 1;
				}
			}
		}
	}
}

// populate atom list with correctly initalized atom objects
void SimCell::fillAtomList(vector<vector<float>>& _pos_list, vector<int>& _species_list, vector<float> dist_list, string phase_init, string spin_init, string species_init) {
	// declare or init variables
	int atom_index = 0;
	int spin;
	int phase = 0; // just a place holder
	double spin_rand;
	int index_rand;
	bool use_rand = false;
	vector<int> species_list;

	// Random number for if spin init is random
	// THIS USES A RANDOM NUMB GEN AND MAY CAUSE MAC PROBLEMS !!!!! //
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// loop through each atom in the pos list
	for (int i = 0; i < _pos_list.size(); i++) {
		spin_rand = unif(rng);
		// if spin init is FM
		if (spin_init[0] == 'F') {
			spin = 1;
		}
		// if spin init is random, randomly use -1,0,1 ---------------------------------------------------------------------------This need re-tooling to accomidate SPIN_STATES file
		else if (spin_init[0] == 'R') {
			//if (spin_rand >= (.6666666666666666)) { spin = -1; }
			//if (spin_rand >= (.3333333333333333) and spin_rand < (.6666666666666666)) { spin = 0; }
			//if (spin_rand < (.3333333333333333)) { spin = 1; }
			if (spin_rand >= (0.5)) { spin = -1; }
			else { spin = 1; }
		}
		// AFM init is more complicated. I don't realy need it and can do it later
		else if (spin_init[0] == 'A') {
			spin = 1; // just a place holder
		}
		atom_list.push_back(Atom(i, _species_list[i], spin, phase, _pos_list[i]));
	}
	numb_atoms = atom_list.size();
	// Random number generater for random species placement
	// THIS USES A RANDOM NUMB GEN AND MAY CAUSE MAC PROBLEMS !!!!! //
	//std::mt19937_64 rng_int;
	//uint64_t timeSeed_int = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	//std::seed_seq ss_int{ uint32_t(timeSeed_int & 0xffffffff), uint32_t(timeSeed_int >> 32) };
	//rng_int.seed(ss_int);
	//std::uniform_int_distribution<int> unif_int(0, numb_atoms - 1);
	//
	// if species init is rand... else use init from poscar
	//if (species_init[0] == 'R' || use_rand == true) {
	//	int numb_comp = 0;
	//	while (numb_comp < species_numbs[2]) {
	//		index_rand = unif_int(rng_int);
	//		if (atom_list[index_rand].getSpecies() != 0) {
	//			if (atom_list[index_rand].getSpecies() != 2) {
	//				atom_list[index_rand].setSpecies(2);
	//				numb_comp += 1;
	//			}
	//		}
	//	}
	//}
	setNeighborDists(dist_list); // add neighbor distance for each atom
}

float SimCell::findAtomDists(int atom1, int atom2) {
	auto it = std::find(atom_list[atom1].neighbors.begin(), atom_list[atom1].neighbors.end(), atom2);
	return distance(atom_list[atom1].neighbors.begin(), it);
}

// create the sim cell object and initalize it accordingly
void SimCell::initSimCell(string POSCAR_file, Session& session) {

	int _sup_cell[3] = { session.shape[0], session.shape[1], session.shape[2] };
	vector<int>& _species_numbs = session.atom_numbs;
	float _cutoff = 1; // Double check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	string _sim_type = session.sim_type;
	string _phase_init = session.phase_init;
	string spin_init = session.spin_init;
	string species_init = session.species_init;
	bool use_poscar = session.use_poscar;
	sim_type = _sim_type;
	phase_init = _phase_init;
	cutoff = _cutoff;
	for (int i = 0; i < _species_numbs.size(); i++) {
		species_types.push_back(i);
		species_numbs.push_back(_species_numbs[i]);
	}
	fillUnitCell(POSCAR_file,session.species, use_poscar);
	for (int i = 0; i < 3; i++) {
		sup_cell[i] = _sup_cell[i];
		cell_dim[i] = unit_LC[i] * sup_cell[i];
	}
	vector<vector<float>> _pos_list;
	vector<int> _species_list;
	make_supercell(_pos_list, _species_list, phase_init);
	fillAtomList(_pos_list, _species_list, session.unique_dists, phase_init, spin_init, species_init);
}

// atom object inside SimCell class. used as a container for information about single atoms. Standard constructor and get/set functions are used
SimCell::Atom::Atom(void) {
	species = 10;
	spin = 10;
	phase = 10;
	index = 10;
}

SimCell::Atom::Atom(int _index, int _species, int _spin, int _phase, vector<float> _pos) {
	index = _index;
	species = _species;
	spin = _spin;
	phase = _phase;
	pos[0] = _pos[0];
	pos[1] = _pos[1];
	pos[2] = _pos[2];
}

int SimCell::Atom::getNeighborSpin(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_spin = sim_cell.atom_list[neighbor_index].getSpin();
	return neighbor_spin;
}

int SimCell::Atom::getNeighborSpecies(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_species = sim_cell.atom_list[neighbor_index].getSpecies();
	return neighbor_species;
}

int SimCell::Atom::getNeighborPhase(int _neighbor, SimCell& sim_cell) {
	int neighbor_index = neighbors[_neighbor];
	int neighbor_phase = sim_cell.atom_list[neighbor_index].getPhase();
	return neighbor_phase;
}

int SimCell::Atom::getNeighborIndex(int _neighbor, SimCell& sim_cell) {
	return neighbors[_neighbor];
}

float SimCell::Atom::getNeighborDist(int _neighbor, SimCell& sim_cell) {
	return neighbor_dists[_neighbor];
}

int SimCell::Atom::getNumbNeighbors(int _site, SimCell& sim_cell) {
	return sim_cell.atom_list[_site].neighbors.size();
}

int SimCell::Atom::getSpin(void) {
	return spin;
}

int SimCell::Atom::getSpecies(void) {
	return species;
}

int SimCell::Atom::getPhase(void) {
	return phase;
}

void SimCell::Atom::setSpin(int _spin) {
	spin = _spin;
}

void SimCell::Atom::setSpecies(int _species) {
	species = _species;
}

void SimCell::Atom::setPhase(int _phase) {
	phase = _phase;
}
