#include "algo2.h"

Algo2::Algo2(void) {}

Algo2::Algo2(Session& _session, SimCell& _sim_cell) {
	session = _session;
	sim_cell = _sim_cell;
}

float Algo2::eval_site_chem(int site) {
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	for (int i = 0; i < chem_motif_groups[site].size(); i++) {
		string rule_key = "0.";
		rule_key += to_string(i);
		for (int j : chem_motif_groups[site][i]) {
			rule_key += "." + to_string(chem_list[j]);
		}
		rule_itr = rule_map_chem.find(rule_key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0;
	}
	return enrg;
}

float Algo2::eval_site_spin(int site) {
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	for (int i = 0; i < spin_motif_groups[site].size(); i++) {
		float spin_prod = 1;
		string rule_key = "1.";
		rule_key += to_string(i);
		for (int j : spin_motif_groups[site][i]) {
			rule_key += "." + to_string(chem_list[j]);
			spin_prod *= spin_list[j];
		}
		rule_itr = rule_map_spin.find(rule_key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod : 0.0;
	}
	return enrg;
}

float Algo2::eval_spin_flip(int site, float old_spin) {
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	for (int i = 0; i < spin_motif_groups[site].size(); i++) {
		float spin_prod = 1;
		string rule_key = "1.";
		rule_key += to_string(i);
		for (int j : spin_motif_groups[site][i]) {
			rule_key += "." + to_string(chem_list[j]);
			if (j != site) { spin_prod *= spin_list[j]; }
		}
		rule_itr = rule_map_spin.find(rule_key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod : 0.0;
	}
	return 2 * (enrg * spin_list[site] - enrg * old_spin);
}

float Algo2::eval_chem_flip(int site1, int site2) {
	int site1_old_chem = chem_list[site1];
	int site2_old_chem = chem_list[site2];
	int site1_old_spin = spin_list[site1];
	int site2_old_spin = spin_list[site2];
}

float Algo2::eval_lat() {
	float enrg = 0;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		enrg += eval_site_chem(site);
		enrg += eval_site_spin(site);
	}
	return enrg + session.intercept;
}

float Algo2::eval_lat_spin() {
	float enrg = 0;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		enrg += eval_site_spin(site);
	}
	return enrg;
}

void Algo2::run() {
	//declair variables
	float e_flip = 0;
	float spin_flip = 0;
	float spin_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	float e_avg = 0.0;
	float spin_avg = 0.0;
	int flip_count = 0;
	int flip_count2 = 0;
	float Cmag = 0.0;
	float Xmag = 0.0;
	int passes = session.numb_passes;
	int eq_passes = session.eq_passes;
	float sro_target = session.sro_target;
	float temp1 = session.start_temp;
	float temp2 = session.end_temp;
	float temp_inc = session.temp_inc;
	float keep_rand;
	float keep_prob;
	vector<vector<float>> spin_states = session.spin_states;
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neigh_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<int> spin_atoms;

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

	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: ";
	for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
	Output << "\n";	Output << "MC passes: " << session.numb_passes << "\n";
	Output << "Beginning MC EQ run using Algo2\n";

	// make atom_list more acessable for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		chem_list.push_back(sim_cell.atom_list[i].getSpecies());
		spin_list.push_back(sim_cell.atom_list[i].getSpin());
		pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
		for (int j = 0; j < numb_neighbors; j++) {
			neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neigh_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}

	// make rule_maps for easy lookup
	string rule_key;
	for (Rule rule : session.chem_rule_list) {
		rule_key = to_string(rule.GetType()) + "." + to_string(rule.motif_ind);
		for (int i = 0; i < rule.deco.size(); i++) { rule_key += "." + to_string(rule.deco[i]); }
		rule_map_chem.insert(pair<string, float>(rule_key, rule.GetEnrgCont()));
	}
	for (Rule rule : session.spin_rule_list) {
		rule_key = to_string(rule.GetType()) + "." + to_string(rule.motif_ind);
		for (int i = 0; i < rule.deco.size(); i++) { rule_key += "." + to_string(rule.deco[i]); }
		rule_map_spin.insert(pair<string, float>(rule_key, rule.GetEnrgCont()));
		for (int atom : rule.deco) {
			if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); }
		}
	}


	// fill motif group lists
	fill_CMG(neigh_ind_list);
	fill_SMG(neigh_ind_list);

	// initalize system with desired SRO
	Output << "EQ passes: " << session.eq_passes << ", EQ Temp: " << session.sro_temp << "\n";
	Output << "SRO Target: " << session.sro_target << "\n";
	cout << "SRO Target: " << session.sro_target << "\n";
	//float sro_final = init_SRO(neigh_ind_list, neigh_dist_list);
	//cout << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";
	//Output << "SRO: " << 1 - (sro_final) / (float(sim_cell.species_numbs[1]) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]))) << "\n";

	cout << "Starting Real MC\n";
	// Begin MC

	float init_enrg = eval_lat();
	cout << "evaluated lattice\n";
	float init_spin_cont = eval_lat_spin();
	cout << "evalueated spin contribution\n";
	Output << init_enrg / numb_atoms + 0 << ", " << init_spin_cont / numb_atoms << "\n";
	cout << init_enrg / numb_atoms + 0 << ", " << init_spin_cont / numb_atoms << "\n";
	Output << "temp, enrg, mag, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
	float init_spin = 0.0;
	float var_spin = 0.0;
	float var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	//RunningStat init_enrg_R;
	cout << "counting spins...\n";
	for (int site = 0; site < numb_atoms; site++) {
		if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
			init_spin += spin_list[site];
		}
	}
	cout << "spins counted\n";
	cout << "spin is " << init_spin / numb_atoms << "per atom\n";
	float inc_dir = 1;
	if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }

	// setup rng
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// start MC loop
	cout << "entering main loop\n";
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
					// Flip Spin
					old_spin = spin_list[site];
					spin_same = true;
					while (spin_same == true) {
						spin_rand = unif(rng);
						for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
							if (spin_rand > float(it_spin_state) * 1.0 / float(spin_states[chem_list[site]].size())) { new_spin = spin_states[chem_list[site]][it_spin_state]; }
						}
						if (new_spin != old_spin) { spin_same = false; }
					}
					spin_list[site] = new_spin;
					e_flip = eval_spin_flip(site, old_spin);
					spin_flip = new_spin - old_spin;
					if (e_flip < 0) { flip_count += 1; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * temp) * (e_flip));
						if (keep_rand < keep_prob) { flip_count2 += 1; }
						else { spin_list[site] = old_spin; e_flip = 0.0; spin_flip = 0; }

					}
					init_enrg += e_flip;
					init_spin += spin_flip;
					if (pass >= passes * .2) {
						e_avg += init_enrg; // / (pow(numb_atoms, 2) * 0.8 * passes);
						rs_C.Push(init_enrg);// / numb_atoms);
						spin_avg += init_spin / (pow(numb_atoms, 2) * 0.8 * passes);
						rs_X.Push(init_spin);// / numb_atoms);
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
			<< spin_avg << ", "// / (pow(numb_atoms, 2) * passes * .8) << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	cout << " MC Finished\n";
	print_state();
	Output.close();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float Algo2::calc_struct(int site, vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list) {
	int site_species = chem_list[site];
	int count = 0;
	if (site_species == 2) {
		for (int i = 0; i < neigh_dist_list[site].size(); i++) {
			if (neigh_dist_list[site][i] == 0.5 or neigh_dist_list[site][i] == 0.75) {
				if (chem_list[neigh_ind_list[site][i]] == 1) {
					count += 1;
				}
			}
		}
	}
	return count / 6.0;
}

float Algo2::init_SRO(vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list) {
	// setup rng
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// begin SRO MC run //
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	float keep_rand = 0;
	float keep_prob = 0;


	//sro_target = (1 - sro_target) * (float(sim_cell.species_numbs[2])) / (float(sim_cell.species_numbs[1] + sim_cell.species_numbs[2]));
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		sro_initial += calc_struct(i, neigh_ind_list, neigh_dist_list) / sim_cell.species_numbs[2];
	}
	if (session.eq_passes > 0) {
		cout << "Setting SRO\n";
		for (int i = 0; i < session.eq_passes; i++) {
			for (int site = 0; site < sim_cell.numb_atoms; site++) {
				// Flip Species
				if (chem_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or chem_list[rand_index] == 0) { rand_index = rand() % sim_cell.numb_atoms; }
					int old_species_site = chem_list[site];
					int old_species_rand = chem_list[rand_index];
					if (old_species_site != old_species_rand) {
						sro_site_old = calc_struct(site, neigh_ind_list, neigh_dist_list);
						sro_site_old += calc_struct(rand_index, neigh_ind_list, neigh_dist_list);
						chem_list[site] = old_species_rand;
						chem_list[rand_index] = old_species_site;
						sro_site_new = calc_struct(site, neigh_ind_list, neigh_dist_list);
						sro_site_new += calc_struct(rand_index, neigh_ind_list, neigh_dist_list);
						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
						sro_flip = abs((sro_site_flip + sro_initial) - session.sro_target) - abs(sro_initial - session.sro_target);
						if (sro_flip < 0) { sro_initial += sro_site_flip; }
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb * session.sro_temp) * (sro_flip));
							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
							else {
								chem_list[site] = old_species_site;
								chem_list[rand_index] = old_species_rand;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		sro_final += calc_struct(i, neigh_ind_list, neigh_dist_list) / sim_cell.species_numbs[2];
	}
	return sro_final;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Algo2::bc_check(vector<float> check_vect, vector<float>& pos) {
	bool bc_test = false;
	vector<int> dir { -1, 1 };
	vector<float> bc_pos { 0, 0, 0 };
	vector<float> lc_shift { 0, 0, 0 };
	vector<vector<float>> bc_shifts { pos };

	for (int a : dir) {
		for (int i = 0; i < 3; i++) {
			lc_shift = scale_vect(sim_cell.lat_vect[i], a);
			bc_pos = pos_shift(pos, lc_shift);
			bc_shifts.push_back(bc_pos);
		}
		for (int b : dir) {
			lc_shift = scale_vect(sim_cell.lat_vect[0], a);
			bc_pos = pos_shift(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[1], b);
			bc_pos = pos_shift(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			lc_shift = scale_vect(sim_cell.lat_vect[0], a);
			bc_pos = pos_shift(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[2], b);
			bc_pos = pos_shift(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			lc_shift = scale_vect(sim_cell.lat_vect[1], a);
			bc_pos = pos_shift(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[2], b);
			bc_pos = pos_shift(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			for (int c : dir) {
				lc_shift = scale_vect(sim_cell.lat_vect[0], a);
				bc_pos = pos_shift(pos, lc_shift);
				lc_shift = scale_vect(sim_cell.lat_vect[1], b);
				bc_pos = pos_shift(bc_pos, lc_shift);
				lc_shift = scale_vect(sim_cell.lat_vect[2], c);
				bc_pos = pos_shift(bc_pos, lc_shift);
				bc_shifts.push_back(bc_pos);
			}
		}
	}
	for (vector<float> check : bc_shifts) {
		if (pos_comp(check, check_vect)) {
			bc_test = true;
			break;
		}
	}
	return bc_test;
}

void Algo2::fill_SMG(vector<vector<int>>& neigh_ind_list) {
	vector<int> deco_group;
	vector<float> new_pos { 0.0, 0.0, 0.0 };
	vector<float> self_site { 0.0, 0.0, 0.0 };
	vector<vector<int>> motifs;
	int last_ind = -1;
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		for (Rule rule : session.spin_rule_list) {
			if (last_ind != rule.motif_ind) {
				for (vector<float> shift : rule.motif) {
					if (pos_comp(shift, self_site)) { deco_group.push_back(i); }
					else {
						new_pos = pos_shift(pos_list[i], shift);
						for (int neigh : neigh_ind_list[i]) {
							int x = 0;
							if (bc_check(pos_list[neigh], new_pos)) {
								deco_group.push_back(neigh);
							}
						}
					}
				}
				motifs.push_back(deco_group);
				deco_group.clear();
				last_ind = rule.motif_ind;
			}
		}
		spin_motif_groups.push_back(motifs);
		motifs.clear();
	}
}

void Algo2::fill_CMG(vector<vector<int>>& neigh_ind_list) {
	vector<int> deco_group;
	vector<float> new_pos { 0.0, 0.0, 0.0 };
	vector<float> self_site { 0.0, 0.0, 0.0 };
	vector<vector<int>> motifs;
	int last_ind = -1;
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		for (Rule rule : session.chem_rule_list) {
			if (last_ind != rule.motif_ind) {
				for (vector<float> shift : rule.motif) {
					if (pos_comp(shift, self_site)) { deco_group.push_back(i); }
					else {
						new_pos = pos_shift(pos_list[i], shift);
						for (int neigh : neigh_ind_list[i]) {
							int x = 0;
							if (bc_check(pos_list[neigh], new_pos)) {
								deco_group.push_back(neigh);
							}
						}
					}
				}
				motifs.push_back(deco_group);
				deco_group.clear();
				last_ind = rule.motif_ind;
			}
		}
		chem_motif_groups.push_back(motifs);
		motifs.clear();
	}
}

void Algo2::print_state() {
	vector<int> perm;
	vector<int> temp_spec;
	vector<float> temp_spin;
	vector<vector<float>> temp_pos;
	temp_spec.assign(chem_list.begin(), chem_list.end());
	temp_spin.assign(spin_list.begin(), spin_list.end());
	temp_pos = pos_list;// insert(temp_pos.end(), pos_list.begin(), pos_list.end());
	perm = vect_permut(temp_spin);
	sort_vect(temp_spin, perm);
	sort_vect(temp_pos, perm);
	sort_vect(temp_spec, perm);
	perm = vect_permut(temp_spec);
	sort_vect(temp_spin, perm);
	sort_vect(temp_pos, perm);
	sort_vect(temp_spec, perm);
	ofstream OUT_file;
	OUT_file.open("CONTCAR");
	if (OUT_file.is_open()) {
		OUT_file << "Alloy of";
		for (string spec : session.species_str) { OUT_file << " " << spec; }
		OUT_file << "\n 1 \n";
		for (int i = 0; i < 3; i++) {
			vector<float>vect = sim_cell.unit_lat_vect[i];
			for (int j = 0; j < 3; j++) { OUT_file << vect[j] * session.shape[i] << " "; }
			OUT_file << "\n";
		}
		for (string spec : session.species_str) { OUT_file << " " << spec; }
		OUT_file << "\n";
		for (int i = 0; i < sim_cell.species_numbs.size(); i++) { OUT_file << sim_cell.species_numbs[i] << " "; }
		OUT_file << "\nCartesian\n";
		for (int i = 0; i < sim_cell.numb_atoms; i++) {
			OUT_file << temp_pos[i][0] << " " << temp_pos[i][1] << " " << temp_pos[i][2] << " ";
			OUT_file << " # " << temp_spec[i] << " " << temp_spin[i] << "\n";
		}
	}
	OUT_file.close();
}