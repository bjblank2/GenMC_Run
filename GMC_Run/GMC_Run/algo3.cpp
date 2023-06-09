#include "algo3.h"

Algo3::Algo3(void) {}

Algo3::Algo3(Session& _session, SimCell& _sim_cell) {
	session = _session;
	sim_cell = _sim_cell;
	session.fill_sro_list();
	int max_sro_type = 0;
	bool zero_ind_check = false;
	for (int i = 0; i < session.sro_rule_list.size(); i++) {
		max_sro_type = (session.sro_rule_list[i].GetType() > max_sro_type) ? session.sro_rule_list[i].GetType() : max_sro_type;
		if (session.sro_rule_list[i].GetType() == 0) {
			zero_ind_check = true;
		}
	}
	for (int i = 0; i < max_sro_type + 1; i++) {
		sro_list.push_back(0);
		sro_flip.push_back(0);
	}

	if (session.numb_subpasses < 1) {
		cout << "_______________________________________________________________________________" << endl;
		cout << "Possible Error: Algo3 has been given 0 subpasses per main pass" << endl;
		cout << "This implies that no seperate spin flips are made which is likely not phisical." << endl;
		cout << "This is probably not what you want..." << endl;
		cout << "_______________________________________________________________________________" << endl;
	}

	if (zero_ind_check != true) {
		cout << "_______________________________________________________________________________" << endl;
		cout << "Possible Error: Algo3 has been given SRO types that are not zero-indexed" << endl;
		cout << "This will result in inccorect tracking of all SRO parameters" << endl;
		cout << "_______________________________________________________________________________" << endl;
	}
}

void Algo3::eval_lat_sro() {
	vector<int> site_sro;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		site_sro = eval_site_sro(site);
		for (int i = 0; i < site_sro.size(); i++) {
			sro_list[i] += site_sro[i];
		}
	}
}

float Algo3::eval_site_chem(int site) {
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

float Algo3::eval_site_spin(int site) {
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

float Algo3::eval_spin_flip(int site, float old_spin) {
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
		enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod) / spin_motif_groups[site].size() : 0.0;
	}
	return (enrg * spin_list[site] - enrg * old_spin);
}

float Algo3::eval_atom_flip(int site1, int site2) {
	int new_spec1 = chem_list[site1];
	int new_spec2 = chem_list[site2];
	float new_enrg = 0;
	float old_enrg = 0;
	float new_spin1 = spin_list[site1];
	float new_spin2 = spin_list[site2];
	vector<int> old_site_sro(sro_list.size(), 0);
	vector<int> new_site_sro(sro_list.size(), 0);
	map<string, float>::iterator rule_itr;
	// new and old enrgey for site1
	// for chemestry motifs
	for (int i = 0; i < chem_motif_groups[site1].size(); i++) {
		string new_chem_key = "0." + to_string(i);
		string old_chem_key = "0." + to_string(i);
		for (int j : chem_motif_groups[site1][i]) {
			if (j != site1) {
				new_chem_key += "." + to_string(chem_list[j]);
				old_chem_key += "." + to_string(chem_list[j]);
			}
			else {
				new_chem_key += "." + to_string(new_spec1);
				old_chem_key += "." + to_string(new_spec2);
			}
		}
		rule_itr = rule_map_chem.find(new_chem_key);
		new_enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / chem_motif_groups[site1].size() : 0.0;
		if (rule_itr != rule_map_sro.end()) { new_site_sro[rule_itr->second] += 1; }
		rule_itr = rule_map_chem.find(old_chem_key);
		old_enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / chem_motif_groups[site1].size() : 0.0;
		if (rule_itr != rule_map_sro.end()) { old_site_sro[rule_itr->second] += 1; }
	}
	// for spin motifs
	for (int i = 0; i < spin_motif_groups[site1].size(); i++) {
		float spin_prod = 1;
		string new_spin_key = "1." + to_string(i);
		string old_spin_key = "1." + to_string(i);
		for (int j : spin_motif_groups[site1][i]) {
			if (j != site1) {
				new_spin_key += "." + to_string(chem_list[j]);
				old_spin_key += "." + to_string(chem_list[j]);
				spin_prod *= spin_list[j];
			}
			else {
				new_spin_key += "." + to_string(new_spec1);
				old_spin_key += "." + to_string(new_spec2);
			}
		}
		rule_itr = rule_map_spin.find(new_spin_key);
		new_enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod * new_spin1 / spin_motif_groups[site1].size() : 0.0;
		rule_itr = rule_map_spin.find(old_spin_key);
		old_enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod * new_spin2 / spin_motif_groups[site1].size() : 0.0;
	}

	// new and old energy for site2
	// for chemestry motifs
	for (int i = 0; i < chem_motif_groups[site2].size(); i++) {
		string new_chem_key = "0." + to_string(i);
		string old_chem_key = "0." + to_string(i);
		for (int j : chem_motif_groups[site2][i]) {
			if (j != site2) {
				new_chem_key += "." + to_string(chem_list[j]);
				old_chem_key += "." + to_string(chem_list[j]);
			}
			else {
				new_chem_key += "." + to_string(new_spec2);
				old_chem_key += "." + to_string(new_spec1);        
			}
		}
		rule_itr = rule_map_chem.find(new_chem_key);
		new_enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / chem_motif_groups[site2].size() : 0.0;
		if (rule_itr != rule_map_sro.end()) { new_site_sro[rule_itr->second] += 1; }
		rule_itr = rule_map_chem.find(old_chem_key);
		old_enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / chem_motif_groups[site2].size() : 0.0;
		if (rule_itr != rule_map_sro.end()) { old_site_sro[rule_itr->second] += 1; }

	}
	// for spin motifs
	for (int i = 0; i < spin_motif_groups[site2].size(); i++) {
		float spin_prod = 1;
		string new_spin_key = "1." + to_string(i);
		string old_spin_key = "1." + to_string(i);
		for (int j : spin_motif_groups[site2][i]) {
			if (j != site2) {
				new_spin_key += "." + to_string(chem_list[j]);
				old_spin_key += "." + to_string(chem_list[j]);
				spin_prod *= spin_list[j];
			}
			else {
				new_spin_key += "." + to_string(new_spec2);
				old_spin_key += "." + to_string(new_spec1);
			}
		}
		rule_itr = rule_map_spin.find(new_spin_key);
		new_enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod * new_spin2 / spin_motif_groups[site2].size() : 0.0;
		rule_itr = rule_map_spin.find(old_spin_key);
		old_enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * spin_prod * new_spin1 / spin_motif_groups[site2].size() : 0.0;
	}
	old_site_sro = scale_vect(old_site_sro, -1);
	sro_flip = vect_add(new_site_sro,  old_site_sro);
	return (new_enrg - old_enrg);
}

float Algo3::eval_lat() {
	float enrg = 0;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		enrg += eval_site_chem(site);
		enrg += eval_site_spin(site);
	}
	return enrg + session.intercept;
}

float Algo3::eval_lat_spin() {
	float enrg = 0;
	for (int site = 0; site < sim_cell.numb_atoms; site++) {
		enrg += eval_site_spin(site);
	}
	return enrg;
}

void Algo3::run() {
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
	int sub_passes = session.numb_subpasses;
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
	while (file_exists == true) {
		const char* c_file = file_name.c_str();
		int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		if (fd < 0) {
			// file exists or otherwise uncreatable
			outfile_count += 1;
			file_name = "OUTPUT" + to_string(outfile_count);
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
	Output << "\n";	Output << "MC atom-flip passes: " << session.numb_passes << ", ";
	Output << "\n";	Output << "MC spin-flip passes: " << session.numb_subpasses << "\n";
	Output << "Beginning MC EQ run using Algo3\n";

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
	for (SRO rule : session.sro_rule_list) {
		rule_key = to_string(0) + "." + to_string(rule.motif_ind);
		for (int i = 0; i < rule.deco.size(); i++) { rule_key += "." + to_string(rule.deco[i]); }
		rule_map_sro.insert(pair<string, float>(rule_key, rule.GetType())); 
	}


	// fill motif group lists
	fill_CMG(neigh_ind_list);
	fill_SMG(neigh_ind_list);
	fill_SROMG(neigh_ind_list);

	// initalize system with desired SRO
	Output << "EQ passes: " << session.eq_passes << ", EQ Temp: " << session.sro_temp << "\n";
	Output << "SRO Target: " << session.sro_target << "\n";
	cout << "SRO Target: " << session.sro_target << "\n";
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

	// setup rng for spin and acceptance probability
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	std::random_device dev;
	std::mt19937 rng2(dev());
	std::uniform_int_distribution<std::mt19937::result_type> rand_atom(0, sim_cell.numb_atoms);

	// start MC loop
	cout << "entering main loop\n";
	for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
		e_avg = 0.0;
		spin_avg = 0.0;
		flip_count = 0.0;
		flip_count2 = 0.0;

		// Do passes for atom flips  
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				bool same_atom = true;
				int rand_site;
				while (same_atom == true) {
					rand_site = rand_atom(rng2);
					if (rand_site != site) {
						if (find(sim_cell.atom_list[rand_site].allowed_species.begin(), sim_cell.atom_list[rand_site].allowed_species.end(), chem_list[rand_site]) != sim_cell.atom_list[rand_site].allowed_species.end()) { same_atom = false; }
					}
				}
				// Flip atoms
				int old_site_chem = chem_list[site];
				int old_site_spin = spin_list[site];
				int old_rand_site_chem = chem_list[rand_site];
				int old_rand_site_spin = spin_list[rand_site];
				chem_list[site] = old_rand_site_chem;
				spin_list[site] = old_rand_site_spin;
				chem_list[rand_site] = old_site_chem;
				spin_list[rand_site] = old_site_spin;
				// Evaluate energy change
				e_flip = eval_atom_flip(site, rand_site);
				if (e_flip < 0) { flip_count += 1; vect_add(sro_list, sro_flip); }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; vect_add(sro_list, sro_flip); }
					else { 
						chem_list[site] = old_site_chem;
						spin_list[site] = old_site_spin;	
						chem_list[rand_site] = old_rand_site_chem;
						spin_list[rand_site] = old_rand_site_spin;
						e_flip = 0.0; 
						spin_flip = 0; }
				}
				init_enrg += e_flip;
				if (pass >= passes * .2) {
					e_avg += init_enrg;
					rs_C.Push(init_enrg);
					rs_X.Push(init_spin);
				}
			}
		}

		// Do sub_passes for spin flips
		for (int sub_pass = 0; sub_pass < sub_passes; sub_pass++) {
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
					if (sub_pass >= sub_passes * .2) {
						e_avg += init_enrg;
						rs_C.Push(init_enrg);
						spin_avg += init_spin / (pow(numb_atoms, 2) * 0.8 * sub_passes);
						rs_X.Push(init_spin);
					}
				}
			}
		}
		e_avg /= double(numb_atoms * numb_atoms * 0.8 * (sub_passes + passes));
		var_e = rs_C.Variance();
		var_spin = rs_X.Variance();
		Cmag = var_e / (Kb * double(pow(temp, 2)));
		Xmag = var_spin / (Kb * double(pow(temp, 2)));
		Output << " # "
			<< temp << ", "
			<< e_avg << ", "
			<< spin_avg << ", "
			<< var_e << ", "
			<< var_spin << ", "
			<< Cmag << ", "
			<< Xmag << ", "
			<< sro_list[0] << ", "
			<< flip_count << ", "
			<< flip_count2 << "\n";
		rs_C.Clear();
		rs_X.Clear();
	}
	cout << " MC Finished\n";
	print_state(temp2);
	Output.close();
}

vector<int> Algo3::eval_site_sro(int site) {
	vector<int> site_sro(sro_list.size(), 0);
	map<string, float>::iterator rule_itr;
	for (int i = 0; i < sro_motif_groups[site].size(); i++) {
		string rule_key = "0.";
		rule_key += to_string(i);
		for (int j : sro_motif_groups[site][i]) {
			rule_key += "." + to_string(chem_list[j]);
		}
		rule_itr = rule_map_sro.find(rule_key);
		if (rule_itr != rule_map_sro.end()) { sro_list[rule_itr->second] += 1; }
	}
	return site_sro;
}

//float Algo3::init_SRO(vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list) {
//	// setup rng
//	std::mt19937_64 rng;
//	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
//	rng.seed(ss);
//	std::uniform_real_distribution<double> unif(0, 1);
//
//	std::random_device dev;
//	std::mt19937 rng2(dev());
//	std::uniform_int_distribution<std::mt19937::result_type> rand_atom(0, sim_cell.numb_atoms);
//
//	// begin SRO MC run //
//	float sro_final = 0;
//	float sro_initial = 0;
//	float sro_site_new = 0;
//	float sro_site_old = 0;
//	float sro_site_flip = 0;
//	float sro_flip = 0;
//	float keep_rand = 0;
//	float keep_prob = 0;
//
//	for (int i = 0; i < sim_cell.numb_atoms; i++) {
//		sro_initial += calc_struct(i, neigh_ind_list, neigh_dist_list) / sim_cell.species_numbs[2];
//	}
//	if (session.eq_passes > 0) {
//		cout << "Setting SRO\n";
//		for (int i = 0; i < session.eq_passes; i++) {
//			for (int site = 0; site < sim_cell.numb_atoms; site++) {
//				// Flip Species
//				if (chem_list[site] != 0) {
//					int rand_index = site;
//					while (rand_index == site or chem_list[rand_index] == 0) { rand_index = rand_atom(rng2); }
//					int old_species_site = chem_list[site];
//					int old_species_rand = chem_list[rand_index];
//					if (old_species_site != old_species_rand) {
//						sro_site_old = calc_struct(site, neigh_ind_list, neigh_dist_list);
//						sro_site_old += calc_struct(rand_index, neigh_ind_list, neigh_dist_list);
//						chem_list[site] = old_species_rand;
//						chem_list[rand_index] = old_species_site;
//						sro_site_new = calc_struct(site, neigh_ind_list, neigh_dist_list);
//						sro_site_new += calc_struct(rand_index, neigh_ind_list, neigh_dist_list);
//						sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
//						sro_flip = abs((sro_site_flip + sro_initial) - session.sro_target) - abs(sro_initial - session.sro_target);
//						if (sro_flip < 0) { sro_initial += sro_site_flip; }
//						else {
//							keep_rand = unif(rng);
//							keep_prob = exp(-1 / (Kb * session.sro_temp) * (sro_flip));
//							if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
//							else {
//								chem_list[site] = old_species_site;
//								chem_list[rand_index] = old_species_rand;
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	for (int i = 0; i < sim_cell.numb_atoms; i++) {
//		sro_final += calc_struct(i, neigh_ind_list, neigh_dist_list) / sim_cell.species_numbs[2];
//	}
//	return sro_final;

bool Algo3::bc_check(vector<float> check_vect, vector<float>& pos) {
	bool bc_test = false;
	vector<int> dir { -1, 1 };
	vector<float> bc_pos { 0, 0, 0 };
	vector<float> lc_shift { 0, 0, 0 };
	vector<vector<float>> bc_shifts { pos };

	for (int a : dir) {
		for (int i = 0; i < 3; i++) {
			lc_shift = scale_vect(sim_cell.lat_vect[i], a);
			bc_pos = vect_add(pos, lc_shift);
			bc_shifts.push_back(bc_pos);
		}
		for (int b : dir) {
			lc_shift = scale_vect(sim_cell.lat_vect[0], a);
			bc_pos = vect_add(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[1], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			lc_shift = scale_vect(sim_cell.lat_vect[0], a);
			bc_pos = vect_add(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[2], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			lc_shift = scale_vect(sim_cell.lat_vect[1], a);
			bc_pos = vect_add(pos, lc_shift);
			lc_shift = scale_vect(sim_cell.lat_vect[2], b);
			bc_pos = vect_add(bc_pos, lc_shift);
			bc_shifts.push_back(bc_pos);

			for (int c : dir) {
				lc_shift = scale_vect(sim_cell.lat_vect[0], a);
				bc_pos = vect_add(pos, lc_shift);
				lc_shift = scale_vect(sim_cell.lat_vect[1], b);
				bc_pos = vect_add(bc_pos, lc_shift);
				lc_shift = scale_vect(sim_cell.lat_vect[2], c);
				bc_pos = vect_add(bc_pos, lc_shift);
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

void Algo3::fill_SMG(vector<vector<int>>& neigh_ind_list) {
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
						new_pos = vect_add(pos_list[i], shift);
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

void Algo3::fill_CMG(vector<vector<int>>& neigh_ind_list) {
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
						new_pos = vect_add(pos_list[i], shift);
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

void Algo3::fill_SROMG(vector<vector<int>>& neigh_ind_list) {
	vector<int> deco_group;
	vector<float> new_pos { 0.0, 0.0, 0.0 };
	vector<float> self_site { 0.0, 0.0, 0.0 };
	vector<vector<int>> motifs;
	int last_ind = -1;
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		for (SRO rule : session.sro_rule_list) {
			if (last_ind != rule.motif_ind) {
				for (vector<float> shift : rule.motif) {
					if (pos_comp(shift, self_site)) { deco_group.push_back(i); }
					else {
						new_pos = vect_add(pos_list[i], shift);
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
		sro_motif_groups.push_back(motifs);
		motifs.clear();
	}
}

void Algo3::print_state(float temp) {
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
	string temp_str = to_string(temp);
	temp_str = replace_char(temp_str, '.', '_');
	string file_name = "CONTCAR" + to_string(outfile_count) + temp_str;
	OUT_file.open("state.txt");
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
