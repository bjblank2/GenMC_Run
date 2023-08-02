#include "algo3.h"

Algo3::Algo3(void) {}

Algo3::Algo3(Session& _session, SimCell& _sim_cell) {
	session = _session;
	sim_cell = _sim_cell;
	if (session.numb_passes < 1) {
		cout << "_______________________________________________________________________________" << endl;
		cout << "Possible Error: Algo3 has been given 0 passes per main pass" << endl;
		cout << "This implies that no seperate spin flips are made" << endl;
		cout << "This is probably not what you want..." << endl;
		cout << "_______________________________________________________________________________" << endl;
	}
}

float Algo3::eval_site_chem(int site) {
    float enrg = 0;
    map<string, vector<float>>::iterator rule_itr;
    for (int i = 0; i < chem_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = chem_motif_groups[site][i];
        for (int j =0; j < motif.size(); j++) {
            string rule_key = "0."; // chem ind
            rule_key += to_string(i); // clust_ind
            vector<int> group = motif[j];
            for ( int k : group ) {
                rule_key += "." + to_string(chem_list[k]); // sites ind
            }
            rule_itr = rule_map_chem.find(rule_key);
            if (rule_itr != rule_map_chem.end()) {
                enrg += rule_itr->second[0] / group.size();
                lat_rule_count_list[round(rule_itr->second[1])] += 1.0 / group.size();
            }
        }
    }
    return enrg;
}

float Algo3::eval_site_spin(int site) {
    float enrg = 0;
    map<string, float>::iterator rule_itr;
    for (int i = 0; i < spin_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = spin_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            string rule_key = "1."; // spin ind
            rule_key += to_string(i); // clust_ind
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_key += "." + to_string(chem_list[k]); // sites ind
                spin_prod *= spin_list[k];
            }
            rule_itr = rule_map_spin.find(rule_key);
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod / group.size()) : 0.0;
        }
    }
    return enrg;
}

float Algo3::eval_spin_flip(int site, float old_spin) {
    float enrg = 0;
    map<string, float>::iterator rule_itr;
    for (int i = 0; i < spin_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = spin_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            string rule_key = "1.";
            rule_key += to_string(i);
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_key += "." + to_string(chem_list[k]);
                if (k != site) { spin_prod *= spin_list[k]; }
            }
            rule_itr = rule_map_spin.find(rule_key);
            enrg += (rule_itr != rule_map_spin.end()) ? (rule_itr->second * spin_prod) : 0.0;
        }
    }
    return (enrg * spin_list[site] - enrg * old_spin);
}

float Algo3::eval_atom_flip(int site) {
    map<string, vector<float>>::iterator rule_itr_chem;
    map<string, float>::iterator rule_itr_spin;
    float enrg = 0.0;
    fill(site_rule_count_list.begin(), site_rule_count_list.end(), 0);
    for (int i = 0; i < chem_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = chem_motif_groups[site][i];
        for (int j =0; j < motif.size(); j++) {
            string rule_key = "0."; // chem ind
            rule_key += to_string(i); // clust_ind
            vector<int> group = motif[j];
            for (int k : group ) {
                rule_key += "." + to_string(chem_list[k]); // sites ind
            }
            rule_itr_chem = rule_map_chem.find(rule_key);
            if (rule_itr_chem != rule_map_chem.end()) {
                enrg += rule_itr_chem->second[0];
                site_rule_count_list[round(rule_itr_chem->second[1])] += 1.0;
            }
        }
    }
    for (int i = 0; i < spin_motif_groups[site].size(); i++) {
        vector<vector<int>> motif = spin_motif_groups[site][i];
        for (int j = 0; j < motif.size(); j++) {
            string rule_key = "1."; // spin ind
            rule_key += to_string(i); // clust_ind
            vector<int> group = motif[j];
            float spin_prod = 1;
            for (int k : group) {
                rule_key += "." + to_string(chem_list[k]); // sites ind
                spin_prod *= spin_list[k];
            }
            rule_itr_spin = rule_map_spin.find(rule_key);
            enrg += (rule_itr_spin != rule_map_spin.end()) ? (rule_itr_spin->second * spin_prod ) : 0.0;
        }
    }
    return enrg;
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
    vector<float> new_pos { 0.0, 0.0, 0.0 };
    vector<float> self_site { 0.0, 0.0, 0.0 };
    vector<int> sites;
    vector<vector<int>> groups;
    vector<vector<vector<int>>> motifs;
    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
        for (vector<vector<vector<float>>> clust : session.spin_motif_list) {
            for (vector<vector<float>> motif : clust) {
                for (vector<float> shift : motif) {
                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
                    else {
                        new_pos = vect_add(pos_list[atom], shift);
                        for (int neigh : neigh_ind_list[atom]) {
                            if (bc_check(pos_list[neigh], new_pos)) {
                                sites.push_back(neigh);
                            }
                        }
                    }
                }
                groups.push_back(sites);
                sites.clear();
            }
            motifs.push_back(groups);
            groups.clear();
        }
        spin_motif_groups.push_back(motifs);
        motifs.clear();
    }
}

void Algo3::fill_CMG(vector<vector<int>>& neigh_ind_list) {
    vector<float> new_pos { 0.0, 0.0, 0.0 };
    vector<float> self_site { 0.0, 0.0, 0.0 };
    vector<int> sites;
    vector<vector<int>> groups;
    vector<vector<vector<int>>> motifs;
    for (int atom = 0; atom < sim_cell.numb_atoms; atom++) {
        for (vector<vector<vector<float>>> clust : session.chem_motif_list) {
            for (vector<vector<float>> motif : clust) {
                for (vector<float> shift : motif) {
                    if (pos_comp(shift, self_site)) { sites.push_back(atom); }
                    else {
                        new_pos = vect_add(pos_list[atom], shift);
                        for (int neigh : neigh_ind_list[atom]) {
                            if (bc_check(pos_list[neigh], new_pos)) {
                                sites.push_back(neigh);
                            }
                        }
                    }
                }
                groups.push_back(sites);
                sites.clear();
            }
            motifs.push_back(groups);
            groups.clear();
        }
        chem_motif_groups.push_back(motifs);
        motifs.clear();
    }
}

void Algo3::print_state(string contcar_name, int temp) {
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
    string file_name = contcar_name + "-" + to_string(temp);
    OUT_file.open(file_name);
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

void Algo3::run() {
    // declare variables
    int rand_site = 1;
    float rand_spin = 0.0;
    bool same_spin;
    bool same_atom;
    float Cmag = 0.0;
    float Xmag = 0.0;
    int passes = session.numb_passes;
    int eq_passes = session.eq_passes;
    float temp1 = session.start_temp;
    float temp2 = session.end_temp;
    float temp_inc = session.temp_inc;
    float keep_rand;
    float keep_prob;
    vector<vector<float>> spin_states = session.spin_states;
    int numb_atoms = sim_cell.numb_atoms;
    int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
    vector<vector<int>> neigh_ind_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
    vector<int> spin_atoms; // atomic species that can have spin

    // Create seperate output file to avoid race condition
    string file_name = "OUTPUT";
    string sro_file_name = "OUTPUT_SRO";
    string contcar_name = "CONTCAR";
    bool file_exists = true;
    while (file_exists == true) {
        const char* c_file = file_name.c_str();
        int fd = open(c_file, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
        if (fd < 0) {
            // file exists or otherwise uncreatable
            outfile_count += 1;
            file_name = "OUTPUT" + to_string(outfile_count);
            sro_file_name = "OUTPUT_SRO" + to_string(outfile_count);
            contcar_name = "CONTCAR" + to_string(outfile_count);
        }
        else {
            file_exists = false;
            close(fd);
        }
    }
    const char* c_file = file_name.c_str();
    const char* sro_file = sro_file_name.c_str();
    ofstream Output;
    ofstream SROout;
    Output.open(c_file);
    SROout.open(sro_file);
    
    // Output energy and spin for convergence test
    ofstream Output_converge;
    if ( session.do_conv_output) { Output_converge.open("OUTPUT_CONVERG"); }

    Output << "Phase: " << sim_cell.phase_init;
    Output << "Composition: ";
    for (int i = 0; i < sim_cell.species_numbs.size(); i++) { Output << sim_cell.species_numbs[i] << ", "; }
    Output << "\n";    Output << "MC passes: " << session.numb_passes << ", ";
    Output << "Beginning MC EQ run using Algo3\n";
    
    cout << "Making atom list and neighbor index list\n";
    // Make atom_list more acessable for species and spin and neighbors
    for (int i = 0; i < sim_cell.numb_atoms; i++) {
        chem_list.push_back(sim_cell.atom_list[i].getSpecies());
        spin_list.push_back(sim_cell.atom_list[i].getSpin());
        pos_list.push_back({ sim_cell.atom_list[i].pos[0], sim_cell.atom_list[i].pos[1], sim_cell.atom_list[i].pos[2] });
        for (int j = 0; j < numb_neighbors; j++) {
            neigh_ind_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
        }
    }
    
    cout << "Making rule maps\n";
    // Make rule_maps for easy lookup and initialize rule_count_list
    float ind = 0;
    string rule_key;
    map<string, vector<float>>::iterator rule_itr;
    for (Rule rule : session.chem_rule_list) {
        rule_key = to_string(rule.GetType()) + "." + to_string(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_key += "." + to_string(rule.deco[i]); }
        rule_itr =
        rule_itr = rule_map_chem.find(rule_key);
        if (rule_itr == rule_map_chem.end()) {
            vector<float> enrg_ind_pair = { rule.GetEnrgCont(), ind };
            rule_map_chem.insert(pair<string, vector<float>>(rule_key, enrg_ind_pair));
            lat_rule_count_list.push_back(0);
            site_rule_count_list.push_back(0);
            ind += 1;
        }
    }
    for (Rule rule : session.spin_rule_list) {
        rule_key = to_string(rule.GetType()) + "." + to_string(rule.clust_ind);
        for (int i = 0; i < rule.deco.size(); i++) { rule_key += "." + to_string(rule.deco[i]); }
        rule_map_spin.insert(pair<string, float>(rule_key, rule.GetEnrgCont()));
        for (int atom : rule.deco) {
            if (find(spin_atoms.begin(), spin_atoms.end(), atom) == spin_atoms.end()) { spin_atoms.push_back(atom); } // initialize spin_atoms
        }
    }
    
    // Fill motif group lists
    cout << "Making chem motif group lists\n";
    fill_CMG(neigh_ind_list);
    cout << "Making spin motif group lists\n";
    fill_SMG(neigh_ind_list);
	
    // Begin MC
    float init_enrg = eval_lat();
    vector<float> init_sro = lat_rule_count_list; // record initial SRO/rule count list
    cout << "Evaluated total energy: " << init_enrg / numb_atoms << "\n";
    float init_spin_enrg = eval_lat_spin();
    cout << "Evaluated spin energy: " << init_spin_enrg / numb_atoms << "\n";
    Output << "initial total energy, spin energy\n";
    Output << init_enrg / numb_atoms << ", " << init_spin_enrg / numb_atoms << "\n";
    Output << "temp, enrg, mag, var_e, var_spin, Cmag, Xmag, flip_count, flip_count2 \n";
    float init_spin = 0.0;
    float var_spin = 0.0;
    float var_e = 0.0;
    RunningStat rs_C;
    RunningStat rs_X;
    cout << "Counting spins...\n";
    for (int site = 0; site < numb_atoms; site++) {
        if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
            init_spin += spin_list[site];
        }
    }
    cout << "Initial spin is " << init_spin / numb_atoms << " per atom\n";
    float inc_dir = 1;
    if (signbit(temp2 - temp1) == 1) { inc_dir = -1; }

    // setup rng for random spin choice and acceptance probability
    std::mt19937_64 rng;
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);
    std::uniform_int_distribution<int> rand_atom(0, sim_cell.numb_atoms - 1);
    std::uniform_int_distribution<int> rand_method(0, passes - 1);
    
    // Start MC loop
    cout << "Entering main loop\n";
    int temp_count = 0;
    for (float temp = temp1; (temp2 - temp) * inc_dir >= 0; temp += temp_inc) {
        float e_avg = 0.0;
        float spin_avg = 0.0;
        int flip_count = 0;
        int flip_count2 = 0;
        vector<float> count_avg = lat_rule_count_list;
        fill(count_avg.begin(), count_avg.end(), 0.0);
        for (int pass = 0; pass < passes; pass++) {
            for (int site = 0; site < numb_atoms; site++) {
                float e_flip = 0.0;
                float spin_flip = 0.0;
                vector<float> sro_flip;
                vector<float> sro_flip1;
                vector<float> sro_flip2;
                int method_index = rand_method(rng);
                // Do the pass for spin flips
                if ( method_index < passes * 0.33) {
                    if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
                        // Flip Spin
                        float old_spin = spin_list[site];
                        float new_spin = 0.0;
                        same_spin = true;
                        while (same_spin == true) {
                            rand_spin = unif(rng);
                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[site]].size())) {
                                    new_spin = spin_states[chem_list[site]][it_spin_state]; }
                            }
                            if (new_spin != old_spin) { same_spin = false; }
                        }
                        spin_list[site] = new_spin;
                        e_flip = eval_spin_flip(site, old_spin);
                        spin_flip = new_spin - old_spin;
                        if (e_flip < 0) { flip_count += 1; }
                        else {
                            keep_rand = unif(rng);
                            keep_prob = exp(-1 / (Kb * temp) * (e_flip));
                            if (keep_rand < keep_prob) { flip_count2 += 1; }
                            else {spin_list[site] = old_spin; e_flip = 0.0; spin_flip = 0.0; }
                        }
                        // Record the enrg and spin changes
                        init_enrg += e_flip;
                        init_spin += spin_flip;
                        if (session.do_conv_output ) {
                            Output_converge << "method1 " << eval_lat() << "; " << init_enrg << ", " << e_flip << "; " << init_spin << ", " << spin_flip << "\n";
                        }
                        if (pass >= passes * 0.5) {
                            e_avg += init_enrg;
                            rs_C.Push(init_enrg);
                            spin_avg += init_spin;
                            rs_X.Push(init_spin);
                            count_avg = vect_add(count_avg, init_sro);
                        }
                    }
                }
                // Do the pass for atom swaps
                else if (method_index < passes * 0.67) {
                    same_atom = true;
                    while (same_atom == true) {
                        rand_site = rand_atom(rng);
                        if (rand_site != site) {
                            if (find(sim_cell.atom_list[rand_site].allowed_species.begin(), sim_cell.atom_list[rand_site].allowed_species.end(), chem_list[site]) != sim_cell.atom_list[rand_site].allowed_species.end() && find(sim_cell.atom_list[site].allowed_species.begin(), sim_cell.atom_list[site].allowed_species.end(), chem_list[rand_site]) != sim_cell.atom_list[site].allowed_species.end()) {
                                if (chem_list[site] != chem_list[rand_site]) { same_atom = false;}
                            }
                        }
                    }
                    int old_site_chem = chem_list[site];
                    float old_site_spin = spin_list[site];
                    int old_rand_site_chem = chem_list[rand_site];
                    float old_rand_site_spin = spin_list[rand_site];
                    // Flip atom for site
                    float old_enrg = eval_atom_flip(site);
                    vector<float> old_sro = site_rule_count_list;
                    chem_list[site] = old_rand_site_chem;
                    spin_list[site] = old_rand_site_spin;
                    float new_enrg = eval_atom_flip(site);
                    vector<float> new_sro = site_rule_count_list;
                    e_flip += new_enrg - old_enrg;
                    sro_flip1 = vect_subtract(new_sro, old_sro);
                    // Flip atom for rand_site
                    old_enrg = eval_atom_flip(rand_site);
                    old_sro = site_rule_count_list;
                    chem_list[rand_site] = old_site_chem;
                    spin_list[rand_site] = old_site_spin;
                    new_enrg = eval_atom_flip(rand_site);
                    new_sro = site_rule_count_list;
                    e_flip += new_enrg - old_enrg;
                    sro_flip2 = vect_subtract(new_sro, old_sro);
                    sro_flip = vect_add(sro_flip1, sro_flip2);
                    if (e_flip < 0) { flip_count += 1; }
                    else {
                        keep_rand = unif(rng);
                        keep_prob = exp(-1 / (Kb * temp) * (e_flip));
                        if (keep_rand < keep_prob) { flip_count2 += 1; }
                        else {
                            chem_list[site] = old_site_chem;
                            spin_list[site] = old_site_spin;
                            chem_list[rand_site] = old_rand_site_chem;
                            spin_list[rand_site] = old_rand_site_spin;
                            e_flip = 0.0;
                            fill(sro_flip.begin(), sro_flip.end(), 0.0);
                        }
                    }
                    // Record the enrg and spin changes
                    init_enrg += e_flip;
                    init_spin += spin_flip;
                    init_sro = vect_add(init_sro, sro_flip);
                    if (session.do_conv_output ) {
                        Output_converge << "method2 " << eval_lat() << "; " << init_enrg << ", " << e_flip << "; " << init_spin << ", " << spin_flip << "\n";
                    }
                    if (pass >= passes * 0.5) {
                        e_avg += init_enrg;
                        rs_C.Push(init_enrg);
                        spin_avg += init_spin;
                        rs_X.Push(init_spin);
                        count_avg = vect_add(count_avg, init_sro);
                    }
                }
                // Do the pass for atom swaps and randomly set spins
                else {
                    same_atom = true;
                    while (same_atom == true) {
                        rand_site = rand_atom(rng);
                        if (rand_site != site) {
                            if (find(sim_cell.atom_list[rand_site].allowed_species.begin(), sim_cell.atom_list[rand_site].allowed_species.end(), chem_list[site]) != sim_cell.atom_list[rand_site].allowed_species.end() && find(sim_cell.atom_list[site].allowed_species.begin(), sim_cell.atom_list[site].allowed_species.end(), chem_list[rand_site]) != sim_cell.atom_list[site].allowed_species.end()) {
                                if (chem_list[site] != chem_list[rand_site]) { same_atom = false;}
                            }
                        }
                    }
                    int old_site_chem = chem_list[site];
                    float old_site_spin = spin_list[site];
                    int old_rand_site_chem = chem_list[rand_site];
                    float old_rand_site_spin = spin_list[rand_site];
                    // Flip atom for site
                    float old_enrg = eval_atom_flip(site);
                    vector<float> old_sro = site_rule_count_list;
                    chem_list[site] = old_rand_site_chem;
                    spin_list[site] = old_rand_site_spin;
                    // Flip spin for site
                    if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[site]) != spin_atoms.end()) {
                        float old_spin1 = spin_list[site];
                        float new_spin1 = spin_list[site];
                        same_spin = true;
                        while (same_spin == true) {
                            rand_spin = unif(rng);
                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[site]].size(); it_spin_state++) {
                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[site]].size())) {
                                    new_spin1 = spin_states[chem_list[site]][it_spin_state]; }
                            }
                            if (new_spin1 != old_spin1) { same_spin = false; }
                        }
                        spin_list[site] = new_spin1;
                        spin_flip += new_spin1 - old_spin1;
                    }
                    float new_enrg = eval_atom_flip(site);
                    vector<float> new_sro = site_rule_count_list;
                    e_flip += new_enrg - old_enrg;
                    sro_flip1 = vect_subtract(new_sro, old_sro);
                    // Flip atom for rand_site
                    old_enrg = eval_atom_flip(rand_site);
                    old_sro = site_rule_count_list;
                    chem_list[rand_site] = old_site_chem;
                    spin_list[rand_site] = old_site_spin;
                    // Flip spin for rand_site
                    if (find(spin_atoms.begin(), spin_atoms.end(), chem_list[rand_site]) != spin_atoms.end()) {
                        float old_spin2 = spin_list[rand_site];
                        float new_spin2 = spin_list[rand_site];
                        same_spin = true;
                        while (same_spin == true) {
                            rand_spin = unif(rng);
                            for (int it_spin_state = 0; it_spin_state < spin_states[chem_list[rand_site]].size(); it_spin_state++) {
                                if (rand_spin > float(it_spin_state) / float(spin_states[chem_list[rand_site]].size())) {
                                    new_spin2 = spin_states[chem_list[rand_site]][it_spin_state]; }
                            }
                            if (new_spin2 != old_spin2) { same_spin = false; }
                        }
                        spin_list[rand_site] = new_spin2;
                        spin_flip += new_spin2 - old_spin2;
                    }
                    new_enrg = eval_atom_flip(rand_site);
                    new_sro = site_rule_count_list;
                    e_flip += new_enrg - old_enrg;
                    sro_flip2 = vect_subtract(new_sro, old_sro);
                    sro_flip = vect_add(sro_flip1, sro_flip2);
                    if (e_flip < 0) { flip_count += 1; }
                    else {
                        keep_rand = unif(rng);
                        keep_prob = exp(-1 / (Kb * temp) * (e_flip));
                        if (keep_rand < keep_prob) { flip_count2 += 1; }
                        else {
                            chem_list[site] = old_site_chem;
                            spin_list[site] = old_site_spin;
                            chem_list[rand_site] = old_rand_site_chem;
                            spin_list[rand_site] = old_rand_site_spin;
                            e_flip = 0.0;
                            spin_flip = 0.0;
                            fill(sro_flip.begin(), sro_flip.end(), 0.0);
                        }
                    }
                    // Record the enrg and spin changes
                    init_enrg += e_flip;
                    init_spin += spin_flip;
                    init_sro = vect_add(init_sro, sro_flip);
                    if (session.do_conv_output ) {
                        Output_converge << "method3 " << eval_lat() << "; " << init_enrg << ", " << e_flip << "; " << init_spin << ", " << spin_flip << "\n";
                    }
                    if (pass >= passes * 0.5) {
                        e_avg += init_enrg;
                        rs_C.Push(init_enrg);
                        spin_avg += init_spin;
                        rs_X.Push(init_spin);
                        count_avg = vect_add(count_avg, init_sro);
                    }
                }
            }
        }
        double scale = 1/(pow(numb_atoms, 2) * 0.5 * passes);
        e_avg *= scale;
        spin_avg *= scale;
        count_avg = scale_vect(count_avg, scale * numb_atoms);
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
            << flip_count << ", "
            << flip_count2 << "\n";
        rs_C.Clear();
        rs_X.Clear();
        SROout << "\n" << temp << " ";
        for (float x : count_avg) { SROout << x << " "; }
        SROout << "\n";
        print_state(contcar_name, temp_count);
        temp_count += 1;
    }
    cout << " MC Finished\n";
    Output.close();
    SROout.close();
}
