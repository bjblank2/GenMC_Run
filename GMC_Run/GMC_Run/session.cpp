#include "session.h"

Session::Session(void) {}

Session::Session(string input_file) {
	vector<string> input_lines;
	string input_settings;
	ifstream session_input;
	vector<string> setting;
	session_input.open(input_file, ifstream::in);
	if (session_input.is_open()) {
		while (getline(session_input, input_settings))
		{
			input_lines.push_back(input_settings);
		}
		session_input.close();
	}
	else cout << "Unable to open input file\n";
	mag_ext = 0; //Default value in case it isnt listed in the INPUT file.
	for (int i = 0; i < input_lines.size(); i++) {
		setting = split(input_lines[i], " ");
		if (setting[0].compare("STRCTURE") == 0) { structure_file = setting[2]; }
		else if (setting[0].compare("ALGO") == 0) { algo = stoi(setting[2]); }
		else if (setting[0].compare("RULES_FILE") == 0) { rules_file = setting[2]; }
		else if (setting[0].compare("SIM_TYPE") == 0) { sim_type = setting[2]; }
		else if (setting[0].compare("PHASE_INIT") == 0) { phase_init = setting[2]; }
		else if (setting[0].compare("SPIN_INIT") == 0) { spin_init = setting[2]; }
		else if (setting[0].compare("SPECIES_INIT") == 0) { species_init = setting[2]; }
		else if (setting[0].compare("NUMB_PASSES") == 0) { numb_passes = stoi(setting[2]); }
		else if (setting[0].compare("NUMB_SUBPASSES") == 0) { numb_subpasses = stoi(setting[2]); }
		else if (setting[0].compare("START_TEMP") == 0) { start_temp = stof(setting[2]); }
		else if (setting[0].compare("END_TEMP") == 0) { end_temp = stof(setting[2]); }
		else if (setting[0].compare("TEMP_INC") == 0) { temp_inc = stof(setting[2]); }
		else if (setting[0].compare("EQ_PASSES") == 0) { eq_passes = stoi(setting[2]); }
		else if (setting[0].compare("SRO_TEMP") == 0) { sro_temp = stof(setting[2]); }
		else if (setting[0].compare("MAG_EXT") == 0) { mag_ext = stof(setting[2]); }
		else if (setting[0].compare("SRO_TARGET") == 0) { sro_target = stof(setting[2]); }
		else if (setting[0].compare("USE_STATES") == 0) {
			if (setting[2][0] == 'T') { use_states = true; }
			else { use_states = false; }
		}
		else if (setting[0].compare("USE_POSCAR") == 0) {
			//cout << setting[0] << "_" << setting[1] << "_" << setting[2] << "\n";
			setting[2].erase(std::remove(setting[2].begin(), setting[2].end(), '\n'), setting[2].end());
			setting[2].erase(std::remove(setting[2].begin(), setting[2].end(), ' '), setting[2].end());
			//cout << setting[0] << "_" << setting[1] << "_" << setting[2] << "\n";
			if (setting[2][0] == 'T') { use_poscar = true; }// .compare("TRUE") == 0) { use_poscar = true; }
			else { use_poscar = false; }
			//cout << "useposcar " << use_poscar << "\n";
		}
		else if (setting[0].compare("ATOM_NUMBS") == 0) {
			tot_atoms = 0;
			for (int j = 2; j < setting.size(); j++) {
				atom_numbs.push_back(stoi(setting[j]));
				tot_atoms += stoi(setting[j]);
			}
		}
		else if (setting[0].compare("SPECIES") == 0) {
			for (int j = 2; j < setting.size(); j++) {
				string str = setting[j];
				str.erase(remove(str.begin(), str.end(), '\r'), str.end());
				str.erase(remove(str.begin(), str.end(), ' '), str.end());
				species_str.push_back(str);
				species_inds.push_back(j - 2);
			}
		}
		else if (setting[0].compare("SHAPE") == 0) {
			shape[0] = stoi(setting[2]);
			shape[1] = stoi(setting[3]);
			shape[2] = stoi(setting[4]);
		}
		else {
			cout << "Unrecognised input flag: " << setting[0] << "\n";
		}
		rules_file = "Rule_file.in";
	}
	if (moments.size() == 0) {
		moments.clear();
		for (int i = 0; i < atom_numbs.size(); i++) {
			moments.push_back(1);
		}
	}
	if (species_inds.size() != atom_numbs.size()) { 
		cout << "*ERROR* number of simulated atoms is unclear\n" << " SPECIES should reflect ATOM_NUMBS\n"; 
	}
}

Session::Session(Session& _session) {
	algo = _session.algo;
	structure_file = _session.structure_file;
	rules_file = _session.rules_file;
	sim_type = _session.sim_type;
	phase_init = _session.phase_init;
	spin_init = _session.spin_init;
	species_init = _session.species_init;
	use_poscar = _session.use_poscar;
	numb_passes = _session.numb_passes;
	numb_subpasses = _session.numb_subpasses;
	eq_passes = _session.eq_passes;
	start_temp = _session.start_temp;
	end_temp = _session.end_temp;
	temp_inc = _session.temp_inc;
	sro_temp = _session.sro_temp;
	sro_target = _session.sro_target;
	mag_ext = _session.mag_ext;
	copy(_session.shape, _session.shape + 3, shape);
	atom_numbs = _session.atom_numbs;
	use_states = _session.use_states;
	spin_states = _session.spin_states;
	moments = _session.moments;
}

void Session::_copy(Session& _session) {
	algo = _session.algo;
	structure_file = _session.structure_file;
	rules_file = _session.rules_file;
	sim_type = _session.sim_type;
	phase_init = _session.phase_init;
	spin_init = _session.spin_init;
	species_init = _session.species_init;
	use_poscar = _session.use_poscar;
	numb_passes = _session.numb_passes;
	numb_subpasses = _session.numb_subpasses;
	eq_passes = _session.eq_passes;
	start_temp = _session.start_temp;
	end_temp = _session.end_temp;
	temp_inc = _session.temp_inc;
	sro_temp = _session.sro_temp;
	sro_target = _session.sro_target;
	mag_ext = _session.mag_ext;
	copy(_session.shape, _session.shape + 3, shape);
	atom_numbs = _session.atom_numbs;
	use_states = _session.use_states;
	spin_states = _session.spin_states;
	moments = _session.moments;
}

void Session::add_spin_states(string input_file) {
	vector<string> input_lines;
	string input_settings;
	ifstream session_input;
	vector<string> setting;
	vector<float> atom_spin_states;
	session_input.open(input_file, ifstream::in);
	if (session_input.is_open()) {
		while (getline(session_input, input_settings))
		{
			input_lines.push_back(input_settings);
		}
		session_input.close();
	}
	else cout << "Unable to open input file\n";
	for (int i = 0; i < input_lines.size(); i++) {
		setting = split(input_lines[i], " ");
		for (int j = 0; j < setting.size(); j++) { atom_spin_states.push_back(stof(setting[j])); }
		spin_states.push_back(atom_spin_states);
		atom_spin_states.clear();
	}
	if (spin_states.size() != atom_numbs.size()) { cout << "Error: atom_numbs not equal to number of lines in spin_states file\n"; }
}

void Session::fill_sro_list() {
	vector<float> distances;
	vector<int> species;
	vector<float> type;
	string sro_line;
	vector<string> sro_lines;
	vector<string> setting;
	vector<string> line;
	vector<vector<float>> motif;
	vector<vector<int>> deco;
	ifstream sro_list_file;
	sro_file = "SRO.in";
	sro_list_file.open(sro_file, ifstream::in);
	// Parce rule list txt file
	if (sro_list_file.is_open()) {
		while (getline(sro_list_file, sro_line))
		{
			sro_lines.push_back(sro_line);
		}
		sro_list_file.close();
		cout << "Read sro file\n";
		int sro_motif_ind = 0;
		for (int i = 0; i < sro_lines.size(); i++) {
			if (sro_lines[i].find('#') != std::string::npos) {
				for (int j = 0; j < deco.size(); j++) {
					sro_rule_list.push_back(SRO(type[j], deco[j], motif, sro_motif_ind));
				}
				motif.clear();
				deco.clear();
				type.clear();
				sro_motif_ind +=1;
			}
			setting = split(sro_lines[i], ":");
			if (setting[0].find("otif") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				for (int j = 0; j < line.size(); j++) {
					vector<string> pos = split(line[j], ",");
					motif.push_back({ stof(pos[0]), stof(pos[1]), stof(pos[2]) });
				}

			}
			else if (setting[0].find("eco") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				for (int j = 0; j < line.size(); j++) {
					vector<string> specs = split(line[j], ",");
					vector<int>spec = {};
					for (int k = 0; k < specs.size(); k++) {
						spec.push_back(stoi(specs[k]));
					}
					deco.push_back(spec);
				}
			}
			else if (setting[0].find("ype") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				for (int j = 0; j < line.size(); j++) {
					type.push_back(stof(line[j]));
				}
			}
		}
		cout << "Filled Rule_list\n";
	}
	else cout << "*ERROR* Unable to open rule file\n";

}

void Session::fill_rule_list(){
	vector<float> distances;
	vector<int> spins;
	vector<int> species;
	float energy_contribution = 0;
	int rule_type = 0;
	int rule_length = 0;
	string rule_line;
	vector<string> rule_lines;
	vector<string> setting;
	vector<string> line;
	int type;
	vector<float> enrg;
	vector<int> phase;
	vector<vector<float>> motif;
	vector<vector<int>> deco;
	ifstream rule_list_file;
	rule_list_file.open(rules_file, ifstream::in);
	// Parce rule list txt file
	if (rule_list_file.is_open()) {
		while (getline(rule_list_file, rule_line))
		{
			rule_lines.push_back(rule_line);
		}
		rule_list_file.close();
		cout << "Reading rule file\n";
		int chem_motif_ind = 0;
		int spin_motif_ind = 0;
		bool intercept_flag = false;
		for (int i = 0; i < rule_lines.size(); i++) {
			if (rule_lines[i].find('#') != std::string::npos and intercept_flag == false) {
				if (phase.size() == 0) { for (int j = 0; j < enrg.size(); j++) { phase.push_back(0); } }
				//cout << phase.size() << "  " << type.size() << "\n";
				for (int j = 0; j < enrg.size(); j++) {
					if (type == 0 and enrg[j] != 0.0) { chem_rule_list.push_back(Rule(enrg[j], type, phase[j], deco[j], motif, chem_motif_ind)); }
					if (type == 1 and enrg[j] != 0.0) { spin_rule_list.push_back(Rule(enrg[j], type, phase[j], deco[j], motif, spin_motif_ind)); }
				}
				motif.clear();
				deco.clear();
				if (type == 0) { chem_motif_ind += 1; }
				else if (type == 1) { spin_motif_ind += 1; }
				type=0;
				enrg.clear();
				phase.clear();
			}
			setting = split(rule_lines[i], ":");
			if (setting[0].find("Motif") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				if (line[0].find("nter") != std::string::npos) {
					intercept_flag = true;
				}
				else {
					intercept_flag = false;
					for (int j = 0; j < line.size(); j++) {
						vector<string> pos = split(line[j], ",");
						motif.push_back({ stof(pos[0]), stof(pos[1]), stof(pos[2]) });
					}
				}
			}
			else if (setting[0].find("Deco") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				for (int j = 0; j < line.size(); j++) {
					vector<string> specs = split(line[j], ",");
					vector<int>spec = {};
					for (int k = 0; k < specs.size(); k++) {
						spec.push_back(stoi(specs[k]));
					}
					deco.push_back(spec);
				}
			}
			else if (setting[0].find("Type") != std::string::npos) {
					type = stoi(setting[1]);
			}
			else if (setting[0].find("Enrg") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				if (intercept_flag == true) { intercept = stof(line[0]); }
				else {
					for (int j = 0; j < line.size(); j++) {
						enrg.push_back(stof(line[j]));
					}
				}
			}
			else if (setting[0].find("Phase") != std::string::npos) {
				setting.erase(setting.begin());
				line = setting;
				for (int j = 0; j < line.size(); j++) {
					phase.push_back(stoi(line[j]));
				}
			}
		}
		cout << "Rule list filled\n";
	}
	else cout << "*ERROR* Unable to open rule file\n";
	
}

void Session::find_unique_dists(){
	// Loop through all mc_rules
	for (int i = 0; i < chem_rule_list.size(); i++) {
		for (int j = 0; j < chem_rule_list[i].GetDists().size(); j++) {
			if (find(unique_dists.begin(), unique_dists.end(), chem_rule_list[i].GetDists()[j]) == unique_dists.end()) {
				unique_dists.push_back(chem_rule_list[i].GetDists()[j]);
			}
		}
	}
	for (int i = 0; i < spin_rule_list.size(); i++) {
		for (int j = 0; j < spin_rule_list[i].GetDists().size(); j++) {
			if (find(unique_dists.begin(), unique_dists.end(), spin_rule_list[i].GetDists()[j]) == unique_dists.end()) {
				unique_dists.push_back(spin_rule_list[i].GetDists()[j]);
			}
		}
	}
}
