#include "rule.h"


// This is the rule class (mc_rules object). It servs as a container for the MC rules and is used to create the rule map
Rule::Rule(void) {
}

Rule::Rule(long double _enrg_cont, int _type, int _phase, vector<int> _species, vector<vector<float>> _motif) {
	motif = _motif;
	enrg_cont = _enrg_cont;
	type = _type;
	phase = _phase;
	species = _species;
	for (int i = 0; i < motif.size(); i++) {
		for (int j = i + 1; j < motif.size(); j++) {
			float new_dist = 0;
			for (int k = 0; k < 3; k++) {
				new_dist += pow((motif[i][k] - motif[j][k]), 2.0);
				new_dist = sqrt(new_dist);
				if (std::find(dists.begin(), dists.end(), new_dist) == dists.end()) {
					dists.push_back(new_dist);
				}
			}
		}
	}
}

int Rule::GetPhase() {
	return phase;
}
long double Rule::GetEnrgCont() {
	return enrg_cont;
}
int Rule::GetType() {
	return type;
}
int Rule::GetLength() {
	return species.size();
}
vector<int> Rule::GetSpecies() {
	return species;
}
vector<float> Rule::GetDists() {
	return dists;
}

bool Rule::IsRuleChem(vector<int> test_species, vector<float> test_dists) {
	if (type != 0) { return false; }
	else if (test_dists.size() != dists.size() or test_species.size() != test_species.size()) { return false; }
	else {
		for (int i = 0; i < test_dists.size(); i++) { test_dists[i] = round(test_dists[i] * 100000) / 100000; }
		for (int i = 0; i < dists.size(); i++) { dists[i] = round(dists[i] * 100000) / 100000; }
		if (species.size() == 1) {
			if (test_species[0] != species[0]) { return false; }
			else { return true; }
		}
		else if (species.size() == 2) {
			if (test_dists[0] == dists[0]) {
				if (test_species[0] == species[0] and test_species[1] == species[1]) { return true; }
				else if (test_species[0] == species[1] and test_species[1] == species[0]) { return true; }
				else { return false; }
			}
			else { return false; }
		}
		else if (species.size() == 3) {
			if (test_species[0] == species[0] and test_species[1] == species[1] and test_species[2] == species[2] and test_dists[0] == dists[0] and test_dists[1] == dists[1] and test_dists[2] == dists[2]) { return true; }
			else if (test_species[0] == species[0] and test_species[2] == species[1] and test_species[1] == species[2] and test_dists[1] == dists[0] and test_dists[0] == dists[1] and test_dists[2] == dists[2]) { return true; }
			else if (test_species[1] == species[0] and test_species[0] == species[1] and test_species[2] == species[2] and test_dists[0] == dists[0] and test_dists[2] == dists[1] and test_dists[1] == dists[2]) { return true; }
			else if (test_species[1] == species[0] and test_species[2] == species[1] and test_species[0] == species[2] and test_dists[2] == dists[0] and test_dists[0] == dists[1] and test_dists[1] == dists[2]) { return true; }
			else if (test_species[2] == species[0] and test_species[0] == species[1] and test_species[1] == species[2] and test_dists[1] == dists[0] and test_dists[2] == dists[1] and test_dists[0] == dists[2]) { return true; }
			else if (test_species[2] == species[0] and test_species[1] == species[1] and test_species[0] == species[2] and test_dists[2] == dists[0] and test_dists[1] == dists[1] and test_dists[0] == dists[2]) { return true; }
			else { return false; }
		}
	}
}

bool Rule::IsRuleSpin(vector<int> test_species, vector<float> test_dists) {
	if (type != 1) { return false; }
	else if (test_dists.size() != dists.size() or test_species.size() != test_species.size()) { return false; }
	else {
		if (species.size() == 1) {
			if (test_species[0] != species[0]) { return false; }
			else { return true; }
		}
		else if (species.size() == 2) {
			if (test_dists[0] == dists[0]) {
				if (test_species[0] == species[0] and test_species[1] == species[1]) { return true; }
				else if (test_species[0] == species[1] and test_species[1] == species[0]) { return true; }
				else { return false; }
			}
			else { return false; }
		}
		else if (species.size() == 3) {
			if (test_species[0] == species[0] and test_species[1] == species[1] and test_species[2] == species[2] and test_dists[0] == dists[0] and test_dists[1] == dists[1] and test_dists[2] == dists[2]) { return true; }
			else { return false; }
		}
	}
}
