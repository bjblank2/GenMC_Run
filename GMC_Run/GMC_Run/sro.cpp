#include "sro.h"


// This is the SRO class (mc_SROs object). It servs as a container for the MC SROs and is used to create the SRO map
SRO::SRO(void) {
}

SRO::SRO( int _type, vector<int> _deco, vector<vector<float>> _motif, int _motif_ind) {
	motif = _motif;
	type = _type;
	deco = _deco;
	motif_ind = _motif_ind;
	for (int i = 0; i < motif.size(); i++) {
		for (int j = i + 1; j < motif.size(); j++) {
			float new_dist = 0;
			for (int k = 0; k < 3; k++) { new_dist += pow((motif[i][k] - motif[j][k]), 2.0); }
			new_dist = sqrt(new_dist);
			if (std::find(dists.begin(), dists.end(), new_dist) == dists.end()) {
				dists.push_back(new_dist);
			}
		}
	}
}


int SRO::GetType() {
	return type;
}

int SRO::GetLength() {
	return motif.size();
}

vector<int> SRO::GetDeco() {
	return deco;
}

vector<float> SRO::GetDists() {
	return dists;
}
