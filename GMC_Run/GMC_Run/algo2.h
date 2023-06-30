#pragma once
#ifndef algo2_h
#define algo2_h

#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <cmath>
#include <random>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "rule.h"
#include "sim_cell.h"
#include "running_stat.h"
#include "session.h"
#include "utils.h"
using namespace std;

class Algo2 {
public:
	int outfile_count = 0;
	SimCell sim_cell;
	Session session;
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = 0.00005788381806; // Bohr magnaton
	vector<int> chem_list;
	vector<float> spin_list;
	vector<vector<float>> pos_list;
	vector<vector<vector<int>>> spin_motif_groups;
	vector<vector<vector<int>>> chem_motif_groups;
	map <string, float> rule_map_chem;
	map <string, float> rule_map_spin;

	Algo2(void);
	Algo2(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void print_state(float temdp);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
	float eval_lat();
	float eval_lat_spin();
	float eval_site_spin(int site);
	float eval_site_chem(int site);
	float eval_spin_flip(int site, float old_spin);
	float eval_atom_flip(int site1, float old_spin1, int site2, float old_spin2);
};

#endif
