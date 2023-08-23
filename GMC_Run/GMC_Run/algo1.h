#pragma once
#ifndef algo1_h
#define algo1_h

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
#define DONE -1

class Algo1 {
public:
    const double Kb = 0.00008617333262; // Boltzmann constant
    const double uB = 0.00005788381806; // Bohr magnaton
	int outfile_count = 0;
	SimCell sim_cell;
	Session session;
	vector<int> chem_list;
	vector<float> spin_list;
	vector<vector<float>> pos_list;
	vector<vector<vector<vector<int>>>> spin_motif_groups;
	vector<vector<vector<vector<int>>>> chem_motif_groups;
	map <string, double> rule_map_chem;
	map <string, double> rule_map_spin;

	Algo1(void);
	Algo1(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
    void print_state(string contcar_name, int temp);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
    double eval_lat();
    double eval_lat_spin();
    double eval_site_spin(int site);
    double eval_site_chem(int site);
    double eval_spin_flip(int site, float old_spin);
    double eval_atom_flip(int site);
};

#endif
