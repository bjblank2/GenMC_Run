#pragma once
#ifndef algo3_h
#define algo3_h

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

class Algo3 {
public:
	int outfile_count = 0;
	SimCell sim_cell;
	Session session;
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = .000057883818012; // Bhor magnaton
    vector<int> chem_list;
    vector<float> spin_list;
    vector<float> site_rule_count_list;
    vector<float> lat_rule_count_list;
    vector<vector<float>> pos_list;
    vector<vector<vector<vector<int>>>> spin_motif_groups;
    vector<vector<vector<vector<int>>>> chem_motif_groups;
    map <string, vector<float>> rule_map_chem;
    map <string, float> rule_map_spin;

	Algo3(void);
	Algo3(Session& _session, SimCell& _sim_cell);
	void run();
	//void eval_lat_sro();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void fill_SROMG(vector<vector<int>>& neigh_ind_list);
	void print_state(string contcar_name, int temp);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
	float init_SRO(vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list);
	float eval_lat();
	float eval_lat_spin();
	float eval_site_spin(int site);
	float eval_site_chem(int site);
	float eval_spin_flip(int site, float old_spin);
	float eval_atom_flip(int site);
    //vector<int> eval_site_sro(int site);
    vector<int> eval_lat_sro();
    vector<int> eval_site_sro(int site);
};

#endif
