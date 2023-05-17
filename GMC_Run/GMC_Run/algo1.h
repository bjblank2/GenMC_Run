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

class Algo1{
public:
	SimCell sim_cell;
	Session session;
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = .000057883818012; // Bhor magnaton 
	vector<int> chem_list;
	vector<float> spin_list;
	vector<vector<float>> pos_list;
	vector<vector<vector<int>>> spin_motif_groups;
	vector<vector<vector<int>>> chem_motif_groups;
	map <string, float> rule_map_chem;
	map <string, float> rule_map_spin;

	Algo1(void);
	Algo1(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	float init_SRO(vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list);
	float calc_struct(int site, vector<vector<int>>& neigh_ind_list, vector<vector<float>>& neigh_dist_list);
	float eval_lat();
	float eval_lat_spin();
	float eval_site_spin(int site);
	float eval_site_chem(int site);
	void print_state();
};

#endif