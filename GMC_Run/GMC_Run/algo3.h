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
#define ONE_SPIN 0
#define NO_SPIN 0
#define ONE_ATOM 0
#define METHOD_0 0
#define METHOD_1 1
#define METHOD_2 2
#define METHOD_3 3

class Algo3 {
public:
	int flip_count = 0;
	int flip_count2 = 0;
	int outfile_count = 0;
	int passes;
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = 0.00005788381806; // Bhor magnaton
	float spin_flip = 0.0;
	float e_flip = 0.0;
	float e_avg = 0.0;
	float spin_avg = 0.0;
	float init_enrg = 0.0;
	float init_spin = 0.0;
	float var_spin = 0.0;
	float var_e = 0.0;
	RunningStat rs_C;
	RunningStat rs_X;
	SimCell sim_cell;
	Session session;
	vector<int> chem_list;
	vector<float> count_avg;
	vector<float> init_sro;
	vector<float> spin_list;
	vector<float> site_rule_count_list;
	vector<float> lat_rule_count_list;
	vector<vector<float>> pos_list;
	vector<vector<vector<vector<int>>>> spin_motif_groups;
	vector<vector<vector<vector<int>>>> chem_motif_groups;
	map <string, vector<float>> rule_map_chem;
	map <string, float> rule_map_spin;
	// setup rng for random spin choice and acceptance probability
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	std::uniform_real_distribution<double> unif;
	std::uniform_int_distribution<int> rand_atom;
	std::uniform_int_distribution<int> rand_method;



	Algo3(void);
	Algo3(Session& _session, SimCell& _sim_cell);
	void run();
	void fill_CMG(vector<vector<int>>& neigh_ind_list);
	void fill_SMG(vector<vector<int>>& neigh_ind_list);
	void print_state(string contcar_name, int temp);
	void spin_move(int site, int pass, float temp, float new_spin, ofstream& Output_converge);
	void spec_move(int site, int rand_site, int pass, float temp, float new_spin, ofstream& Output_converge);
	void atom_move(int site, int rand_site, float new_spin1, float new_spin2, int pass, float temp, float new_spin, ofstream& Output_converge);
	bool bc_check(vector<float> check_vect, vector<float>& pos);
	float eval_lat();
	float eval_lat_spin();
	float eval_site_spin(int site);
	float eval_site_chem(int site);
	float eval_spin_flip(int site, float old_spin);
	float eval_atom_flip(int site);
	vector<int> eval_site_sro(int site);
};

#endif
