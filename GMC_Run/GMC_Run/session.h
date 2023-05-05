#pragma once
#ifndef session_h
#define session_h
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "file_io.h"
#include "rule.h"
#include "sim_cell.h"
using namespace std;

class Session
{
public:
	bool use_poscar;
	bool use_states;
	int algo = 0;
	int numb_passes = 1;
	int eq_passes = 1;
	int shape[3] = { 0,0,0 };
	float start_temp = 1;
	float end_temp = 0;
	float temp_inc = 0.5;
	float sro_temp = 0;
	float sro_target;
	float mag_ext = 0;
	string structure_file;
	string rules_file;
	string sim_type;
	string phase_init;
	string spin_init;
	string species_init;
	vector<int> atom_numbs;
	vector<int> species_inds;
	vector<float> moments;
	vector<float> unique_dists;
	vector<string> species;
	vector<Rule> rule_list;
	vector<vector<float>> spin_states;
	vector<vector<vector<float>>> motif_list;


	Session(void);
	Session(Session& _session);
	Session(string input_file);

	void _copy(Session& _session);
	void add_spin_states(string input_file);
	void fill_rule_list();
	void find_unique_dists();
};
#endif
