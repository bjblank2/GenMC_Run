#pragma once
#ifndef mc_h
#define mc_h
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
#include "Session.h"
using namespace std;

class MCrun
{
public:
	MCrun(void);
	MCrun(Session& _session, SimCell& _sim_cell);

	Session session;
	SimCell sim_cell;
	const double Kb = 0.00008617333262; // Boltzmann constant
	const double uB = .000057883818012; // Bhor magnaton 
	int site_species;

	void start(vector<Rule>& mc_rules);
	void runMetropolis3(vector<Rule>& mc_rules);
	void runMetropolis4(vector<Rule>& mc_rules);
	void runMetropolis5(vector<Rule>& mc_rules);
	void runMetropolis6(vector<Rule>& mc_rules);
	void runMetropolis7(vector<Rule>& mc_rules);
	void runMetropolis8(vector<Rule>& mc_rules);
	void runMetropolis9(vector<Rule>& mc_rules);
	void runMetropolis10(vector<Rule>& mc_rules);
	void outputSimCell(vector<Rule>& mc_rules);
	void writeSuperCell(string out_path, vector<int>& atom_species, vector<int>& atom_spins, SimCell& sim_cell);
	long double evalLattice(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalLattice3(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalLattice_tri(map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalLatticeSpin(map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalLatticeSpin3(map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalSiteEnergyAll(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalSiteEnergyAll_tri(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalSiteEnergySpin(int site, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalSiteEnergyAll3(int site, map<string, long double>& rule_map_spin, map<string, long double>& rule_map_chem, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double evalSiteEnergySpin3(int site, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	long double delSiteEnergySpin(int site, int old_spin, map<string, long double>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float evalLatticeFast(vector<int> spins, vector<int> species, vector<vector<vector<int>>> JcountChem, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatChem, vector<vector<float>> JMatSpin);
	float evalLatticeSpinFast(vector<int> spins, vector<int> species, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatSpin);
	float evalLattiveChemFast(vector<vector<vector<int>>> JcountChem, vector<vector<float>> JMatChem);
	float evalSiteEnergyAllFast(int spin, vector<int> spins, vector<vector<int>> JcountChem, vector<vector<int>> JcountSpin, vector<vector<float>> JMatChem, vector<vector<float>> JMatSpin);
	float evalSiteEnergySpinFast(int spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin);
	float evalSiteEnergyChemFast(vector<vector<int>> JcountChem, vector<vector<float>> JMatChem);
	float evalDelEnergySpinFast(int old_spin, int new_spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin);
	float evalPottsLatticeSpinFast(vector<int> spins, vector<vector<vector<int>>> JcountSpin, vector<vector<float>> JMatSpin);
	float evalPottsSiteEnergySpinFast(int spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin);
	float evalPottsDelEnergySpinFast(int old_spin, int new_spin, vector<int> spins, vector<vector<int>> JcountSpin, vector<vector<float>> JMatSpin);
	float calcMag2(int site, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float calcMag2Diff(int site, int old_spin, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float calc_struct(int site, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float calcWarrenCrowley(int site, float numb_Mn, float numb_In, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float calcSpecies(int site, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
	float stagMag(int site, int spin, SimCell& sim_cell);
	float stagMagDiff(int site, int spin, int old_spin, SimCell& sim_cell);
};
#endif