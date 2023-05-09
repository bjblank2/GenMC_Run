#pragma once
#ifndef rule_h
#define rule_h
#include "utils.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;


class Rule {
private:
	long double enrg_cont;
	int type; // spin or chem
	int phase;
public:
	vector<vector<float>> motif;
	vector<int> deco;
	vector<float> dists;
	Rule(void);
	Rule(long double _enrg_cont, int _type, int _phase, vector<int> _deco, vector<vector<float>> _motif);
	long double GetEnrgCont();
	int GetType();
	int GetPhase();
	int GetLength();
	bool IsRuleSpin(vector<int> test_species, vector<float> test_dists);
	bool IsRuleChem(vector<int> test_species, vector<float> test_dists);
	vector<int> GetDeco();
	vector<float> GetDists();
};

#endif