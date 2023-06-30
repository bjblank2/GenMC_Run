#pragma once
#ifndef sro_h
#define sro_h
#include "utils.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;


class SRO{
private:
	int type; // SRO group type
public:
	vector<vector<float>> motif;
	int motif_ind;
	vector<int> deco;
	vector<float> dists;
	SRO(void);
	SRO(int _type, vector<int> _deco, vector<vector<float>> _motif, int motif_ind);
	int GetType();
	int GetLength();
	vector<int> GetDeco();
	vector<float> GetDists();
};

#endif