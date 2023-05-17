#pragma once
#ifndef utils_h
#define utils_h
#include <stdio.h>
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
using namespace std;

vector<string> split(string str, const string delim);
vector<string> split(string str);
vector<float> pos_shift(vector<float>& vect1, vector<float>& vect2);
int get_index(vector<int> vect, int elem);
int get_index(vector<float> vect, float elem);
int get_index(vector<string> vect, string elem);
int vect_max(vector<int> vect);
float vect_max(vector<float> vect);
vector<int> vect_permut(vector<int>& vect);
vector<int> vect_permut(vector<float>& vect);
void sort_vect(vector<int>& vect, vector<int>& perm);
void sort_vect(vector<float>& vect, vector<int>& perm);
void sort_vect(vector<vector<int>>& vect, vector<int>& perm);
void sort_vect(vector<vector<float>>& vect, vector<int>& perm);



#endif