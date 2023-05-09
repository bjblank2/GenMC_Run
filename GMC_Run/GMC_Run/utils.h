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
int get_index(vector<int> vect, int elem);
int get_index(vector<float> vect, float elem);
int get_index(vector<string> vect, string elem);

#endif