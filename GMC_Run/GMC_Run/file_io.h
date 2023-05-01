#pragma once
#ifndef file_io_h
#define file_io_h
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

#endif