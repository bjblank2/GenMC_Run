#include "utils.h"

vector<string> split(string str, const string delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}

vector<string> split(string str)
{
	istringstream iss(str);
	vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
	return results;
}

int get_index(vector<int> vect, int elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<string> vect, string elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<float> vect, float elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}