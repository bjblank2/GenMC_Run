#include "utils.h"

vector<string> split(string str, const string delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
    str.erase(remove(str.begin(), str.end(), '\r'), str.end());
    str.erase(remove(str.begin(), str.end(), '\n'), str.end());

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

vector<float> pos_shift(vector<float>& vect1, vector<float>& vect2){
    vector<float> temp;
    for (int i = 0; i < 3; i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

vector<float> pos_shift(vector<float>& vect1, float vect2[3]) {
    vector<float> temp;
    for (int i = 0; i < 3; i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

int get_index(vector<int>& vect, int elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<string>& vect, string elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<float>& vect, float elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int vect_max(vector<int>& vect) {
    int max = vect[0];
    for (int i : vect) { if (i > max) { max = i; } }
    return max;
}

float vect_max(vector<float>& vect) {
    float max = vect[0];
    for (float i : vect) { if (i > max) { max = i; } }
    return max;
}

vector<int> vect_permut(vector<int>& vect) {
    vector<int> vect_temp;
    vect_temp.assign(vect.begin(), vect.end());
    vector<pair<int, int>> indices;
    for (int i = 0; i < vect_temp.size(); ++i) {
        indices.push_back(make_pair(vect_temp[i], i));
    }
    sort(indices.begin(), indices.end());

    vector<int> permut(vect.size());
    for (int i = 0; i < indices.size(); ++i) {
        permut[i] = indices[i].second;
    }
    return permut;
}

vector<int> vect_permut(vector<float>& vect) {
    vector<int> vect_temp;
    vect_temp.assign(vect.begin(), vect.end());
    vector<pair<int, int>> indices;
    for (int i = 0; i < vect_temp.size(); ++i) {
        indices.push_back(make_pair(vect_temp[i], i));
    }
    sort(indices.begin(), indices.end());

    vector<int> permut(vect.size());
    for (int i = 0; i < indices.size(); ++i) {
        permut[i] = indices[i].second;
    }
    return permut;
}

void sort_vect(vector<int>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<int> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<float>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<float> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<vector<float>>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<vector<float>> temp;
    temp = vect;// .insert(temp.end(), vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<vector<int>>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<vector<int>> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

vector<float> pos_transform(vector<float>& pos, vector<vector<float>>& trans) {
    //cartesian coordinate vector x*a+y*b+z*c
    vector<float> vect(pos.size());
    vect[0] += pos[0] * trans[0][0] + pos[1] * trans[1][0] + pos[2] * trans[2][0];
    vect[1] += pos[0] * trans[0][1] + pos[1] * trans[1][1] + pos[2] * trans[2][1];
    vect[2] += pos[0] * trans[0][2] + pos[1] * trans[1][2] + pos[2] * trans[2][2];
    return vect;
}
