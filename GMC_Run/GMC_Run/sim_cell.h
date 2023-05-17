#pragma once
#ifndef sim_cell_h
#define sim_cell_h
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "utils.h"
#include "session.h"
using namespace std;


class SimCell {
public:
	int sup_cell[3];
	int numb_atoms;
	int numb_cells[3];
	float cell_dim[3];
	float cutoff;
	float unit_LC[3];
	string sim_type;
	string phase_init;
	vector<int> species_types;
	vector<int> species_numbs;
	vector<int> poscar_comp;
	class Atom {
	private:
		int species;
		float spin;
		int phase;
	public:
		int index;
		float pos[3];
		vector<int> neighbors;
		vector<float> neighbor_dists;
		vector<int> allowed_species; //list of species that can occupy this position on the lattice (used as a way to discribe any sublattice) 
		Atom(void);
		Atom(int _index, int _species, float _spin, int _phase, vector<float> _pos, vector<int> _allowed_species);
		void setSpin(float _spin);
		void setSpecies(int _species);
		void setPhase(int _phase);
		float getSpin(void);
		int getSpecies(void);
		int getPhase(void);
		int getNeighborSpin(int _neighbor, SimCell& sim_cell);
		int getNeighborSpecies(int _neighbor, SimCell& sim_cell);
		int getNeighborPhase(int _neighbor, SimCell& sim_cell);
		int getNeighborIndex(int _neighbor, SimCell& sim_cell);
		int getNumbNeighbors(int _site, SimCell& sim_cell);
		float getNeighborDist(int _neighbor, SimCell& sim_cell);
	};
	vector<Atom> atom_list;
	vector<Atom> unit_cell;

	SimCell(void);
	//SimCell(string POSCAR_file, int _sup_cell[3], vector<int>& _species_numbs, float _cutoff, string _sim_type, string phase_init, string spin_init, string species_init);
	SimCell(SimCell& sc_copy);
	void _copy(SimCell& sc_copy);
	void initSimCell(string POSCAR_file, Session& session);
	void fillUnitCell(string POSCAR_file, Session& session);
	void make_supercell(Session& sess);
	void setNeighborDists(vector<float>& dist_list);
	float findAtomDists(int atom1, int atom2);
};
#endif
