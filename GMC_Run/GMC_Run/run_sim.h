#pragma once
#ifndef run_sim_h
#define run_sim_h

#include "session.h"
#include "sim_cell.h"
#include "algo1.h"
#include "algo2.h"
#include "algo3.h"
#include "function_test.h"

class Sim {
public:
	SimCell sim_cell;
	Session session;

	Sim(void);
	Sim(Session& _session, SimCell& _sim_cell);
	void start();
};

#endif
