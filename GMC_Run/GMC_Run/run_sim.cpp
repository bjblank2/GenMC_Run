#include "run_sim.h"

Sim::Sim(void) {}

Sim::Sim(Session& _session, SimCell& _sim_cell) {
	sim_cell = _sim_cell;
	session = _session;
}

void Sim::start() {
	cout << "using algo " << session.algo << "\n";
	if (session.algo == 0) { return; }
	else if (session.algo == 1) { Algo1 algo(session, sim_cell); algo.run(); }
	//else if (session.algo == -2) { debug_run(); }
}