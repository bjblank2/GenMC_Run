#include "function_test.h"

test::test(void) {}

test::test(Session& _session, SimCell& _sim_cell) {
    session = _session;
    sim_cell = _sim_cell;
}

void test::run(){
    cout << "1";
}
