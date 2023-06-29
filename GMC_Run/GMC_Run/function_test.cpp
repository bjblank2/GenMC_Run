#include "function_test.h"

test::test(void) {}

test::test(Session& _session, SimCell& _sim_cell) {
    session = _session;
    sim_cell = _sim_cell;
}

void test::run(){
    //setup rng for random atom/spin method
    std::random_device dev2;
    std::mt19937 rng3(dev2());
    std::uniform_int_distribution<std::mt19937::result_type> rand_method(0, 100000);
    int x = 0;
    int y = 0;
    int z = 0;
    for (int i = 0; i < 1000000; i++) {
        if (rand_method(rng3) < 33333) {
            x++;
        }
        else if (rand_method(rng3) < 50000) {
            y++;
        }
        else {
            z++;
        }
    }
    cout << x << " " << y << " " << z;
}
