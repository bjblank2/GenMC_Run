#ifndef function_test_h
#define function_test_h

#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <cmath>
#include <random>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "rule.h"
#include "sim_cell.h"
#include "running_stat.h"
#include "session.h"
#include "utils.h"
using namespace std;

class test{
public:
    SimCell sim_cell;
    Session session;
    test(void);
    test(Session& _session, SimCell& _sim_cell);
    void run();
};

#endif /* function_test_h */
