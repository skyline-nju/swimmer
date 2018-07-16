#include <iostream>
#include "benchmark.h"

using namespace std;

int main(int argc, char* argv[]) {
  lbm::D2Q9::benchmark_taylor_green(200, 200);
}