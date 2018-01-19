#include <iostream>
#include "force2d.h"
#include "particle2d.h"
#include "integrate2d.h"

using namespace std;

int main() {
  int N = 1000;
  double Lx = 100;
  double Ly = Lx;
  double h = 1e-4;
  double Pe = 10;

  BP_with_ori_2 *par = new BP_with_ori_2[N];
  Ran *myran = new Ran(1);
  
  create_rand_2(par, N, 1, Lx, Ly, myran);
  
  F_WCA_2 fWCA(Lx, Ly, 1);
  Int_EM_2 euler(h, Lx, Ly);

  for (int i = 0; i <= 1000; i++) {
    cal_pair_force_simply(par, N, fWCA);
    euler.int_SP_all(par, N, myran, Pe);
  }
}