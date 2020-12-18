// Stochastic control of prismatic assemblies
//
// Reference:
// G. P. T. Choi, S. Chen, L. Mahadevan, 
// ``Control of connectivity and rigidity in prismatic assemblies.''
// Proceedings of the Royal Society A, 476(2244), 20200485, 2020.

#include <iostream>
#include <math.h>
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <fstream>
#include "prismatic.hpp"
#include <omp.h>
#include <cstdio>
#include "your-directory-to-SuiteSparse/SuiteSparse-5.7.1/include/SuiteSparseQR.hpp"



int main(int argc, char* argv[]) {
  srand(time(NULL));
  double r = ((double) rand() / (RAND_MAX));
  int L_init = 20;
  int L_max = 20;
  int cst_bin = 20;  // number of sampling points
  int repeats = 100; // number of random initializations Default: 200
  for (int L_cubic = L_init; L_cubic <= L_max; L_cubic+=5){
    std::ofstream myfile, myfile2, myfile3;
    myfile.open("results/rank_L"  + std::to_string(L_cubic)+".txt", std::ofstream::out | std::ofstream::app);
    myfile2.open("results/link_L" + std::to_string(L_cubic)+".txt", std::ofstream::out | std::ofstream::app);
    // myfile3.open("results/alllinks_L" + std::to_string(L_cubic)+".txt", std::ofstream::out | std::ofstream::app);
    long n_cst_link_max = (L_cubic-1)*L_cubic*L_cubic*4*3 + (L_cubic-1)*(L_cubic-1)*4*L_cubic*3 + 4*(L_cubic-1)*(L_cubic-1)*(L_cubic-1);
    long step_size = round(n_cst_link_max/cst_bin);
    #pragma omp parallel for num_threads(25) collapse(2) schedule(dynamic)
    for (long n_cst_link = 0; n_cst_link <= n_cst_link_max; n_cst_link+=step_size) {
      for (int tt = 0; tt < repeats; tt++) {
          prismatic Prismatic(L_cubic, L_cubic, L_cubic, n_cst_link);
          Prismatic.gen_random_cst_list();
          long rank = Prismatic.gen_DoF(Prismatic.constraint_all);
    #pragma omp critical
          {
          myfile << n_cst_link << " " << rank << std::endl;
          for (auto c:Prismatic.constraint_all)
            myfile2 << c << " ";
          myfile2 << std::endl;
          // for (auto c:Prismatic.link_list)
          //   myfile3 << c[0] << " " << c[1] << std::endl;
          }
        }
      }
      myfile.close();
      myfile2.close();
      // myfile3.close();
  }
  return 0;
}
