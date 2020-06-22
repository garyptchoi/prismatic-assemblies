// Stochastic control of prismatic assemblies
//
// Reference:
// G. P. T. Choi, S. Chen, L. Mahadevan, 
// ``Control of connectivity and rigidity in prismatic assemblies.''
// Preprint, arXiv:2005.05845, 2020.

#include "prismatic.hpp"
#include <iostream>
#include <math.h>
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <cstdio>
#include "your-directory-to-SuiteSparse/SuiteSparse-5.7.1/include/SuiteSparseQR.hpp"

prismatic::prismatic(long input_n_row_cubic, long input_n_col_cubic, long input_n_layer_cubic, long input_n_link){
  n_row_cubic = input_n_row_cubic;     // in the y direction
  n_col_cubic = input_n_col_cubic;     // in the x direction
  n_layer_cubic = input_n_layer_cubic; // in the z direction
  n_link = input_n_link;

  n_cubic = n_row_cubic * n_col_cubic * n_layer_cubic;
  n_cubic_per_layer = n_row_cubic * n_col_cubic;
  n_node = 8 * n_cubic;

  n_cst_edges = n_cubic * 18;
  n_cst_link = n_link * 3;

  row_num = 0;
  shuffle_done = 0;

  gen_All_Link_List();
  n_cst_link_all = link_list.size();

  gen_random_cst_list();


}
void prismatic::add_A_entry(long r, long c, double x)
{
    ((long*)A->i)[A->nnz] = r;
    ((long*)A->j)[A->nnz] = c;
    ((double*)A->x)[A->nnz] = x;
    (A->nnz)++;
}

void prismatic::add_entry_two(long e1, long e2) {

  add_A_entry(row_num, e1,   -x_2el);
  add_A_entry(row_num, e2,   +x_2el);
  row_num++;
}

void prismatic::add_entry_four(long e1, long e2, long e3, long e4) {

  add_A_entry(row_num, e1,   -x_4el);
  add_A_entry(row_num, e2,   +x_4el);
  add_A_entry(row_num, e3,   -x_4el);
  add_A_entry(row_num, e4,   +x_4el);
  row_num++;
}

//For a given cubic, add all the 18 constraints
// Order of the 8 vertices
//   8 - - 7
//  /     /|
// 5 - - 6 |
// | 4 - | 3
// |/    |/
// 1 - - 2

// i along the x, j along y.
// i from 0 to n_col_cubic (number of columns); j from 0 to n_row_cubic
// Calculation: j*n_col_node+i is the current
void prismatic::gen_edge_Constraint_Triplets()
{

  for (long i=0;i<n_col_cubic;i++){
    for (long j=0;j<n_row_cubic;j++){
      for (long k=0; k<n_layer_cubic; k++) {
        long f = k*n_cubic_per_layer + j * n_col_cubic + i; // First node at the bottom left corner;
        //std::cout << "f" << f << std::endl;

        //
        add_entry_two(24*f+0, 24*f+3);   //1-2 x
        add_entry_two(24*f+4, 24*f+7);   //2-3 y
        add_entry_two(24*f+6, 24*f+9);   //3-4 x
        add_entry_two(24*f+1, 24*f+10);  //4-1 y
        add_entry_two(24*f+2, 24*f+14);  //1-5 z
        add_entry_two(24*f+5, 24*f+17);  //2-6 z
        add_entry_two(24*f+8, 24*f+20);  //3-7 z
        add_entry_two(24*f+11, 24*f+23); //4-8 z
        add_entry_two(24*f+12, 24*f+15); //5-6 x
        add_entry_two(24*f+16, 24*f+19); //6-7 y
        add_entry_two(24*f+18, 24*f+21); //7-8 x
        add_entry_two(24*f+13, 24*f+22); //5-8 y

        add_entry_four(24*f+0,  24*f+6,  24*f+1,  24*f+7);  //1-3 xy
        add_entry_four(24*f+21, 24*f+15, 24*f+16, 24*f+22);  //6-8 xy
        add_entry_four(24*f+0,  24*f+15, 24*f+2,  24*f+17);  //1-6 xz
        add_entry_four(24*f+4,  24*f+19, 24*f+5,  24*f+20);  //2-7 yz
        add_entry_four(24*f+21,  24*f+6, 24*f+8,  24*f+23);  //3-8 xz
        add_entry_four(24*f+13, 24*f+10, 24*f+11, 24*f+14);  //4-5 yz

        //Alternatively
        // add_entry_four(24*f+0,  24*f+6,  24*f+1,  24*f+7);  //1-3 xy
        // add_entry_four(24*f+12, 24*f+18, 24*f+13, 24*f+19);  //1-6 xz
        // add_entry_four(24*f+0,  24*f+15, 24*f+2,  24*f+17);  //2-7 yz
        // add_entry_four(24*f+9,  24*f+18, 24*f+11, 24*f+20);  //3-8 xz
        // add_entry_four(24*f+1,  24*f+22, 24*f+2,  24*f+23);  //4-5 yz
        // add_entry_four(24*f+4,  24*f+19, 24*f+5,  24*f+20);  //6-8 xy

      }
    }
  }

  //return row_num;
}
void prismatic::add_link(long e1, long e2) {
  std::vector<long> pr={e1, e2};
  link_list.push_back(pr);
}
// Order of the cubes
//    [6] - [7]
//   /     / |
// [4] - [5] |
//  |     |  |
//  | [2] - [3]
//  |/    |/
// [0] - [1]
void prismatic::gen_All_Link_List() {
  // std::vector<std::vector<long>> link_list;
  // std::vector<long> pr = {0,0};
  for (int i=0; i<n_col_cubic; i++) {
    for (int j=0; j<n_row_cubic; j++) {
      for (int k=0; k<n_layer_cubic; k++) {
        long cub_i = k*n_cubic_per_layer + j*n_col_cubic + i;
        if (i < n_col_cubic-1) {
          // cub_i to cub_i+1 (to the right)
          add_link(8*cub_i+2, 8*cub_i+8+1);
          add_link(8*cub_i+3, 8*cub_i+8+4);
          add_link(8*cub_i+6, 8*cub_i+8+5);
          add_link(8*cub_i+7, 8*cub_i+8+8);
        }
        if (j < n_row_cubic-1) {
          // cub_i to cub_i+n_col_cubic (to the back)
          add_link(8*cub_i+3, 8*cub_i+8*n_col_cubic+2);
          add_link(8*cub_i+4, 8*cub_i+8*n_col_cubic+1);
          add_link(8*cub_i+7, 8*cub_i+8*n_col_cubic+6);
          add_link(8*cub_i+8, 8*cub_i+8*n_col_cubic+5);
        }
        if (k < n_layer_cubic-1) {
          // cub_i to cub_i+n_cubic_per_layer
          add_link(8*cub_i+5, 8*cub_i+8*n_cubic_per_layer+1);
          add_link(8*cub_i+6, 8*cub_i+8*n_cubic_per_layer+2);
          add_link(8*cub_i+7, 8*cub_i+8*n_cubic_per_layer+3);
          add_link(8*cub_i+8, 8*cub_i+8*n_cubic_per_layer+4);
        }
        if (i < n_col_cubic-1 && j < n_row_cubic-1) {
          // cub_i to cub_i+n_col_cubic+1 and cub_i+1 to cub_i+input_n_col_cubic
          add_link(8*cub_i+3, 8*cub_i+8*n_col_cubic+8+1);
          add_link(8*cub_i+7, 8*cub_i+8*n_col_cubic+8+5);
          add_link(8*cub_i+8+4, 8*cub_i+8*n_col_cubic+2);
          add_link(8*cub_i+8+8, 8*cub_i+8*n_col_cubic+6);
        }
        if (i < n_col_cubic-1 && k < n_layer_cubic-1) {
          // cub_i to cub_i + n_cubic_per_layer + 1 and cub_i+n_cubic_per_layer to cub_i+1
          add_link(8*cub_i+6, 8*cub_i+8*n_cubic_per_layer+8+1);
          add_link(8*cub_i+7, 8*cub_i+8*n_cubic_per_layer+8+4);
          add_link(8*cub_i+8+5, 8*cub_i+8*n_cubic_per_layer+2);
          add_link(8*cub_i+8+8, 8*cub_i+8*n_cubic_per_layer+3);
        }
        if (j < n_row_cubic-1 && k < n_layer_cubic-1) {
          // cub_i to cub_i + n_cubic_per_layer + n_col_cubic and cub_i + n_col_cubic to cub_i + n_cubic_per_layer
          add_link(8*cub_i+7, 8*cub_i+8*n_cubic_per_layer+8*n_col_cubic+2);
          add_link(8*cub_i+8, 8*cub_i+8*n_cubic_per_layer+8*n_col_cubic+1);
          add_link(8*cub_i+8*n_col_cubic+5, 8*cub_i+8*n_cubic_per_layer+4);
          add_link(8*cub_i+8*n_col_cubic+6, 8*cub_i+8*n_cubic_per_layer+3);
        }

        // Body diagonal (only goes up)
        if (k < n_layer_cubic-1) {
          if (i > 0) {
            if (j > 0)
              add_link(8*cub_i+5, 8*cub_i+8*n_cubic_per_layer-8*n_col_cubic-8+3);
            if (j < n_row_cubic-1)
              add_link(8*cub_i+8, 8*cub_i+8*n_cubic_per_layer+8*n_col_cubic-8+2);
          }
          if (i < n_col_cubic-1) {
            if (j > 0)
              add_link(8*cub_i+6, 8*cub_i+8*n_cubic_per_layer-8*n_col_cubic+8+4);
            if (j < n_row_cubic-1)
              add_link(8*cub_i+7, 8*cub_i+8*n_cubic_per_layer+8*n_col_cubic+8+1);
          }
        }
      }
    }
  }
}
/*
points:
4  3
1  2
*/
void prismatic::gen_Link_Constraint_Triplets(std::vector<long> constraint_list) {

  for (long t=0; t<n_link; t++){
    long node_i = link_list[constraint_list[t]][0]; // jth row (j is along y)
    long node_j = link_list[constraint_list[t]][1];
    //std::cout << row_num << "     " << node_i << "      " << node_j << std::endl;
    add_A_entry(row_num, 3*node_i-3, -x_2el);
    add_A_entry(row_num, 3*node_j-3, +x_2el);
    row_num++;
    add_A_entry(row_num, 3*node_i-2, -x_2el);
    add_A_entry(row_num, 3*node_j-2, +x_2el);
    row_num++;
    add_A_entry(row_num, 3*node_i-1, -x_2el);
    add_A_entry(row_num, 3*node_j-1, +x_2el);
    row_num++;
  }
}


// Generate DoF of a given prismatic assembly
SuiteSparse_long prismatic::gen_DoF(std::vector<long> constraint_list)
{
  //std::cout << "row_num " << row_num << " " << n_cst_link_all << " " << n_cst_edges+n_cst_link << std::endl;
  cholmod_common Common, * com;
  com = &Common;
  cholmod_l_start(com);
  //std::cout << "n_node" << n_node << std::endl;

  A = cholmod_l_allocate_triplet(n_cst_edges+n_cst_link, 3*n_node, n_cst_edges * 4 + n_cst_link * 4, 0, CHOLMOD_REAL, com);
  //Generate all edge constraints
  gen_edge_Constraint_Triplets();
  // Generate all coplanar constraints
  gen_Link_Constraint_Triplets(constraint_list);


  // Calculate Rank
  SuiteSparse_long econ = 0;//min(rgd_Matrix->nrow, rgd_Matrix->ncol);
  //cholmod_l_print_triplet(A, "A", com);
  cholmod_sparse * rgd_Matrix = cholmod_l_triplet_to_sparse(A, n_cst_edges * 4 + n_cst_link * 4, com);
  SuiteSparse_long rank = SuiteSparseQR <double> (SPQR_ORDERING_GIVEN, SPQR_DEFAULT_TOL, econ, rgd_Matrix, NULL, NULL, com);

  // Cleaning
  cholmod_l_free_triplet(&A, com);
  cholmod_l_free_sparse(&rgd_Matrix, com);
  cholmod_l_finish(com);
  SuiteSparse_long dof = 3 * n_node - rank - 6;
  return dof;
}


// Randomly generate n_cst_coplanar coplanar constraints.
// Result is updated in constraint_all and constraint_all_extended
// constraint_all_extended will be updated by find_Rigid_Quads() later
void prismatic::gen_random_cst_list() {
  std::vector<long> constraint_list(n_cst_link_all);
  std::iota (std::begin(constraint_list), std::end(constraint_list), 0);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(constraint_list.begin(), constraint_list.end(), gen);
  constraint_list.resize(n_link);

  constraint_all = constraint_list;
  //Make sure both are sorted initially
  std::sort (constraint_all.begin(),constraint_all.end());
  shuffle_done = 1;
}
