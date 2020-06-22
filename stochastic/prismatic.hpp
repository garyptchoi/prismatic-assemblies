// Stochastic control of prismatic assemblies
//
// Reference:
// G. P. T. Choi, S. Chen, L. Mahadevan, 
// ``Control of connectivity and rigidity in prismatic assemblies.''
// Preprint, arXiv:2005.05845, 2020.

#ifndef PRISMATIC_HPP
#define PRISMATIC_HPP

#define _USE_MATH_DEFINES
#define FOLD_CLOSE_PAIRS 0
#define FOLD_RANDOM 1
#define CURR_FOLD_OPTION 1

#include <math.h>
//#include <iostream>
#include <math.h>
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include "your-directory-to-SuiteSparse/SuiteSparse-5.7.1/include/SuiteSparseQR.hpp"
class prismatic
{
public:
  prismatic(long input_n_row_cubic, long input_n_col_cubic, long input_n_layer_cubic, long input_n_link);
  void add_A_entry(long r, long c, double x);
  void gen_All_Link_List();
  void add_entry_two(long e1, long e2);
  void add_entry_four(long e1, long e2, long e3, long e4);
  void gen_edge_Constraint_Triplets();
  void gen_Link_Constraint_Triplets(std::vector<long> constraint_list);
  void add_link(long e1, long e2);
  SuiteSparse_long gen_DoF(std::vector <long> constraint_list);
  void gen_random_cst_list();

  void find_largest_cluster();
  void DFS(long start_point, long cluster_label, bool visited[], long * cluster_list, long * rigid_list);
  void print_cluster(long * rigid_list, long* cluster_list);

  long rank;
  std::vector <long> constraint_all;
  std::vector <std::vector<long>> link_list;

  long n_free;
  long n_nodes_in_largest_cluster;
  long n_nodes_in_largest_cluster_square;
  long n_cluster;
  long n_cst_link_all;

private:
  cholmod_triplet * A;

  double x_2el = sqrt(2)/2.;
  double x_4el = 1./2.;

  long n_row_cubic, n_col_cubic, n_layer_cubic;
  long n_cubic, n_node, n_cubic_per_layer;
  long n_cst_edges, n_cst_link;
  long n_link;

  long row_num; //Tracking how many lines have been added to the triplet A
  long shuffle_done;



};

#endif
