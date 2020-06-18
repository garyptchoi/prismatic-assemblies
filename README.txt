Codes and link patterns for the deterministic and stochastic control of connectivity and rigidity in prismatic assemblies 

Reference:
G. P. T. Choi, S. Chen, and L. Mahadevan, 
"Control of connectivity and rigidity in prismatic assemblies."
Preprint, arXiv:2005.05845, 2020.

Copyright (c) 2020, Gary Pui-Tung Choi, Siheng Chen, L. Mahadevan

==============================================
Deterministic control of prismatic assemblies:
The following MATLAB files are used for constructing a minimum rigidifying link pattern (MRP) or a minimum connecting link pattern (MCP) for a given prismatic assembly.

deterministic/MRP_rectangular_LxMxN.m: Construct a MRP for an LxMxN rectangular prismatic assembly
deterministic/MCP_rectangular_LxMxN.m: Construct a MCP for an LxMxN rectangular prismatic assembly
deterministic/MRP_rectangular_LxMxN.m: Construct a MRP for an LxMxN triangular prismatic assembly
deterministic/patterns: The MATLAB 3D plots of the patterns

---------------------------------------------
Stochastic control of prismatic assemblies:
stochastic/prismatic.cpp and stochastic/prismatic.hpp build a rectangular prismatic assembly with a prescribed number of random links and calculate the degree of freedom. 

prismatic Prismatic(L_cubic, L_cubic, L_cubic, n_cst_link);
Prismatic.gen_random_cst_list();
long rank = Prismatic.gen_DoF(Prismatic.constraint_all);

The first line builds an L_cubic*L_cubic*L_cubic rectangular prismatic assembly with n_cst_link links.
The second line assigns the n_cst_link links randomly among all possible locations.
The third line calculates the degree of freedom using the sparse matrix rank determination from SuiteSparse, which can be downloaded at http://faculty.cse.tamu.edu/davis/suitesparse.html.

The demo stochastic/dofvscst.cpp calculates the degree of freedom at 20 different link densities for a 20*20*20 assembly. At each link density, 100 random link patterns are generated.
