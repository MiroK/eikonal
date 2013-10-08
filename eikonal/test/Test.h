/*

Test("name")

run(N_min, N_max)

{
  for N in N_min, N_max  //doubling
  
  mesh = ...
  V = ...

  Function u(V)
  Function u_exact(V)
  std::set<la_index> fixed_dofs;

  Problem problem(!Seeder!)

  problem.init(...)
  problem.exact_solution(...)

  Solver !solver!(v)
  iters = solver.solve(u, fixed_dofs);

  // compute L1, L2 norms of (u-u_exact)
  // write to file mesh.hmin(), L1, L2, iters
}















*/
