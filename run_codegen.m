cfg = coder.config('mex');
cfg.CustomInclude = fullfile("/opt/gurobi912/linux64/include/");
cfg.CustomLibrary = fullfile("/opt/gurobi912/linux64/lib/libgurobi91.so");
codegen -config cfg solve_qp -args {0,0,0,H,G,Omega,omega,constraint_weights, Uopt_mat} dense_qp.c -v