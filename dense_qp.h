int run_gurobi(int nRows, // Number of constraints
               int nCols, // Number of optimisation variables
               int nInputs, // Number of inputs - not used
               double *Q,
               double *c,
               double *A,
               double *rhs,
               double *w,
               double *UOpt // MPC output
               );