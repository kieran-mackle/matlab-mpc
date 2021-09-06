/*
  Solve an LP/QP/MILP/MIQP represented using dense matrices.  This
  routine assumes that A and Q are both stored in row-major order.
  It returns 1 if the optimization succeeds.  When successful,
  it returns the optimal objective value in 'objvalP', and the
  optimal solution vector in 'solution'.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gurobi_c.h>
#include "dense_qp.h"

static int
dense_optimize(GRBenv   *env,
               int      rows,
               int      cols,
               double   *c,        // linear part of objective function
               double   *Q,        // quadratic part of objective function
               double   *A,        // constraint matrix
               char     *sense,    // constraint senses
               double   *rhs,      // RHS vector
               double   *w,        // constraint weights
               double   *lb,       // variable lower bounds
               double   *ub,       // variable upper bounds
               char     *vtype,    // variable types (continuous, binary, etc.)
               double   *solution,
               double   *objvalP)
{
    GRBmodel    *model = NULL;
    int         i, j, optimstatus;
    int         error = 0;
    int         success = 0;

    // Create an empty model
    error = GRBnewmodel(env, &model, "dense", cols, c, lb, ub, vtype, NULL);
    if (error) goto QUIT;

    error = GRBaddconstrs(model, rows, 0, NULL, NULL, NULL, sense, rhs, NULL);
    if (error) goto QUIT;

    // Update constraints
    error = GRBupdatemodel(model);
    if (error) goto QUIT;

    // Populate A matrix
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (A[i*cols+j] != 0) {
                error = GRBchgcoeffs(model, 1, &i, &j, &A[i*cols+j]);
                if (error) goto QUIT;
            }
        }
    }

    // Populate Q matrix
    if (Q) {
        for (i = 0; i < cols; i++) {
            for (j = 0; j < cols; j++) {
                if (Q[i*cols+j] != 0) {
                    error = GRBaddqpterms(model, 1, &i, &j, &Q[i*cols+j]);
                    if (error) goto QUIT;
                }
            }
        }
    }

    // Update constraints
    error = GRBupdatemodel(model);
    if (error) goto QUIT;

    // Settings to handle numerical difficulties
    error = GRBsetintparam(GRBgetenv(model), "NumericFocus", 3);
    error = GRBsetintparam(GRBgetenv(model), "ScaleFlag", 2);

    // Add slack variables (currently handling all linear constraint senses)
    int vind[1];
    double vval[1];
    for (i = 0; i < rows; i++) {
        if (!isnan(w[i])) {
            vind[0] = i;
            vval[0] = -1.0;
            error = GRBaddvar(model, 1, vind, vval, w[i], 0.0, GRB_INFINITY,
                                GRB_CONTINUOUS, NULL);
            if (error) goto QUIT;
        }
    }

    // Optimize model
    error = GRBoptimize(model);
    if (error) goto QUIT;

    // Write model to 'dense.lp'
    error = GRBwrite(model, "dense_qp.lp");
    if (error) goto QUIT;

    // Capture solution info
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;

    if (optimstatus == GRB_OPTIMAL) {
        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objvalP);
        if (error) goto QUIT;
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, solution);
        if (error) goto QUIT;
        success = 1;
    }

QUIT:
    // Error Reporting
    if (error) {
        printf("ERROR %s\n", GRBgeterrormsg(env));
        exit(1);
    }

    // Free model
    GRBfreemodel(model);
    return success;
}

int
run_gurobi(int nRows, // Number of constraints
           int nCols, // Number of optimisation variables
           int nInputs, // Number of inputs - not used
           double *Q,
           double *c,
           double *A,
           double *rhs,
           double *w,
           double *UOpt // MPC output
           )
{
    GRBenv  *env    = NULL;
    int     error   = 0;
    int     solved;
    double  objval;

    // Initialise 'sense' array
    char sense[nRows];
    for (int i = 0; i < nRows; i++) {
        sense[i] = '<';
    }

    // Initialise lower bound array
    double lb[nRows];
    for (int i = 0; i < nRows; i++) {
        lb[i] = -GRB_INFINITY;
    }

    // Create environment
    error = GRBloadenv(&env, "dense_qp.log");
    if (error) goto QUIT;

    // // Suppress output if desired
    // error = GRBsetintparam(env, GRB_INT_OUTPUT_FLAG, 0);
    // if (error) goto QUIT;

    // Solve the model
    solved = dense_optimize(env, nRows, nCols, c, Q, A, sense, rhs, w, lb, NULL,
            NULL, UOpt, &objval);

    // // Print outputs
    // if (solved)
    // {
        // // Write MPC inputs
        // printf("Solved: u0=%.8f, u1=%.8f, u2=%.8f, u3=%.8f\n", UOpt[0], UOpt[1],
               // UOpt[2], UOpt[3]);
    // }

    GRBfreeenv(env);
    return 1;

QUIT:
    printf("ERROR %s\n", GRBgeterrormsg(env));
    GRBfreeenv(env);
    return 0;
}
