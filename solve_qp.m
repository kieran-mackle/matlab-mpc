function [output, Uopt] = solve_qp(nRows, nCols, nInputs, H, G, Omega, omega, constraint_weights, Uopt_mat)
%#codegen

refQ        = reshape(H.',1,[])';        % Quadratic objective matrix   % Should it be 2H?
refc        = reshape(-G.',1,[])';
refA        = reshape(Omega.',1,[])';
refRHS      = reshape(omega.',1,[])';
refw        = constraint_weights;
refUopt     = Uopt_mat;

% pre-initialize output
output = 0;
coder.cinclude('stdlib.h');
coder.cinclude('stdio.h');
coder.cinclude('math.h');
coder.cinclude('gurobi_c.h');
coder.cinclude('dense_qp.h');

output = coder.ceval('run_gurobi', nRows, nCols, nInputs, ...
                       coder.ref(refQ), coder.ref(refc), ...
                       coder.ref(refA), coder.ref(refRHS), ...
                       coder.ref(refw), coder.wref(refUopt));

Uopt = refUopt;
end