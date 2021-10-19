function output = mpc_updateV2(mpc_input)
% Calculates a single MPC update.
% This is a modified objective function version

% ----------------------------------------------------------------------- %
% Deconstruct input
% ----------------------------------------------------------------------- %
cost            = mpc_input.cost;
constraints     = mpc_input.constraints;
if isfield(mpc_input, 'penalty_method') 
    penalty_method  = mpc_input.penalty_method; 
end
if isfield(mpc_input, 'penalty_weight') 
    penalty_weight  = mpc_input.penalty_weight; 
end

control_model   = mpc_input.control_model;

Ts              = mpc_input.params.timestep;      % Sampling interval
Hp              = mpc_input.params.horizon;       % Prediction horizon

initial         = mpc_input.initial;
t               = mpc_input.t;

% ----------------------------------------------------------------------- %
% Define LTI model
% ----------------------------------------------------------------------- %
if strcmpi(control_model, 'statespace')
    plant           = mpc_input.statespace;
    f0d             = zeros(size(1));
else
    control_dynamics    = control_model.dynamics;
    control_output      = control_model.output;

    input.phase.time    = t;
    input.phase.state   = initial.state;
    input.phase.control = initial.control;
    input.auxdata       = control_model.auxdata;

    % Call myJac.m to evaluate the Jacobians at the input time
    A = myJac(control_dynamics, input, 'state');
    B = myJac(control_dynamics, input, 'control');

    % Construct C and D
    C = numerical_jac(control_output, input.phase.state);
    D = zeros(size(C,1), size(B,2));

    plant       = ss(A, B, C, D);
    n           = size(A, 1);     % Number of states
    
    f0          = control_dynamics(input).dynamics';
    augm        = expm([A, eye(size(A));zeros(size(A)), zeros(size(A))] * Ts);
    gamma_coefficient = augm(1:n, n+1:end);
    f0d         = gamma_coefficient * f0;
end

% Define plant and convert to discrete-time domain
% ------------------------------------------------
plant       = c2d(plant, Ts);
Ad          = plant.A;
Bd          = plant.B;
C           = plant.C;
r           = mpc_input.reference;

% Define state space parameters
n           = size(Ad, 1);     % Number of states
m           = size(Bd, 2);     % Number of control inputs
p           = size(C, 1);     % Number of outputs


% Define current state and transform system
% -----------------------------------------
x0          = input.phase.state;       % x(k)
u0          = input.phase.control;     % u(k-1)
xbar_k      = zeros(size(x0))';
ubar_km1    = zeros(size(u0))';

% Calculate current operating point
% -----------------------------------------
g0          = control_output(x0)';

% Define cost weighting matrices
% -----------------------------------------
Q           = C' * cost.output * C; % Convert output costs to  state costs
R           = cost.control;

% Define constraint matrices
% -----------------------------------------
E           = build_constraints(constraints.hard.rate);
F           = build_constraints(constraints.hard.input - u0');
G_const     = build_constraints(constraints.hard.output);

if strcmpi(mpc_input.solver, 'quadprog') && ...
   strcmpi(mpc_input.constraints.type, 'hard') ~= 1
    E_soft      = build_constraints(constraints.soft.rate);
    F_soft      = build_constraints(constraints.soft.input - u0');
    G_soft      = build_constraints(constraints.soft.output);
end

% Format constraint and weighting matrices
% -----------------------------------------
Q_cell      = repmat({Q}, 1, Hp);           % Repeats Q Hp+1 times
Q_fancy     = blkdiag(Q_cell{:});           % Diagonal matrix of Qc
R_cell      = repmat({R}, 1, Hp);
R_fancy     = blkdiag(R_cell{:});

a           = size(E,1);
E_nom       = E(1:a, 1:m);
w_nom       = E(:, end);
EE          = zeros(Hp*a, Hp*m+1);
E_weights   = zeros(Hp*a, 1);

b           = size(F,1);
F_nom       = F(1:b, 1:m);
f_nom       = F(:, end);
FF          = zeros(Hp*b, Hp*m+1);
F_fancy     = zeros(Hp*b, m*Hp);
F_weights   = zeros(Hp*b, 1);

c           = size(G_const,1);
g_nom       = G_const(:, end);
G_nom       = G_const(1:c, 1:p);
GG          = zeros(Hp*c, Hp*p+1);
G_weights   = zeros(Hp*c, 1);


% Soft constraints
if strcmpi(mpc_input.solver, 'quadprog') && ...
   strcmpi(mpc_input.constraints.type, 'hard') ~= 1
EE_soft     = zeros(Hp*a, Hp*m+1);
E_soft_nom  = E_soft(1:a, 1:m);
w_soft_nom  = E_soft(:, end);

FF_soft     = zeros(Hp*b, Hp*m+1);
F_fancy_soft = zeros(Hp*b, m*Hp);
F_soft_nom  = F_soft(1:b, 1:m);
f_soft_nom  = F_soft(:, end);

GG_soft     = zeros(Hp*c, Hp*p+1);
G_soft_nom  = G_soft(1:c, 1:p);
g_soft_nom  = G_soft(:, end);
end

% Build constraint matrices across horizon
% -----------------------------------------
for i = 1:Hp
    % Hard constraints
    EE(a*(i-1)+1:a*i, m*(i-1)+1:m*i)    = E_nom;
    EE(a*(i-1)+1:a*i, Hp*m+1)           = w_nom;

    FF(b*(i-1)+1:b*i, m*(i-1)+1:m*i)    = F_nom;
    FF(b*(i-1)+1:b*i, Hp*m+1)           = f_nom;

    GG(c*(i-1)+1:c*i, p*(i-1)+1:p*i)    = G_nom;
    GG(c*(i-1)+1:c*i, Hp*p+1)           = g_nom;

    % Soft constraints
    if strcmpi(mpc_input.solver, 'quadprog') && ...
        strcmpi(mpc_input.constraints.type, 'hard') ~= 1
        EE_soft(a*(i-1)+1:a*i, m*(i-1)+1:m*i)    = E_soft_nom;
        EE_soft(a*(i-1)+1:a*i, Hp*m+1)           = w_soft_nom;

        FF_soft(b*(i-1)+1:b*i, m*(i-1)+1:m*i)    = F_soft_nom;
        FF_soft(b*(i-1)+1:b*i, Hp*m+1)           = f_soft_nom;

        GG_soft(c*(i-1)+1:c*i, p*(i-1)+1:p*i)    = G_soft_nom;
        GG_soft(c*(i-1)+1:c*i, Hp*p+1)           = g_soft_nom;
    end
    
    % Build constraint weights
    if strcmpi(mpc_input.solver, 'gurobi')
        eweight = mpc_input.constraints.weights.hard_rate;
        fweight = mpc_input.constraints.weights.hard_input;
        gweight = mpc_input.constraints.weights.hard_output;
        
        E_weights(a*(i-1)+1:a*i, 1) = reshape(eweight.',1,[])';
        F_weights(b*(i-1)+1:b*i, 1) = reshape(fweight.',1,[])';
        G_weights(c*(i-1)+1:c*i, 1) = reshape(gweight.',1,[])';
    end
end

% ------------------------------------------------------------------- %
% Calculate targets
% ------------------------------------------------------------------- %
LHS = [A, B; C, zeros(size(C,1), size(B,2))];
RHS = [zeros(size(f0)); r];

% TODO - RHS must be solver otherwise
% RHS = []

xu_inf = lsqminnorm(LHS, RHS);  % Optimal steady-state
% xu_inf = linsolve(LHS, RHS);  % Optimal steady-state
x_inf = xu_inf(1:n); % [5;0]; %
u_inf = xu_inf(n+1:end); % [0]; % 
X_inf = repmat(x_inf, [Hp,1]);
U_inf = repmat(u_inf, [Hp,1]);

% ------------------------------------------------------------------- %
% Build Prediction Matrices
% ------------------------------------------------------------------- %
AB_big      = zeros(n*Hp, m);
ABdu_big    = zeros(n*Hp, m*Hp);
T_fancy     = zeros(p*Hp, 1);
IA_big      = zeros(n*Hp, n);

for i = 1:Hp
    A_big(n*(i-1)+1:n*i, 1:n)           = Ad^i;
    C_big(p*(i-1)+1:p*i, n*(i-1)+1:n*i) = C;

    AB_sum = 0;
    for j = 0:i-1
        AB_sum                      = AB_sum + Ad^j * Bd;
    end
    AB_big(n*(i-1)+1:i*n, 1:m)      = AB_sum;

    ABdu_big(n*(i-1)+1:n*i, 1:m)    = AB_sum;
    if i > 1
        ABdu_big(n*(i-1)+1:n*i, m+1:m*Hp) = ABdu_big(n*(i-2)+1:n*(i-1), 1:m*(Hp-1));
    end

    T_fancy(p*(i-1)+1:p*i, 1)       = r;

    IA_sum = 0;
    for j = 0:i-1
        IA_sum                      = IA_sum + Ad^j;
    end
    IA_big(n*(i-1)+1:i*n, 1:n) = IA_sum;

    % Constraints - hard
    F_sum = zeros(b*Hp, m);
    for j = i:Hp
        F_sum                       = F_sum + FF(:, m*(j-1)+1:m*j);
    end
    F_fancy(:, m*(i-1)+1:m*i)       = F_sum;

    % Constraints - soft
    if strcmpi(mpc_input.solver, 'quadprog') && ...
        strcmpi(mpc_input.constraints.type, 'hard') ~= 1
        F_soft_sum = zeros(b*Hp, m);
        for j = i:Hp
            F_soft_sum                  = F_soft_sum + FF_soft(:, m*(j-1)+1:m*j);
        end
        F_fancy_soft(:, m*(i-1)+1:m*i)  = F_soft_sum;
    end
    
end

W           = EE(:, 1:end-1);
w           = EE(:, end);
f           = FF(:,end);
Tau         = GG(:, 1:end-1);
g           = GG(:, end);

if strcmpi(mpc_input.solver, 'quadprog') && ...
   strcmpi(mpc_input.constraints.type, 'hard') ~= 1
    W_soft      = EE_soft(:, 1:end-1);
    w_soft      = EE_soft(:, end);
    f_soft      = FF_soft(:,end);
    g_soft      = GG_soft(:, end);
    Tau_soft    = GG_soft(:, 1:end-1);
    F_1_soft    = F_fancy_soft(:, 1:m);
end

% Prediction matrices
LTM = build_LTM(m,Hp);
x_diff = X_inf - A_big*x0' - AB_big*u0';
u_diff = repmat(u_inf - u0, Hp, 1);


% Formulate Quadratic Program
% ----------------------------------
% Prediction matrices
G           = 2 * (ABdu_big' * Q_fancy * x_diff + LTM' * R_fancy * u_diff);
H           = ABdu_big' * Q_fancy * ABdu_big + LTM' * R_fancy * LTM;

% Constraint matrices - for hard constraints
F_1         = F_fancy(:, 1:m);
psi         = C_big*A_big;
gamma       = C_big*AB_big;
theta       = C_big*ABdu_big;
phi         = C_big*IA_big;

Omega       = [F_fancy; 
               Tau*theta; 
               W];
omega       = [-F_1*ubar_km1 - f;
               -Tau*(psi*xbar_k + gamma*ubar_km1 + phi*f0d + repmat(g0, Hp, 1)) - g;
               -w];

% Constraint matrices - for soft constraints
if strcmpi(mpc_input.solver, 'quadprog') && ...
   strcmpi(mpc_input.constraints.type, 'hard') ~= 1
    Omega_soft  = [F_fancy_soft; 
                      Tau_soft*theta; 
                      W_soft];

    omega_soft  = [-F_1_soft*ubar_km1 - f_soft;
                   -Tau_soft*(psi*xbar_k + gamma*ubar_km1 + phi*f0d +  ...
                      repmat(g0, Hp, 1)) - g_soft;
                   -w_soft];
end


% Solve Quadratic Program
% ----------------------------------
if strcmpi(mpc_input.solver, 'gurobi')
    % SOLVE WITH GUROBI
    
    % qp_solve integration - requires codegen of solve_qp to get solve_qp_mex
    nRows       = size(Omega,1);        % Number of constraints?
    nCols       = m*Hp;                 % dimension of A
    nInputs     = 1;                    % un-used int
%     constraint_weights = 1*ones(size(Omega,1),1);
    constraint_weights = [F_weights; G_weights; E_weights];
    Uopt_mat    = zeros(m*Hp,1);

    % IMPORTANT - the line below will only work if solve_qp_mex has been
    % compiled with the correct input argument sizes (see notes.txt)
%     cd '/home/kieran/Documents/MATLAB/MPC'
%     run_codegen
    [output, dUbar] = solve_qp_mex(nRows, nCols, nInputs, H, G, Omega, omega, ...
                                   constraint_weights, Uopt_mat);

    if output == 0
        return
    end

else
    % SOLVE WITH QUADPROG

    if strcmpi(constraints.type, 'none')
        dUbar   = quadprog(2*H, -G);

    elseif strcmpi(constraints.type, 'hard')
        [dUbar, ~, ~, ~, ~] = quadprog(2*H, -G, Omega, omega);

    elseif strcmpi(constraints.type, 'mixed')
            % Constraint matrices
        A1      = [[Omega, zeros(size(Omega, 1), size(omega, 1))]; ...
                   [Omega_soft, -eye(size(Omega_soft, 1), size(omega_soft, 1))]; ...
                   [zeros(size(omega, 1), Hp*m), -eye(size(omega,1))]];
        b1      = [omega; omega_soft; zeros(size(omega))];

        if strcmpi(penalty_method, 'quadratic')
                % Penalty matrices
            V       = penalty_weight*eye(size(omega, 1));
            H1      = blkdiag(2*H, V*eye(size(omega, 1)));
            f1      = [-G; zeros(size(omega, 1), 1)];

                % Solution
            X1      = quadprog(H1, f1, A1, b1);
            dUbar   = X1(1:Hp*m, 1);

        else
            % Soft constraint implementation - linear penalty
                % Penalty matrices
            H1      = blkdiag(2*H, zeros(size(omega, 1)));
            f1      = [-G; penalty_weight*eye(size(omega, 1), 1)];

                % Solution
            X1      = quadprog(H1, f1, A1, b1);
            dUbar   = X1(1:Hp*m, 1);

        end

    else
            % Constraint matrices
        A1      = [[Omega_soft, -eye(size(Omega_soft, 1), size(omega_soft, 1))]; ...
                   [zeros(size(omega, 1), Hp*m), -eye(size(omega,1))]];
        b1      = [omega_soft; zeros(size(omega))];

        if strcmpi(penalty_method, 'quadratic')
                % Penalty matrices
            V       = penalty_weight*eye(size(omega, 1));
            H1      = blkdiag(2*H, V*eye(size(omega, 1)));
            f1      = [-G; zeros(size(omega, 1), 1)];

                % Solution
            X1      = quadprog(H1, f1, A1, b1);
            dUbar   = X1(1:Hp*m, 1);

        else
            % Soft constraint implementation - linear penalty
                % Penalty matrices
            H1      = blkdiag(2*H, zeros(size(omega, 1)));
            f1      = [-G; penalty_weight*eye(size(omega, 1), 1)];

                % Solution
            X1      = quadprog(H1, f1, A1, b1);
            dUbar   = X1(1:Hp*m, 1);
        end

    end

end

% Get optimal control input
% ----------------------------------
Ubar                    = zeros(size(dUbar));
ubar_previous           = ubar_km1;
for h = 1:Hp
    delta_ubar              = dUbar(m*(h-1)+1:h*m,:);
    Ubar(m*(h-1)+1:h*m, 1)  = ubar_previous + delta_ubar;
    ubar_previous           = Ubar(m*(h-1)+1:h*m, 1);   
end
U_k     = Ubar + repmat(u0', Hp, 1);


% PLOTTING
ybar_predicted  = psi*xbar_k + gamma*ubar_km1 + ...
                  theta*dUbar + phi*f0d; % + repmat(g0, Hp, 1);
              
figure(4);
clf;

subplot(2,2,1);
hold on; grid on;
title('Flap Input');
stairs(ybar_predicted(2:p:end)*180/pi);

subplot(2,2,2);
hold on; grid on;
title('Thrust setting');
stairs(ybar_predicted(3:p:end));

subplot(2,2,3);
hold on; grid on;
title('Altitude');
plot(ybar_predicted(1:p:end));

subplot(2,2,4);
hold on; grid on;
title('FPA');
stairs(ybar_predicted(5:p:end));

output = U_k;

