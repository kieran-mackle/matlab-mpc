function output = mpc_control(input)
% Applies MPC to a plant input. The input is a struct containing the
% following keys:
%   input
%     |- mpc_input
%     |    |- control_model
%     |    |- cost
%     |    |    |- output
%     |    |    |- control
%     |    |- constraints
%     |    |    |- hard
%     |    |    |- soft
%     |    |    |- type
%     |    |- penalty_method
%     |    |- penalty_weight
%     |    |- params
%     |    |- initial
%     |- sim_input
%     |    |- dynamics_input
%     |    |- plant_model
%     |    |- Ts
%     |- reference_function
%

% ----------------------------------------------------------------------- %
% Deconstruct input
% ----------------------------------------------------------------------- %
mpc_input       = input.mpc_input;
sim_input       = input.sim_input;
initial         = mpc_input.initial;
control_output  = mpc_input.control_model.output;
get_reference   = input.reference_function;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
Ts              = mpc_input.params.timestep;      % Sampling interval
total_time      = mpc_input.params.sim_time;      % Total time horizon

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0              = initial.state;
u0              = initial.control;
z0              = control_output(x0);

% ----------------------------------------------------------------------- %
% Simulate trajectory
% ----------------------------------------------------------------------- %
N               = round(total_time/Ts);
t               = 0;
X               = x0';
U               = u0';
Z               = z0';
mpc_input       = get_reference(mpc_input);
R               = [mpc_input.reference'];

% figure(1); clf; title('Altitude (km)'); hold on; grid on;
% figure(2); clf; title('Flap deflection (deg)'); hold on; grid on;
% figure(3); clf; title('Thrust setting'); hold on; grid on;

for k = 1:N
    % Update inputs with latest information
    % --------------------------------------
    mpc_input.initial           = initial;
    mpc_input.t                 = t(k);
    mpc_input                   = get_reference(mpc_input); % Updates reference 
                                                            % and cost matrices
%     
%     % Normalise model before MPC update
% %     % ----------------------------------
%     % STATE VARIABLES
%     scaled_state                = scale_input(mpc_input.initial.state, ...
%                                               mpc_input.scaling.state);
%     scaled_initial.state        = scaled_state;
%     
%     % REFERENCE
%     scaled_reference            = scale_input(mpc_input.reference, ...
%                                               mpc_input.scaling.reference);
%     
%     % CONTROL INPUTS
%     scaled_control              = scale_input(mpc_input.initial.control, ...
%                                               mpc_input.scaling.control);
%     scaled_initial.control      = scaled_control;
%     
%     % CONSTRAINTS
%     
%     
%     % ANYTHING ELSE
%     % Auxdata?
%     
%     mpc_input.initial           = scaled_initial;
%     mpc_input.reference         = scaled_reference;
    
    % Run MPC update
    % ----------------------------------
    U_k                         = mpc_update(mpc_input);
    uk                          = U_k(1:length(u0), :);
    
%     % Scale MPC output
%     % ----------------------------------
%     uk                          = normalise_input(uk, mpc_input.scaling.control);
    
    % Run forward simultation
    % ----------------------------------
    dyn_input.phase.time        = t(k);
    dyn_input.phase.state       = initial.state;
    dyn_input.phase.control     = uk';
    sim_input.dynamics_input    = dyn_input;
    sim_update                  = forward_simulate(sim_input);
    
    % Update state struct
    % ----------------------------------
    initial.state               = sim_update.xk1;
    initial.control             = uk';
    
    % Append results to output arrays
    % ----------------------------------
    t                           = [t, Ts*k];
    U                           = [U,uk];
    X                           = [X, sim_update.xk1'];
    Z                           = [Z, sim_update.zk1'];
    R                           = [R, mpc_input.reference'];
    
%     if mod(t(end), 5) == 0
%         figure(1); plot(t(1:k+1), Z(1, 1:k+1)*1e-3, 'k-');
%         figure(2); stairs(t(1:k+1), Z(2, 1:k+1)*180/pi, 'k-');
%         figure(3); stairs(t(1:k+1), Z(3, 1:k+1), 'k-');
%     end
    
end

% Package output
% ----------------------------------
output.time     = t;
output.state    = X;
output.control  = U;
output.Z        = Z;
output.ref      = R;

