function output = forward_simulate(sim_input)

% ----------------------------------------------------------------------- %
% Deconstruct input
% ----------------------------------------------------------------------- %
plant_model             = sim_input.plant_model;
input                   = sim_input.dynamics_input;
input.auxdata           = plant_model.auxdata;

plant_dynamics          = plant_model.dynamics;
plant_output            = plant_model.output;

% ----------------------------------------------------------------------- %
% Initialise loop
% ----------------------------------------------------------------------- %
xk1                     = input.phase.state;
Ts                      = sim_input.Ts;

% ----------------------------------------------------------------------- %
% Forward simulate
% ----------------------------------------------------------------------- %
for j = 1:100
    input.phase.state   = xk1;
    dynamics            = plant_dynamics(input);
    f_dot               = dynamics.dynamics;
    xk1                 = xk1 + f_dot*Ts/100;
end

zk                      = plant_output(xk1);

% ----------------------------------------------------------------------- %
% Package output
% ----------------------------------------------------------------------- %
output.zk1              = zk;
output.xk1              = xk1;
