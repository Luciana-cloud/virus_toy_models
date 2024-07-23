function [output] = calibration_kill_winner(p)

% Options for the solver
abstol = 1e-7;
reltol = 1e-5;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:4); 

% Time vector (h)
tspan   = [0 100];

% Parameters  
% p(1) = 0.7; % Growth rate coefficient for microbes [1/h]
% p(2) = 5;   % Saturation constant [µgC/ml]
% p(3) = 0.5; % Efficiency [-]
% p(4) = 0.15; % Growth rate coefficient for virus [1/h]
% p(5) = 0.1; % Dead rate coefficient for microbes [1/h]
% p(6) = 0.1; % Dead rate coefficient for virus [1/h]
% p(7) = 0.2; % Fraction of dead biomass sorbed onto soils [-]

% Fixed parameters
W_V = (0.055 + 0.2)*1e-9/2; % Weight of virus in C [µgC/particle] [0.055 - 0.2 fg/particle]
W_B = (50 + 250)*1e-9/2;    % Weight of Host in C (Bacteria as reference) [µgC/particle] [50 - 250 fg/cell]

% Set Initial Conditions
c(1) = W_B*10^8;                      % Microbial Biomass [µgC/ml]
c(2) = W_V*10^9;                      % Viral Biomass [µgC/ml] 
c(3) = 1000;                          % Soil organic matter (substrate) [µgC/ml] 1000
c(4) = 0;                             % CO2 [µgC/ml]

% Running the model
try
    warning off
    tic
    [ty,cu] = ode15s(@kill_winner_simple,tspan,c,o_opts,p'); 
 catch ME
     warning off
end
if length(cu) < length(tspan)
    cu = ones(length(tspan),length(c))*1e+99;
end

if isreal(cu)==0
    cu = ones(length(tspan),length(c))*1e+99;    
end

% Model Output

output = horzcat(ty,cu);

end