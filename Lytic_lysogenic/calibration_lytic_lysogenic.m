function [output] = calibration_lytic_lysogenic(p)

% Options for the solver
abstol = 1e-7;
reltol = 1e-5;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:6); 

% Time vector (h)
tspan   = [0 100];

% Parameters  
% p(1)  = 0.7;  % Growth rate coefficient for non-infected microbes [1/h]
% p(2)  = 0.1;  % Growth rate coefficient for Lysogenic microbes [1/h]
% p(3)  = 5;    % Saturation constant for non-infected and Lysogenic microbes [µgC/ml]
% p(4)  = 0.3;  % Efficiency [-]
% p(5)  = 0.1;  % Dead rate coefficient for non-infected microbes [1/h]
% p(6)  = 0.2;  % Dead rate coefficient for Lysogenic microbes [1/h]
% p(7)  = 0.4;  % Dead rate coefficient for Lytic microbes [1/h]
% p(8)  = 0.2;  % Dead rate coefficient for virus [1/h]
% p(9)  = 1e-2; % Rate of absorption of virus population to hosts [1/h]
% p(10) = 1e-2; % Induction rate [1/h]
% p(11) = 0.2;  % Fraction of dead biomass sorbed onto soils [-] 
% p(12) = 1e-2; % Fraction of cells that undergo lysogenic infection [-]

% Fixed parameters
W_V = (0.055 + 0.2)*1e-9/2; % Weight of virus in C [µgC/particle] [0.055 - 0.2 fg/particle]
W_B = (50 + 250)*1e-9/2;    % Weight of Host in C (Bacteria as reference) [µgC/particle] [50 - 250 fg/cell]

% Set Initial Conditions
c(1) = W_B*10^8;                      % Non-infected Microbial Biomass [µgC/ml]
c(2) = 0;                             % Lysogenic Microbial Biomass [µgC/ml]
c(3) = 0;                             % Lytic Microbial Biomass [µgC/ml]
c(4) = W_V*10^9;                      % Viral Biomass [µgC/ml] 
c(5) = 1000;                          % Soil organic matter (substrate) [µgC/ml] 1000
c(6) = 0;                             % CO2 [µgC/ml]

% Running the model
try
    warning off
    tic
    [ty,cu] = ode15s(@lytic_lysogenic,tspan,c,o_opts,p'); 
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