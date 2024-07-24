function [output] = calibration_lytic_lysogenic(p)

% Options for the solver
abstol = 1e-7;
reltol = 1e-5;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:6); 

% Time vector (h)
tspan   = [0 100];

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