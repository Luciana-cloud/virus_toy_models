function [output] = calibration_kill_winner(p)

% Options for the solver
abstol = 1e-7;
reltol = 1e-5;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:4); 

% Time vector (d)
tspan   = [0 1000];

% Set Initial Conditions
c(1) = 100;                           % Microbial Biomass [mgC/ml]
c(2) = 10;                            % Viral Biomass [mgC/ml]
c(3) = 1000;                          % Soil organic matter (substrate) [mgC/ml]
c(4) = 0;                             % CO2 [mgC/ml]

% Running the model
try
    warning off
    tic
    [ty,cu] = ode15s(@kill_winner_simple,tspan,c,o_opts,p'); 
 catch ME
     warning off
end
if length(cu) < length(t)
    cu = ones(length(t),length(c))*1e+99;
end

if isreal(cu)==0
    cu = ones(length(t),length(c))*1e+99;    
end

% Model Output

Microbes = cu(:,1);
Virus    = cu(:,2);
SOM      = cu(:,3);
CO2      = cu(:,4);

Ratio = Virus./Microbes;

plot(ty,Microbes)
hold on
plot(ty,Virus)

end