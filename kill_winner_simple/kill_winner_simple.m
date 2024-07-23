function vf_ = kill_winner_simple(t,x,p)

%% State Variables %%

C_B = x(1); % Microbial Biomass [µgC/ml]
C_V = x(2); % Viral Biomass [µgC/ml]
C_S = x(3); % Soil organic matter (substrate) [µgC/ml]
% CO2 = x(4); % CO2 [µgC/ml]

%% Model Parameters %%

mu_B = p(1); % Growth rate coefficient for microbes [1/h]
K_B  = p(2); % Saturation constant [µgC/ml]
Y    = p(3); % Efficiency [-]
mu_V = p(4); % Growth rate coefficient for virus [1/h]
d_B  = p(5); % Dead rate coefficient for microbes [1/h]
d_V  = p(6); % Dead rate coefficient for virus [1/h]
g    = p(7); % Fraction of dead biomass sorbed onto soils [-] 

%% Functions %%

uptake_B = mu_B*C_B*C_S*(K_B + C_S)^(-1);
growth_B = Y * uptake_B;
growth_V = mu_V*C_B*C_V;
decay_B  = d_B*C_B*C_V;
decay_V  = d_V*C_V;

%% ODE system %%

vf_      = zeros(4,1);
vf_(1)   = growth_B - decay_B; % Microbial Biomass [µgC/ml]
vf_(2)   = growth_V - decay_V; % Viral Biomass [µgC/ml]
vf_(3)   = - uptake_B + g*(decay_B + decay_V); % Soil organic matter (substrate) [µgC/ml]
vf_(4)   = (1 - Y)*uptake_B + (1 - g)*(decay_B + decay_V); % CO2 [µgC/ml]

end