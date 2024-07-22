function vf_ = kill_winner(t,x,p)

%% State Variables %%

C_B = x(1); % Free Biomass [mgC/ml]
C_L = x(2); % Lytic Biomass [mgC/ml]
C_V = x(3); % Viral Biomass [mgC/ml]
C_S = x(4); % Soil organic matter (substrate) [mgC/ml]
% CO2 = x(5); % CO2 [mgC/ml]

%% Model Parameters %%

mu_B = p(1); % Growth rate coefficient for microbes [1/d]
K_B  = p(2); % Saturation constant [mgC/ml]
Y    = p(3); % Efficiency [-]
d_B  = p(4); % Dead rate coefficient for microbes [1/d]
d_V  = p(5); % Dead rate coefficient for virus [1/d]
g    = p(6); % Fraction of dead biomass sorbed onto soils [-] 
phi  = p(7); % Rate of absorption of virus population to hosts [1/d]

%% Functions %%

uptake_B  = mu_B*C_B*C_S*(K_B + C_S)^(-1);
growth_B  = Y * uptake_B;
decay_B   = d_B*C_B;
decay_L   = d_L*C_L*C_V;
decay_V   = d_V*C_V;
infection = phi*C_B*C_V;

%% ODE system %%

vf_      = zeros(5,1);
vf_(1)   = growth_B - infection - decay_B; % Free Microbial Biomass [mgC/ml]
vf_(2)   = infection - decay_L;            % Lytic Microbial Biomass [mgC/ml]
vf_(3)   = decay_L - decay_V; % Viral Biomass [mgC/ml]
vf_(4)   = - uptake_B + g*(decay_B + decay_V); % Soil organic matter (substrate) [mgC/ml]
vf_(5)   = (1 - Y)*uptake_B + (1 - g)*(decay_B + decay_V); % CO2 [mgC/ml]

end