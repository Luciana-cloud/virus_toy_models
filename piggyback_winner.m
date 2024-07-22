function vf_ = piggyback_winner(x,p,q)

%% State Variables %%

C_B_a = x(1); % Active Microbial Biomass [mgC/ml]
C_B_i = x(2); % Dormant Microbial Biomass [mgC/ml]
C_V   = x(3); % Viral Biomass [mgC/ml]
C_S   = x(4); % Soil organic matter (substrate) [mgC/ml]
% CO2 = x(5); % CO2 [mgC/ml]

%% Model Parameters %%

mu_B  = p(1); % Growth rate coefficient for microbes [1/d]
K_B   = p(2); % Saturation constant [mgC/ml]
Y     = p(3); % Efficiency [-]
mu_V  = p(4); % Growth rate coefficient for virus [1/d]
d_B   = p(5); % Dead rate coefficient for microbes [1/d]
d_V   = p(6); % Dead rate coefficient for virus [1/d]
g     = p(7); % Fraction of dead biomass sorbed onto soils [-]
C_V_T = p(8); % Threshold concentration of virus [mgC/ml]

%% Model Fixed Parameters %%

n    = q(1); % Shape parameter [-]
K_a  = q(2); % Coefficient of activation [1/d]
K_i  = q(3); % Coefficient of deactivation [1/d]

%% Functions %%

uptake_B     = mu_B*C_B_a*C_S*(K_B + C_S)^(-1);
growth_B     = Y * uptake_B;
growth_V     = mu_V*C_B_a*C_V;
decay_B      = d_B*C_B_a*C_V;
decay_V      = d_V*C_V;
tau          = 1/(exp((C_V_T - C_V)/(n*C_V_T)) + 1);
activation   = tau*K_a*C_B_i;
deactivation = (1 - tau)*K_i*C_B_a;

%% ODE system %%

vf_      = zeros(5,1);
vf_(1)   = growth_B - decay_B + activation - deactivation; % Active Microbial Biomass [mgC/ml]
vf_(2)   = deactivation - activation; % Dormant Microbial Biomass [mgC/ml]
vf_(3)   = growth_V - decay_V; % Viral Biomass [mgC/ml]
vf_(4)   = uptake_B + g*(decay_B + decay_V); % Soil organic matter (substrate) [mgC/ml]
vf_(5)   = (1 - Y)*uptake_B - (1 - g)*(decay_B + decay_V); % CO2 [mgC/ml]

end