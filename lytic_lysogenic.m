function vf_ = lytic_lysogenic(x,p,q)

%% State Variables %%

C_B = x(1); % Non-infected Microbial Biomass [mgC/ml]
C_L = x(2); % Lysogenic Microbial Biomass [mgC/ml]
C_I = x(3); % Lytic Microbial Biomass [mgC/ml]
C_V = x(4); % Viral Biomass [mgC/ml]
C_S = x(5); % Soil organic matter (substrate) [mgC/ml]
% CO2 = x(6); % CO2 [mgC/ml]

%% Model Parameters %%

mu_B  = p(1); % Growth rate coefficient for non-infected microbes [1/d]
mu_L  = p(2); % Growth rate coefficient for Lysogenic microbes [1/d]
K_B   = p(3); % Saturation constant for non-infected microbes [mgC/ml]
K_L   = p(4); % Saturation constant for Lysogenic microbes [mgC/ml]
Y     = p(5); % Efficiency [-]
mu_V  = p(6); % Growth rate coefficient for virus [1/d]
d_B   = p(7); % Dead rate coefficient for non-infected microbes [1/d]
d_L   = p(8); % Dead rate coefficient for Lysogenic microbes [1/d]
d_I   = p(9); % Dead rate coefficient for Lytic microbes [1/d]
d_V   = p(10); % Dead rate coefficient for virus [1/d]
phi   = p(11); % Rate of absorption of virus population to hosts [1/d]
gamma = p(12); % Induction rate [1/d]
g     = p(13); % Fraction of dead biomass sorbed onto soils [-] 

%% Functions %%

uptake_B  = mu_B*C_B*C_S*(K_B + C_S)^(-1);
uptake_L  = mu_L*C_L*C_S*(K_L + C_S)^(-1);
growth_B  = Y * uptake_B;
growth_L  = Y * uptake_L;
growth_V  = mu_V*C_B*C_V;
decay_B   = d_B*C_B;
decay_L   = d_L*C_L;
decay_I   = d_I*C_I*C_V;
decay_V   = d_V*C_V;
infection = phi*C_B*C_V;
induction = gamma*C_L;

%% ODE system %%

vf_      = zeros(6,1);
vf_(1)   = growth_B - infection - decay_B; % Non-infected Microbial Biomass [mgC/ml]
vf_(2)   = alpha*infection + growth_L - induction - decay_L; % Lysogenic Microbial Biomass [mgC/ml]
vf_(3)   = induction + (1 - alpha)*infection - decay_I; % Lytic Microbial Biomass [mgC/ml]
vf_(4)   = decay_I - decay_V; % Viral Biomass [mgC/ml]
vf_(5)   = - uptake_B - uptake_L + g*(decay_B + decay_L + decay_V); % Soil organic matter (substrate) [mgC/ml]
vf_(6)   = (1 - Y)*uptake_B + (1 - Y)*uptake_L + (1 - g)*(decay_B + ...
    decay_L + decay_V); % CO2 [mgC/ml]

end