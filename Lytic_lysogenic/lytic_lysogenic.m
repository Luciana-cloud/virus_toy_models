function vf_ = lytic_lysogenic(t,x,p)

%% State Variables %%

C_B = x(1); % Non-infected Microbial Biomass [µgC/ml]
C_L = x(2); % Lysogenic Microbial Biomass [µgC/ml]
C_I = x(3); % Lytic Microbial Biomass [µgC/ml]
C_V = x(4); % Viral Biomass [µgC/ml]
C_S = x(5); % Soil organic matter (substrate) [µgC/ml]
% CO2 = x(6); % CO2 [µgC/ml]

%% Model Parameters %%

mu_B  = p(1); % Growth rate coefficient for non-infected microbes [1/h]
mu_L  = p(2); % Growth rate coefficient for Lysogenic microbes [1/h]
K_B   = p(3); % Saturation constant for non-infected and Lysogenic microbes [µgC/ml]
Y     = p(4); % Efficiency [-]
d_B   = p(5); % Dead rate coefficient for non-infected microbes [1/h]
d_L   = p(6); % Dead rate coefficient for Lysogenic microbes [1/h]
d_I   = p(7); % Dead rate coefficient for Lytic microbes [1/h]
d_V   = p(8); % Dead rate coefficient for virus [1/h]
phi   = p(9); % Rate of absorption of virus population to hosts [1/h]
gamma = p(10); % Induction rate [1/h]
g     = p(11); % Fraction of dead biomass sorbed onto soils [-] 
alpha = p(12); % Fraction of cells that undergo lysogenic infection [-]

%% Functions %%

uptake_B  = mu_B*C_B*C_S*(K_B + C_S)^(-1);
uptake_L  = mu_L*C_L*C_S*(K_B + C_S)^(-1);
growth_B  = Y * uptake_B;
growth_L  = Y * uptake_L;
decay_B   = d_B*C_B;
decay_L   = d_L*C_L;
decay_I   = d_I*C_I*C_V;
decay_V   = d_V*C_V;
infection = phi*C_B*C_V;
induction = gamma*C_L;

%% ODE system %%

vf_      = zeros(6,1);
vf_(1)   = growth_B - infection - decay_B; % Non-infected Microbial Biomass [µgC/ml]
vf_(2)   = alpha*infection + growth_L - induction - decay_L; % Lysogenic Microbial Biomass [µgC/ml]
vf_(3)   = induction + (1 - alpha)*infection - decay_I; % Lytic Microbial Biomass [µgC/ml]
vf_(4)   = decay_I - decay_V; % Viral Biomass [µgC/ml]
vf_(5)   = - uptake_B - uptake_L + g*(decay_B + decay_L + decay_V); % Soil organic matter (substrate) [µgC/ml]
vf_(6)   = (1 - Y)*uptake_B + (1 - Y)*uptake_L + (1 - g)*(decay_B + ...
    decay_L + decay_V); % CO2 [µgC/ml]

end