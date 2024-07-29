%% Kill the winner simple model %%

load('kill_winner_simple.mat')

[output] = calibration_kill_winner_simple(p);

W_V = (0.055 + 0.2)*1e-9/2; % Weight of virus in C 
W_B = (50 + 250)*1e-9/2;    % Weight of Host in C 

Time      = output(:,1);
Microbes  = output(:,2)/W_B;
Virus     = output(:,3)/W_V;
SOM       = output(:,4);
CO2       = output(:,5);
decay_B   = p(5).*output(:,2).*output(:,3);
decay_V   = p(6).*output(:,3);
Necromass = p(7).*(decay_B + decay_V);
Ratio     = Virus./Microbes;

fig = figure;
subplot(2,2,1);
plot(Time,Microbes)
ylabel('Microbial cells per ml')
hold on
yyaxis right
plot(Time,Virus)
xlabel('Time (hr)')
ylabel('Virus particles per ml')

subplot(2,2,2);
plot(Time,Ratio)
xlabel('Time (hr)')
ylabel('Virus - to - Microbes Ratio')

subplot(2,2,3);
plot(Time,SOM)
xlabel('Time (hr)')
ylabel('Carbon Substrate (µgC/ml)')
yyaxis right
plot(Time,Necromass)
ylabel('Necromass (µgC/ml)')

subplot(2,2,4);
plot(Time,CO2)
xlabel('Time (hr)')
ylabel('CO2')

saveas(fig,'C:\luciana_datos\UCI\Project_12 (Fire modeling project)\virus_paper\virus_toy_models\Figures\kill_winner_simple.png')

%% Piggyback the winner %%

load('piggyback_winner.mat')

[output] = calibration_piggyback_winner(p);

W_V = (0.055 + 0.2)*1e-9/2; % Weight of virus in C 
W_B = (50 + 250)*1e-9/2;    % Weight of Host in C 

Time             = output(:,1);
Active_Microbes  = output(:,2)/W_B;
Dormant_Microbes = output(:,3)/W_B;
Virus            = output(:,4)/W_V;
SOM              = output(:,5);
CO2              = output(:,6);
decay_B          = p(5).*output(:,2);
decay_V          = p(6).*output(:,4);
Necromass        = p(7).*(decay_B + decay_V);
Ratio            = Virus./(Active_Microbes+Dormant_Microbes);

fig1 = figure;
subplot(2,2,1);
A = plot(Time,Active_Microbes);
hold on
B = plot(Time,Dormant_Microbes);
ylabel('Microbial cells per ml')
yyaxis right
C = plot(Time,Virus);
xlabel('Time (hr)')
ylabel('Virus particles per ml')
hold off
legend([A,B,C],'Active','Dormant','Virus')

subplot(2,2,2);
plot(Time,Ratio)
xlabel('Time (hr)')
ylabel('Virus - to - Microbes Ratio')

subplot(2,2,3);
plot(Time,SOM)
xlabel('Time (hr)')
ylabel('Carbon Substrate (µgC/ml)')
yyaxis right
plot(Time,Necromass)
ylabel('Necromass (µgC/ml)')

subplot(2,2,4);
plot(Time,CO2)
xlabel('Time (hr)')
ylabel('CO2')

saveas(fig1,'C:\luciana_datos\UCI\Project_12 (Fire modeling project)\virus_paper\virus_toy_models\Figures\piggyback_winner.png')

%% Lytic-lysogenic %%

load('Lytic_lysogenic.mat')

[output] = calibration_lytic_lysogenic(p);

W_V = (0.055 + 0.2)*1e-9/2; % Weight of virus in C 
W_B = (50 + 250)*1e-9/2;    % Weight of Host in C 

Time                   = output(:,1);
Non_infected_Microbes  = output(:,2)/W_B;
Lysogenic_Microbes     = output(:,3)/W_B;
Lytic_Microbes         = output(:,4)/W_B;
Virus                  = output(:,5)/W_V;
SOM                    = output(:,6);
CO2                    = output(:,7);
decay_B                = p(5).*output(:,2);
decay_L                = p(6).*output(:,3);
decay_V                = p(8).*output(:,5);
Necromass              = p(7).*(decay_B + decay_L + decay_V);
Ratio                  = Virus./(Non_infected_Microbes+...
    Lysogenic_Microbes+Lytic_Microbes);

fig2 = figure;
subplot(2,2,1);
A = plot(Time,(Non_infected_Microbes));
hold on
B = plot(Time,(Lysogenic_Microbes));
C = plot(Time,(Lytic_Microbes));
ylabel('Microbial cells per ml')
yyaxis right
D = plot(Time,(Virus));
xlabel('Time (hr)')
ylabel('Virus particles per ml')
hold off
legend([A,B,C],'non-Infected','Lysogenic','Lytic')

subplot(2,2,2);
plot(Time,Ratio)
xlabel('Time (hr)')
ylabel('Virus - to - Microbes Ratio')

subplot(2,2,3);
plot(Time,SOM)
xlabel('Time (hr)')
ylabel('Carbon Substrate (µgC/ml)')
yyaxis right
plot(Time,Necromass)
ylabel('Necromass (µgC/ml)')

subplot(2,2,4);
plot(Time,CO2)
xlabel('Time (hr)')
ylabel('CO2')

saveas(fig2,'C:\luciana_datos\UCI\Project_12 (Fire modeling project)\virus_paper\virus_toy_models\Figures\Lytic_lysogenic.png')

%%
