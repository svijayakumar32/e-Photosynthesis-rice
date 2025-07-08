% Load temperature file
load('WeatherTemp.mat');

% Load Einput files 
Eidata = importdata('Einput7.txt');

% Specify alpha values to adjust enzyme activity levels here 
global Vrubusco_adj;
Vrubusco_adj = 1.36; 
global VmaxAdj;
VmaxAdj = 1.12;
global pcfactor;
pcfactor = 1;

% Adjust enzymes to calibrated enzyme activity levels (Vmax) for rice
Eidata.data(1,1) = Eidata.data(1,1) * Vrubusco_adj;
Eidata.data(2:26,1) = Eidata.data(2:26,1) * VmaxAdj;
Ei = Eidata.data;

% Ensure double-counted enzymes have the same activity i.e. V8 = V5 and V10 = V7
Ei(7) = Ei(4); 
Ei(9) = Ei(6);

% Load average optimized Vmax modelled under Cc = 130, 250 and 360 to create stacked enzyme scenarios
file = 'Outputs/rice_params/Results_optimization_rice_new_3.xlsx';

% Load optimized Vmax data from sheets for 130, 250 and 360 optimizations
optimised_Vmax_130 = readmatrix(file, 'Sheet', '130', 'Range', 'G5:G30');
optimised_Vmax_250 = readmatrix(file, 'Sheet', '250', 'Range', 'G5:G30');
optimised_Vmax_360 = readmatrix(file, 'Sheet', '360', 'Range', 'G5:G30');

% Generate 20% increase for comparison of enzyme stacking for Rubisco
%Vmax_plus_20 = Ei*1.2; % Salesse-Smith et al 2024 %OR%

% Generate 25% increase for comparison of enzyme stacking for Rubisco
Vmax_plus_25 = Ei*1.25; % Yoon et al 2020

% Generate twofold increase for comparison of enzyme stacking for all other enzymes
twofold_Vmax = Ei*2; 

% To model low Cc under drought, optimise only Rubisco and SBPase at 130 ppm, keep all other enzyme levels non-optimized
Eidata_low = vertcat(optimised_Vmax_130(1),Ei(2:6),Ei(4),optimised_Vmax_130(8),Ei(6),Ei(10:end));
Eidata_low_2 = vertcat(Vmax_plus_25(1),Ei(2:6),Ei(4),twofold_Vmax(8),Ei(6),Ei(10:end));

% To model ambient Cc, optimise Rubisco, SBPase, Aldolase and PRK at Cc = 250 ppm
Eidata_ambient = vertcat(optimised_Vmax_250(1),Ei(2:3),optimised_Vmax_250(4),Ei(5:6),optimised_Vmax_250(4),optimised_Vmax_250(8),Ei(6),optimised_Vmax_250(10),Ei(11:end));
Eidata_ambient_2 = vertcat(Vmax_plus_25(1),Ei(2:3),twofold_Vmax(4),Ei(5:6),twofold_Vmax(4),twofold_Vmax(8),Ei(6),twofold_Vmax(10),Ei(11:end));

% To model elevated Cc, optimise Rubisco, SBPase, Aldolase, PRK, FBPase and TK at Cc = 360 ppm
Eidata_elevated = vertcat(optimised_Vmax_360(1),Ei(2:3),optimised_Vmax_360(4:6),optimised_Vmax_360(4),optimised_Vmax_360(8),optimised_Vmax_360(6),optimised_Vmax_360(10),Ei(11:end));
Eidata_elevated_2 = vertcat(Vmax_plus_25(1),Ei(2:3),twofold_Vmax(4:6),twofold_Vmax(4),twofold_Vmax(8),twofold_Vmax(6),twofold_Vmax(10),Ei(11:end)); 

Einput=ones(37,1); % 
PPFDi = 2000; % Set light intensity
GRNC=0;

% Create output vectors for A to loop over multiple Cc's
GrossAssimilationRate = zeros(26,1);
GrossAssimilationRate_low = zeros(26,1);
GrossAssimilationRate_ambient = zeros(26,1);
GrossAssimilationRate_elevated = zeros(26,1);
GrossAssimilationRate_low_2 = zeros(26,1);
GrossAssimilationRate_ambient_2 = zeros(26,1);
GrossAssimilationRate_elevated_2 = zeros(26,1);

% Loop function to calculate gross assimilation rate A for a range of Cc values using set of enzyme activity levels 
for i = 1:26 % No. of A values
    CO2i = (130:10:380)'; % Set Cc values ranging from 130-380 with stepsize of 10
    % Calculate assimilation if using Vmax obtained from optimising
    % e-Photosynthesis for a given no. of enzymes (low = 2, ambient = 4, elevated = 6)
    GrossAssimilationRate(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Ei);%0
    GrossAssimilationRate_low(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_low);%2
    GrossAssimilationRate_ambient(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_ambient);%4
    GrossAssimilationRate_elevated(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_elevated);%6

    % Calculate assimilation if using twofold increases in Vmax
    GrossAssimilationRate_low_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_low_2);%2
    GrossAssimilationRate_ambient_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_ambient_2);%4
    GrossAssimilationRate_elevated_2(i) = EPS_Drive_GRNs(Einput,CO2i(i),PPFDi,WeatherTemp,GRNC,0,Eidata_elevated_2);%6
end

% Save assimilation rates to output text file
%file=fopen('non_optimized_A.txt','w');
%fprintf(file, %6.2f %12.8f\r\n', GrossAssimilationRate);
writematrix(GrossAssimilationRate,'non_optimized_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_low,'low_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_ambient,'ambient_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_elevated,'elevated_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_low_2,'low2_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_ambient_2,'ambient2_A.txt','Delimiter','space');
writematrix(GrossAssimilationRate_elevated_2,'elevated2_A.txt','Delimiter','space');