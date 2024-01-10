%Function file for run_optimization
%Eidata_1=importdata('Einput_potato_140_RSA.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_2=importdata('Einput_potato_160.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_3=importdata('Einput_potato_180.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_4=importdata('Einput_potato_200.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_5=importdata('Einput_potato_220.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_6=importdata('Einput_potato_240.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_7=importdata('Einput_potato_260.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_8=importdata('Einput_potato_280.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_9=importdata('Einput_potato_300.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_10=importdata('Einput_potato_320.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_11=importdata('Einput_potato_340.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_12=importdata('Einput_potato_360.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_13=importdata('Einput_potato_380.txt'); % Load calibrated enzyme activity levels for potato
% Eidata_14=importdata('Einput_potato_400.txt'); % Load calibrated enzyme activity levels for potato
Eidata_15=importdata('Einput_potato_420_Rubisco_AGPase_SPS.txt'); % Load calibrated enzyme activity levels for potato
%Eidata=importdata('Einput_potato.txt'); % Load calibrated enzyme activity levels for potato
%Eidata=importdata('Einput_rice.txt'); % Load calibrated enzyme activity levels for rice
Ei=Eidata_15.data;
%Ei([7,9,12],:) = 0;%Remove enzymes omitted during optimization  %18.05.23%

Ei(7) = Ei(4); % Ensure double-counted enzymes have the same activity i.e. V8=V5 and V10=V7
Ei(9) = Ei (6);

% Ei(12) = 0; Leave ATP synthase as original value or set to 0?
%Ei=Eidata;
Einput=ones(37,1); % 
%CO2i = 140; % Set single Ci if not using loop
PPFDi = 2000; % Set light intensity
WeatherTemp = 25; % Set temperature
global Vrubusco_adj; 
Vrubusco_adj = 1.0;
global VmaxAdj;
VmaxAdj = 1.0;
global pcfactor;
pcfactor=1;
% Create output vector for A to loop over multiple Ci's
NetAssimilationRate = zeros(15,1);
% Loop function to calculate net assimilation rate A for a range of Ci values using set of enzyme activity levels 
for i=1:15 % No. of A values
CO2i = (140:20:420)'; % Set Ci values ranging from 140-420 with stepsize of 20
NetAssimilationRate(i)= EPS_Drive_GRNs(Einput,CO2i(i,1),PPFDi,WeatherTemp,0,0,Ei);
%A1 = EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,0,0,Ei); 
end
% Save assimilation rates to output text file
%file=fopen('output_A.txt','w');
%fprintf(file, %6.2f %12.8f\r\n', NetAssimilationRate);
%writematrix(NetAssimilationRate,'output_A.txt','Delimiter','space');