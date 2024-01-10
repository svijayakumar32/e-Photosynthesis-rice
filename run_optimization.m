%Eidata=importdata('Einput_rice_new.txt'); % Load enzyme activity levels optimized at 280 ppm for rice
%Eidata=importdata('Einput7.txt'); % Load enzyme activity levels optimized at 280 ppm
Ei=Eidata.data;
Einput=ones(37,1); % 
CO2i = 140; % Set Ci
PPFDi = 2000; % Set light intensity
WeatherTemp = 25; % Set temperature
global Vrubusco_adj; %reset enzyme activity
Vrubusco_adj = 1.12;%set 1.0
global VmaxAdj;% reset enzyme activity
VmaxAdj = 1.24;%set 1.0
global pcfactor;
pcfactor=1;
% Create output vector for A 
%NetAssimilationRate = zeros(1,1);
% Loop function to calculate net assimilation rate A for a given Ci value using set of enzyme activity levels optimized for 280
% for i=1:9 % No. of A values
% CO2i = (200:20:360)'; % Set Ci values ranging from 200 - 360
NetAssimilationRate= EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,0,0,Ei);
%A =  EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,0,0,Ei); 
%end
% Save assimilation rates to output text file
% file=fopen('output_A.txt','w');
% fprintf(file,'%.4f\n',NetAssimilationRate);
%writematrix does not work for older versions of MATLAB like 2018a on the HEC
%writematrix(NetAssimilationRate,'output_A.txt','Delimiter','space');