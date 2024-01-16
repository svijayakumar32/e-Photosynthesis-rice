clear all;
% REMEMBER TO CHECK IF USING KC OR KCAIR TO CALCULATE ASSIMILATION AT DIFF CI (LINE 98/99)
%Estimated From ACi data from Liana Acevedo-Siaca (2021), 
% Values presented in paper were fitted using Bernacchi and Long (2003) who don't estimate
% Chl level parameters for Rd and Gamma at 25C 
% Therefore O.sativa upper canopy data points were re-estimated from Fig 6,
% using Tleaf average of 6 replicates = 32.35C (Canopy Compiled Traits) and pressure of 100.5kPa 
% These values were used to construct and fit A/Ci and derive new parameters at Tleaf 
% corrected to equivalent values at 25C using constants from the Sharkey et al (2015) tool

% However this didn't work out so tried refitting WT IR64 data provided by IRRI
% Tleaf average of 8 replicates = 28.93104073C
Vcmax_m=134.380530667391;% Avg Estimated using msuRACiFit (Gregory et al 2016) with rice Gamma Star at 28.9C
J=189.091184437385;% Avg Estimated using msuRACiFit (Gregory et al 2016) with rice Gamma Star at 28.9C

%%%%%%%%%%
Lii=2000;%light intensity from IRRI
%Lii=1500;%light intensity from Acevedo-Siaca

%Farquhar model parameters
%Gr=38.6;%µbar von caemmerer 2020 %Gamma Star
Gr=50.2248975;%µbar %Gamma Star calculated from fitting IR64 of IRRI using Hermida-Carrera rice parameters at Tleaf
%Rd=1; %default value for day respiration (in light)
Rd = 1.33898789673117; %Rd estimated using msuRACiFit (Gregory et al 2016) rice parameters at Tleaf 
Gm = 5.395926789339;

I2=Lii/2*0.85*(1-0.15);
Theta=0.7;

Kc=348.318266403677;%ubar, 34.832 Pa, avg of Tleaf values, rice c and dHa
Ko=275.324214275076;%mbar, 27.532 kPa, avg of Tleaf values, rice c and dHa
Kc_air = 61.3901254795337;%ubar, 61.390 Pa avg of Tleaf values, rice c and dHa
O=210;%mbar

%Rearrange
%Kc_air=Kc*(1+(O/Ko))
%Ko=O/((Kc_air/Kc)-1)
%%%%%%%%%%%%%%%%%%%%%
% 05/01/2024 Smaller range to avoid intersection of Rubisco-limited and electron transport
% limited points
%CA=[100,150,200,250,300]; %Red
CA=[100,150,200,250];

global Vrubusco_adj;
global VmaxAdj;
VmaxAdj=1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj from Jmax_adj % default 1.3
global pcfactor;  
ProteinTotalRatio=0;
pcfactor=1/ProteinTotalRatio;
% inE=importdata('MeM_input5_0.txt');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Input conversion
% InputF=importdata('GrCM_output.txt');
% DataCol=InputF.textdata(1,:);
% idx = find(ismember(DataCol, 'Ele:Amb' ));
% DataGeneID=InputF.textdata(:,1);
% GeneNo=size(inE.textdata,1)-1;
% ExpValue=ones(GeneNo,1);
% for i=1:GeneNo
%     if find(ismember(DataGeneID,string(inE.textdata(i+1,1))))
%     idRow=find(ismember(DataGeneID,string(inE.textdata(i+1,1))));
%     ExpValue(i,1)=InputF.data(idRow-1,idx-2);
%     end
% end
% Einput=ExpValue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Einput=ones(37,1);%No gene expression data input
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
MetaOnly=0;% if MetaOnly=1 run only Metabolic model
WeatherTemp=28.9310407291759; %Original = 25C
WeatherRH=0.6;
WeatherWind=5;
Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Radiation_PAR=Lii/Convert*0.85*0.85;%%%%
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=0;

for j=1:25
    j
Vrubusco_adj=1.0+j*0.02;%adjust enzyme activity between 0.6-1.0 %default 0.8
Eio(1)=Edata.data(1,1)*Vrubusco_adj;
Eio(2:27)=Edata.data(2:27,1)*VmaxAdj;
Ci_vals=zeros(4,1);

for i=1:4
Air_CO2=CA(i);
if MetaOnly==1
CO2i=Air_CO2*0.7; % intercellular CO2 
PPFDi=Lii;
NetAssimilation=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eio);
else
LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
Ci=LeafResult(1);
Ci_vals(i)=Ci;
NetAssimilation=LeafResult(2);
Gs=LeafResult(3);
LeafTemperature=LeafResult(4);
Transpirationi=LeafResult(5);
end
%Calculate measured A values at different Ci
%ACI_m=Vcmax_m*(Ci-Gr)/(Ci+Kc*(1+O/Ko))-Rd; %sub gm to calculate Cc
%ACI_m=Vcmax_m*(Ci-Gr)/(Ci+Kc_air)-Rd; %subbed in Kc_air for Kc*(1+O/Ko)

%05/01/2024 modification
b=Vcmax_m-Rd+(Ci+Kc_air)*Gm;
c=((Ci-Gr)*Vcmax_m-(Ci+Kc_air)*Rd)*Gm;

ACI_m=(Vcmax_m-Rd+(Ci+Kc_air)*Gm)-(sqrt((Vcmax_m-Rd+(Ci+Kc_air)*Gm)^2)-(4*(((Ci-Gr)*Vcmax_m-(Ci-Kc_air)*Rd)*Gm)))/2;
ACi_evsm(i)=(ACI_Ac-NetAssimilation)^2;%the squares of the residuals
%ACi_evsm(i)=(ACI_m-NetAssimilation)^2;%the squares of the residuals

end
SSR(j,1)=Vrubusco_adj;
SSR(j,2)=sum(ACi_evsm);%the sum of the squares of the residuals
end
% fileID = fopen('LeafmetaOut.txt','w');
% fprintf(fileID,'%6s\n','A');
% fprintf(fileID,'%6.2f\n',NetAssimilation);
% fclose(fileID);
