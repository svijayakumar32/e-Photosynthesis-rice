clear all;
% Fitting WT IR64 data provided by IRRI (Gas exchange measurement WT plants)
% Tleaf average of 8 replicates = 28.93104073C
Vcmax_m=134.380530667391;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
J=189.091184437385;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C

%%%%%%%%%%
Lii=2000;%light intensity from IRRI

%Farquhar model parameters
%Gr=50.2248975;%µbar %Gamma Star calculated from fitting IR64 of IRRI using Hermida-Carrera rice parameters at Tleaf
Gr = 49.6558492521714;%µmol mol-1

%Rd=1;%default
Rd = 1.33898789673117; %umol m-2 s-1, estimated using msuRACiFit (Gregory et al 2021) rice parameters at Tleaf 

%Gm = 5.395926789339; %umol m-2 s-1 Pa-1
Gm = 0.545813281251705; %umol m-2 s-1

I2=Lii/2*0.85*(1-0.15);
Theta=0.7;

%Kc=348.318266403677;%ubar, 34.832 Pa, avg of Tleaf values, rice c and dHa
%Ko=275.324214275076;%mbar, 27.532 kPa, avg of Tleaf values, rice c and dHa 
%Kc_air = 613.901254795337;%ubar, 61.390 Pa avg of Tleaf values, rice c and dHa
Kc_air = 606.946980366782;%umol mol-1

O=210;%mbar

%J=(I2+Jmax_m-sqrt((I2+Jmax_m)^2-4*Theta*I2*Jmax_m))/(2*Theta);
%%%%%%%%%%%%%%%%%%%%%
CA=[800,1200,1500,1800];
global Vrubusco_adj;
%Vrubusco_adj=1.12; 
Vrubusco_adj=1.0; % keep same rubisco activity as original Einput
global VmaxAdj;%adjust enzyme activity
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
MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp=28.9310407291759; %Avg Tleaf, Original = 25C
WeatherRH=0.6;
WeatherWind=5;
Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Radiation_PAR=Lii/Convert*0.85*0.85;%%%% Converts light in PPFD
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1; %Runs C3 Metabolic, not Farquhar Model
Vcmax25=100;
Jmax25=200;
GRNC=0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

% % Store CA, Ci and NetA in each loop iteration or Ci, GrossA, NetA, Gs, LeafT, Transp if using Leaf
% if MetaOnly==1
%    MetaMatrix=zeros(100,3);
% else 
%    MetaMatrix=zeros(100,5);
% end

for j=1:25
    j
VmaxAdj=1.0+j*0.02;%adjust enzyme activity  %modify ratio 1.2 to 1.4 (or greater) if no fit %default 0.8
Eio(1)=Edata.data(1,1)*Vrubusco_adj;
Eio(2:27)=Edata.data(2:27,1)*VmaxAdj;

for i=1:4
Air_CO2=CA(i);
if MetaOnly==1 %run only Metabolic model
CO2i=Air_CO2*0.7; % intercellular CO2 
Ci=CO2i;
PPFDi=Lii; 
NetAssimilation=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eio);
else
% Adding measured Rd, Gr as inputs to Leaf
LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio,Rd,Gr);
[row_indices, ~] = find(LeafResult);
last_row = max(row_indices); %
Ci=LeafResult(last_row,1); % Ci from last non zero row index in LeafResult
NetAssimilation=LeafResult(last_row,2); % NetAssimilation from last non zero row index in LeafResult
Gs=LeafResult(last_row,3); % Gs from last non zero row index in LeafResult
LeafTemperature=LeafResult(last_row,4); % LeafTemperature from last non zero row index in LeafResult
Transpirationi=LeafResult(last_row,5); % Transpirationi from last non zero row index in LeafResult
end %remove end if MetaOnly loop removed

%Calculate measured A values at different Ci
%ACI_m=J*(Ci-Gr)/(4*Ci+8*Gr)-Rd;

%05/01/2024 modification from Gu et al (2010) - quadratic form
b=(J/4)-Rd+(Ci+2*Gr)*Gm;
c=((Ci-Gr)*J/4-(Ci+2*Gr)*Rd)*Gm;

ACI_m=((b-sqrt(b^2-4*c))/2)+Rd; %%Aj expressed as a function of Ci : ACI_m+Rd if GrossA calculated instead of NetA as a result of EPS_Drive_GRNs
ACi_evsm(i)=(ACI_m-NetAssimilation)^2;%the squares of the residuals
end
SSR(j,1)=VmaxAdj;
SSR(j,2)=sum(ACi_evsm);%the sum of the squares of the residuals
end
% fileID = fopen('LeafmetaOut.txt','w');
% fprintf(fileID,'%6s\n','A');
% fprintf(fileID,'%6.2f\n',NetAssimilation);
% fclose(fileID);