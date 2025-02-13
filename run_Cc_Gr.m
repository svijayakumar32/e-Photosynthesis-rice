%% 19/09/2024
%% Plot A/Cc curve for measured average parameters from IRRI data and using Ac and Aj equations from Vcmax_adj_simple and Jmax_adj_simple
% Vcmax, J and TPU parameters are the averaged values of parameters derived from fitting eight curves 
% using the R package msuRACiFit described in Gregory et al (2021) and is available at: https://github.com/poales/msuRACiFit

%% Farquhar model parameters
Vcmax_m=134.380530667391;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
J=189.091184437385;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
TPU=12.9252127787927;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
Lii=2000;% Light intensity from IRRI

Gr = 49.6558492521714;% CO2 compensation point (µmol mol-1)
Gr_Pa = 5.02248975; % CO2 compensation point (Pa)

Rd = 1.33898789673117; % Day respiration (umol m-2 s-1)
Gm = 0.545813281251705; % Mesophyll conductance (umol m-2 s-1)
Kc_air = 606.946980366782;% Michaelis Menten of Rubisco carboxylation in 21% O2 (umol mol-1)

O=210;%mbar

%%%%%%%%%%%%%%%%%%%%%

%% e-Photosynthesis model parameters
%% For running the untuned e-Photosynthesis model
global Vrubusco;
Vrubusco=1.0;%adjust enzyme activity to 1
global VmaxAdj;
VmaxAdj=1.0;%adjust enzyme activity to 1 

%% For running the tuned e-Photosynthesis model
global Vrubusco_opt;
Vrubusco_opt=1.32;%adjust enzyme activity sub in optimal Vrubisco 
global VmaxAdj_opt;
VmaxAdj_opt=1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj 

global pcfactor;  
ProteinTotalRatio=0;
%pcfactor=1/ProteinTotalRatio;
%21/05/24 Change pcfactor to 1 here to avoid Inf
pcfactor=1;

Einput=ones(37,1);
Edata=importdata('Einput7.txt');
%Eio=Edata.data(:,1);
Eio=Edata.data(:,1);
Eio_opt=Edata.data(:,1);

%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp=28.9310407291759; %Avg Tleaf, Original = 25C
GRNC=0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

%%%%%%%%%%%%%%%%%%%%%
% Calculate Farquhar A
% Work backwards to get Ci from Cc and Gm
% Cc= Ci - Net_Ac/Gm;
% Ci_Ac = Gr+(Net_Ac/Gm);
% Ci_Aj = Gr+(Net_Aj/Gm);

% Ac_b = Vcmax_m-Rd+(Ci(i)+Kc_air) * Gm;
% Ac_c = ((Ci(i)-Gr) * Vcmax_m-(Ci(i)+Kc_air) * Rd) * Gm;
% Aj_b = (J/4) - Rd + (Ci(i)+ 2 * Gr) * Gm;
% Aj_c = ((Ci(i) - Gr) * J/4 - (Ci(i) + 2 * Gr) * Rd) * Gm;

% Set enzyme activities for e_Photosynthesis - nonopt is both alpha values = 1, opt is aRubisco = 1.32, aEnzymes = 1.12 
    Eio(1) = Edata.data(1,1) * Vrubusco;                    %aRubisco = 1
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;               %aEnzymes = 1
    Eio_opt(1) = Edata.data(1,1) * Vrubusco_opt;            %aRubisco = 1.32
    Eio_opt(2:27) = Edata.data(2:27,1) * VmaxAdj_opt;       %aEnzymes = 1.12  

    %  Calculate Gross Assimilation for e-Photosynthesis model at GammaStar
    PPFDi = Lii; 
    %Cc = 51.68068;
    % Define the range and tolerance of Cc
    lowerBound = 51.68;
    upperBound = 51.69;
    tolerance = 1e-4; % Precision to test
    % Use a loop to find the Cc that gives an output of 0
    Cc = lowerBound:tolerance:upperBound; % Generate test values
    Gross_A = zeros(length(Cc),1);
    Gross_A_opt = zeros(length(Cc),1);

for i = 1:length(Cc)
    Gross_A(i) = EPS_Drive_GRNs(Einput,Cc(i),PPFDi,WeatherTemp,GRNC,0,Eio); %aRubisco = 1, aEnzymes = 1
    Gross_A_opt(i) = EPS_Drive_GRNs(Einput,Cc(i),PPFDi,WeatherTemp,GRNC,0,Eio_opt); %aRubisco = 1.32, aEnzymes = 1.12
    if abs(Gross_A) < tolerance
        fprintf('Value %f gives an output of 0 (or close to 0 within tolerance).\n', x(i));
        break;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%