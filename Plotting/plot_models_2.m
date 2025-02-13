%% 18/09/2024
%% Simplified calibration for running only the metabolic model (MetaOnly==1)
%
Vcmax_m=134.380530667391;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
J=189.091184437385;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
Lii=2000;%Light intensity from IRRI

%Farquhar model parameters
Gr = 49.6558492521714;%µmol mol-1
Rd = 1.33898789673117; %umol m-2 s-1
Gm = 0.545813281251705; %umol m-2 s-1
Kc_air = 606.946980366782;%umol mol-1

O=210;%mbar

%%%%%%%%%%%%%%%%%%%%%
CA=[20:20:1000];% include full range of points (including intersection of GammaStar) extending through both Vcmax and J ranges

%For running the untuned e-Photosynthesis model
%global Vrubusco_nonopt;
Vrubusco_nonopt=1.0;%adjust enzyme activity to 1
%global VmaxAdj_nonopt;
VmaxAdj_nonopt=1.0;%adjust enzyme activity to 1 

% For running the tuned e-Photosynthesis model
%global Vrubusco_opt;
Vrubusco_opt=0.98;%adjust enzyme activity sub in optimal Vrubisco from Vcmax_adj_simple_3 
%global Vmax_opt;
VmaxAdj_opt=0.96;%adjust enzyme activity i.e. sub in optimal Vmaxadj from Jmax_adj_simple_3

global pcfactor;  
ProteinTotalRatio=0;
%pcfactor=1/ProteinTotalRatio;
%21/05/24 Change pcfactor to 1 here to avoid Inf
pcfactor=1;
%%%%%%%%%%%%%%%%%%%%%
Einput=ones(37,1);%No gene expression data input
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp=28.9310407291759; %Avg Tleaf, Original = 25C
GRNC=0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

%% Assimilation rates for Farquhar model 
% Other than 15 different Ci values, everything else is the same - 
% Maybe for indexing in later steps its easier to keep the repeated values of assimilation?
Farq_Matrix=zeros(numel(CA),6); %1 set of alpha values tested at 50 CA vals, number of cols depends on no. of output variables we want
k = 1; % Initialize row index
for j = 1:length(CA) %orig 25->50
    % Farquhar model A calculation 
    Air_CO2 = CA(j);
    Ci = Air_CO2 * 0.7; % intercellular CO2
    if j<=15
        b = Vcmax_m-Rd+(Ci+Kc_air)*Gm;
        c = ((Ci-Gr)*Vcmax_m-(Ci+Kc_air)*Rd)*Gm;
    else
        b = (J/4) - Rd + (Ci + 2*Gr) * Gm;
        c = ((Ci - Gr) * J/4 - (Ci + 2*Gr) * Rd) * Gm;
    %Previous Ac calculation, subbing in Kc_air for Kc*(1+O/Ko)
    %ACI_m=Vcmax_m*(Ci-Gr)/(Ci+Kc_air)-Rd; 
    % Net Ac expressed as a function of Cc does not include Rd
    Net_A = ((b - sqrt(b^2 - 4*c)) / 2); 
    % Gross_A is being calculated in EPS_Drive_GRNs instead of Net_A
    % Therefore Gross Ac expressed as a function of Cc should include Rd
    Gross_A = Net_A + Rd; 
    % Use Net_A to calculate Cc
    Cc = (Ci - Net_A/Gm); 
    %Check A when Cc = Gr;
    % Fill in results matrix
    Farq_Matrix(k, 1) = Vrubusco_opt; % just a sanity check / to align the columns with ePhoto
    Farq_Matrix(k, 2) = VmaxAdj_opt; % just a sanity check / to align the columns with ePhoto
    Farq_Matrix(k, 3) = CA(j);
    Farq_Matrix(k, 4) = Ci;
    Farq_Matrix(k, 5) = Cc; 
    Farq_Matrix(k, 6) = Gross_A;
    k = k + 1; % Move to the next row
    end
end

%% Assimilation rates for e-Photosynthesis model
ePhoto_Matrix=zeros(numel(CA),6);%orig 125->750
k = 1; % Initialize row index
for j = 1:length(CA)
    % Adjust enzyme activities
    Eio(1) = Edata.data(1,1) * Vrubusco_opt; % select alpha vals here
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj_opt;
    % e-Photosynthesis model A calculation 
    Air_CO2 = CA(j);
    Ci = Air_CO2 * 0.7; % intercellular CO2 
    Cc = Farq_Matrix(k,5); % check index
    %Check A when Cc = Gr;
    PPFDi = Lii; 
    GrossAssimilation = EPS_Drive_GRNs(Einput,Cc,PPFDi,WeatherTemp,GRNC,0,Eio);
    % Fill in results matrix
        ePhoto_Matrix(k, 1) = Vrubusco_opt; % should be same values all the way down
        ePhoto_Matrix(k, 2) = VmaxAdj_opt; % should be same values all the way down
        ePhoto_Matrix(k, 3) = CA(j);
        ePhoto_Matrix(k, 4) = Ci;
        ePhoto_Matrix(k, 5) = Cc;
        ePhoto_Matrix(k, 6) = GrossAssimilation; 
        k = k + 1; % Move to the next row
end

%If running on HEC, save result to check
save plot_models_2_result.mat;
