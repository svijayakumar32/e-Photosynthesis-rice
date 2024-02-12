%% Simplified Vcmax_adj for running only the metabolic model (MetaOnly==1)
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
CA=[100,150,200,250,300];
global Vrubusco_adj;
global VmaxAdj;%adjust enzyme activity
VmaxAdj=1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj from Jmax_adj % default 1.3
global pcfactor;  
ProteinTotalRatio=0;
pcfactor=1/ProteinTotalRatio;
%%%%%%%%%%%%%%%%%%%%%

Einput=ones(37,1);%No gene expression data input
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
%MetaOnly=1;% if MetaOnly=1 run only Metabolic model
WeatherTemp=28.9310407291759; %Avg Tleaf, Original = 25C
GRNC=0; % In EPS_Drive_GRNs, if GRNC==0 cATPsyn,CPSi,cNADPHsyn and cpsii=1 are all set to 1

% Assimilation rates for e-Photosynthesis model
ePhoto_Matrix=zeros(125,4);
k = 1; % Initialize row index
for j = 1:25
    Vrubusco_adj = 1.0 + j * 0.02; % Adjust enzyme activity - start at 0.5, then 1.0 and check minima of SSR vals
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % e-Photosynthesis model A calculation 
    for i = 1:5
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        PPFDi = Lii; 
        GrossAssimilation = EPS_Drive_GRNs(Einput,Ci,PPFDi,WeatherTemp,GRNC,0,Eio);
        ePhoto_Matrix(k, 1) = VmaxAdj;
        ePhoto_Matrix(k, 2) = CA(i);
        ePhoto_Matrix(k, 3) = Ci;
        ePhoto_Matrix(k, 4) = GrossAssimilation;
        k = k + 1; % Move to the next row
    end
end

% Assimilation rates for Farquhar model 
% Other than 5 different Ci values, everything else is the same - 
% Maybe for indexing in later steps its easier to keep the repeated values of assimilation?
Farq_Matrix=zeros(125,4);
k = 1; % Initialize row index
for j = 1:25
    Vrubusco_adj = 1.0 + j * 0.02; % adjust enzyme activity
    Eio(1) = Edata.data(1,1) * Vrubusco_adj;
    Eio(2:27) = Edata.data(2:27,1) * VmaxAdj;
    % Farquhar model A calculation 
    for i = 1:5
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        b=Vcmax_m-Rd+(Ci+Kc_air)*Gm;
        c=((Ci-Gr)*Vcmax_m-(Ci+Kc_air)*Rd)*Gm;
        ACI_m = (b - sqrt(b^2 - 4*c)) / 2; % Ac expressed as a function of Ci 
        ACI_m = ACI_m + Rd; % If GrossA is being calculated in EPS_Drive_GRNs instead of NetA, ACI_m should be adjusted to match this 
        % Delete if this is not true
        Farq_Matrix(k, 1) = Vrubusco_adj;
        Farq_Matrix(k, 2) = CA(i);
        Farq_Matrix(k, 3) = Ci;
        Farq_Matrix(k, 4) = ACI_m;
        k = k + 1; % Move to the next row
    end
end

%Create vector to store Vmax_adj, differences in assimilation rates and corresponding SSRs
Diff_Matrix = zeros(125,2);

% Compute differences between two photosynthetic models
for k = 1:length(Farq_Matrix)
    for j = 1:25
    Diff_Matrix(k,1)=Farq_Matrix(k,1);
    Diff_Matrix(k,2)=(Farq_Matrix(k,4)-ePhoto_Matrix(k,4))^2;%the squares of the residuals
        while k<125
        k = k + 1; % Move to the next row
        end
    end
end

% Get sums of squared residuals by summing every four rows 
SSR_Matrix = zeros(25,2);
Vrubusco_adj_vals = Diff_Matrix(:,1);
SSR_Matrix(:,1) = Vrubusco_adj_vals(5:5:end);
SSR_Matrix(:,2) = sum(reshape(Diff_Matrix(:,2), 5, []))';
%toc