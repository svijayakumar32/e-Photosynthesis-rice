% 15/01/2024
% Calculating assimilation rates for FvCB model to check transition 
% between Rubisco-limited and electron transport limited phases in A/Ci fit

% Since we averaged Vcmax and J for A/Ci fits,
% It is important to check that an appropriate range of CA are used to calculate A : 
% Currently, CA = [100,150,200,250] for Ac and [800,1200,1500,1800] for Aj
% If using CA in ppm, convert Kcair/Gr/gm to umol/mol and umol/m2s units
% If converting CA to Pa, convert Kcair/Gr to Pa and gm in umol/m2sPa

% FvCB model parameters
Vcmax_m = 134.380530667391;% Avg Estimated using msuRACiFit (Gregory et al 2016) with rice Gamma Star at 28.9C
J = 189.091184437385;% Avg Estimated using msuRACiFit (Gregory et al 2016) with rice Gamma Star at 28.9C

%Gr = 50.2248975;%µbar %Gamma Star calculated from fitting IR64 of IRRI using Hermida-Carrera rice parameters at Tleaf
Gr = 49.6558492521714;%µmol mol-1

Rd = 1.33898789673117; %Rd estimated using msuRACiFit (Gregory et al 2016) rice parameters at Tleaf 

%Gm = 5.395926789339; %umol m-2 s-1 Pa-1
Gm = 0.545813281251705; %umol mol-1

%Kc = 348.318266403677;%ubar, 34.832 Pa, avg of Tleaf values, rice c and dHa
%Ko = 275.324214275076;%mbar, 27.532 kPa, avg of Tleaf values, rice c and dHa 
%Kc_air = 613.901254795337;%ubar, 61.390 Pa avg of Tleaf values, rice c and dHa

Kc_air = 606.946980366782;%umol mol-1

CA=[100,150,200,250,300]; % range of CA for Ac in ppm
%CA_j=[800,1200,1500,1800];% range of CA for Aj

global Vrubusco_adj;
global VmaxAdj;
VmaxAdj=1.12;%adjust enzyme activity i.e. sub in optimal Vmaxadj from Jmax_adj % default 1.3
Vrubusco_adj=1.0;%adjust enzyme activity between 0.6-1.0 %default 0.8
[Ci_vals,ACI_Aj,ACI_Ac,Ac_b,Ac_c,Aj_b,Aj_c,ACI_diffs]=(deal(zeros(4,1)));

for i=1:5
Air_CO2=CA(i);
%Air_CO2_j=CA_j(i); 

%Calculate measured A values at different Ci
CO2i = 0.7 * Air_CO2;
Ci_vals(i)=CO2i;
Ac_b(i)=Vcmax_m-Rd+(CO2i+Kc_air)*Gm;
Ac_c(i)=((CO2i-Gr)*Vcmax_m-(CO2i+Kc_air)*Rd)*Gm;

Aj_b(i)=(J/4)-Rd+(CO2i+2*Gr)*Gm;
Aj_c(i)=((CO2i-Gr)*J/4-(CO2i+2*Gr)*Rd)*Gm;

ACI_Ac(i)=(Ac_b(i)-sqrt(Ac_b(i)^2-4*Ac_c(i)))/2; %Ac expressed as a function of Ci, where
ACI_Aj(i)=(Aj_b(i)-sqrt(Aj_b(i)^2-4*Aj_c(i)))/2; %Aj expressed as a function of Ci, where
ACI_diffs(i)=ACI_Aj(i)-ACI_Ac(i);
end

% plot(Ci_vals,ACI_Ac)
% ylabel("A (µmol m^-2 s^-1)")
% xlabel("C_i (ppm)")