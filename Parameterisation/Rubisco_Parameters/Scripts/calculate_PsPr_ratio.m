function [PrPs_ratio] = calculate_PsPr_ratio(Temp)
% CALC_PRPS_RATIO Calculate Vo/Vc ratio for Rubisco based on temperature

T = Temp + 273.15;                 % Temperature in K 
R = 0.008314;                      % Ideal gas constant in kJ K-1 mol-1
c_Vc = 26.35;                      % Scaling constant, Vcmax             (N. tabacum L. cv. W38, Bernacchi et al 2001)
dHa_Vc = 65.33;                    % Activation energy, Vcmax, kJ mol -1 (N. tabacum L. cv. W38, Bernacchi et al 2001)
c_Vo = 22.98;                      % Scaling constant, Vomax             (N. tabacum L. cv. W38, Bernacchi et al 2001)
dHa_Vo = 60.11;                    % Activation energy, Vomax, kJ mol -1 (N. tabacum L. cv. W38, Bernacchi et al 2001)

Bernacchi_Vc_25 = exp(c_Vc - dHa_Vc / (R * T25)); 
Bernacchi_Vo_25 = exp(c_Vo - dHa_Vo / (R * T25));
Bernacchi_PrPs_ratio_25 = Bernacchi_Vo_25/Bernacchi_Vc_25;

Bernacchi_Vc = exp(c_Vc - dHa_Vc / (R * T)); 
Bernacchi_Vo = exp(c_Vo - dHa_Vo / (R * T));
Bernacchi_PrPs_ratio = Bernacchi_Vo/Bernacchi_Vc;

Makino_Vc_25 = 1.77;                            % RuBP carboxylase activity at pH 8.15 and 25°C from rice leaves (umol(mg enzyme)-1 min-1 (Vmax)
Makino_Vo_25 = 0.58;                            % RuBP oxygenase activity at pH 8.15 and 25°C from rice leaves (umol(mg enzyme)-1 min-1 (Vmax)
Makino_PrPs_ratio_25 = Makino_Vo_25/Makino_Vc_25;

Makino_PrPs_ratio = Makino_PrPs_ratio_25*Bernacchi_PrPs_ratio/Bernacchi_PrPs_ratio_25;
PrPs_ratio = Makino_PrPs_ratio;

end
