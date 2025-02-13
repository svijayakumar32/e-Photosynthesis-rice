function [O2_sol] = calculate_O2_sol(Temp,Press) % Input Temp in C and Press in kPa
% CALCULATE_O2_SOL Calculate molar solubility of oxygen at a given temperature and pressure

% Given constants for oxygen - from JPL Publication 19-5. Chemical Kinetics and Photochemical Data 
% for Use in Atmospheric Studies (Table 5-4, Page 1146/1610)

A = -161.6;   % Empirical constant A - Intercept term representing baseline value of the natural logarithm of H at infinite temperature
B = 8160;     % Empirical constant B - Inverse temperature term describing temperature dependence of H
C = 22.39;    % Empirical constant C - Logarithmic temperature term accounting for non-linear temperature effects on H

T = Temp + 273.15;                 % Temperature in K 

% Calculate ln(H)
ln_HO2 = A + (B / T) + C * log(T); % Units of the result are in M atm -1 - (same as mol/L·atm?)

% Calculate H by taking exponent of logged H
H_atm = exp(ln_HO2);

% Convert to mol/L bar−1 (1 atm = 1.01325 bar), 1 mol L atm/1.01325
H_bar = H_atm*(1/(Press/100));

O2_sol = H_bar;

end