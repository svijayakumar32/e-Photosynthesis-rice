function [PsV1,PrV111,Ru_Act] = find_Rubisco_params(Vcmax,Tempin)
% FIND_RUBISCO_PARAMS A function to find the value of the Rubisco rate parameter PsV1 at different temperatures
% The calculations below were adapted from EPS_Drive_GRNs and CM_Rate to give the necessary outputs
% Q10 is a temperature coefficient - measuring the rate of change when increasing temperature by 10 °C
% ePhotosynthesis uses a Q10 function to model the temperature-dependent response for each enzyme reaction rate in the CBB cycle from their rates at 25 °C
% V1 is the maximum rate of RuBisCO-limited carboxylation (which is RuBP-saturated) in umol m-2 s-1 or mmol L-1 s-1
% V111 is the maximum rate of RuBisCO-limited oxygenation in umol m-2 s-1 or mmol L-1 s-1

Tp = Tempin; % Temperature input
PsV1_0 = Vcmax/30; % Unit conversion of initial Rubisco Vmax outlined in EPS_Drive_GRNs
Q10_1 = 1.93; % Q10 associated with Rubisco carboxylation
Ru_Act = -3E-05*Tp^3 + 0.0013*Tp^2 - 0.0106*Tp + 0.8839; % Rubisco activation state corrected for temperature
PsV1 = PsV1_0*Ru_Act*Q10_1^((Tp-25)/10); % Ps_V1_0 after Q10 adjustment
PrV111 = PsV1 * 0.24; % PrV111 is calculated as a proportion of PsV1
% PrPs_ratio = calculate_PsPr_ratio(Tp);
% PrV111= PsV1* PrPs_ratio;

% N.B. Currently the e-Photosynthesis model sets PrV111 = PsV1 * 0.24 (Whitney, 1999)
% However, this ratio may vary across different crop species in reality.
% To find out what the numbers are for tobacco, use the equations by Bernacchi (2001) and for rice, Makino (1988)
% to check the temperature response of the parameters - see Vmax_temp_adj and calculate_PrPs_ratio

