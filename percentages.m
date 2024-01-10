% ENZYME_PERTURBATION 
% This script generates random percentages of enzyme activities Vmax by drawing random floating point values between 0-1

rng("shuffle");
% Generate a random vector to get list of percentages to increase Vmax of Rubisco Aldolase SBPase

% Random percentages (drawn between 0 and 1 to act as multipliers)
a = 0; 
b = 1;
rand_multiplier = (b-a).*abs(rand(100,1)) + a;
rand_percent = rand_multiplier.*100;
%Check percentage is between limits
rand_range = [min(rand_percent) max(rand_percent)];

% % Aldolase percentages
% aldolase_low = 1; 
% aldolase_high = 100;
% aldolase = (aldolase_high-aldolase_low).*rand(100,1) + aldolase_low;
% a_range = [min(a) max(a)];
% 
% % SBPase percentages 
% sbpase_low = 1; 
% sbpase_high = 100;
% sbpase = (sbpase_low-sbpase_high).*rand(100,1) + sbpase_low;
% s_range = [min(s) max(s)];

% c1 = c(100); % select new percentages in vectors c1, f1 and p1 corresponding to nodes 1-100
% f1 = f(100);
% p1 = p(100);
% A = [c1, f1, p1];
% 
% alluptakes = [carbuptake; fatuptake; protuptake]; % recombine into single vector of uptakes