% 21/06/2024
% Plot A/Cc curve for measured average parameters from IRRI data and using Ac and Aj equations from Vcmax_adj_simple and Jmax_adj_simple
% Vcmax, J and TPU parameters are the averaged values of parameters derived from fitting eight curves 
% using the R package msuRACiFit described in Gregory et al (2021) and is available at: https://github.com/poales/msuRACiFit

Vcmax_m=134.380530667391;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
J=189.091184437385;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
TPU=12.9252127787927;% Mean of 8 curves estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
Lii=2000;% Light intensity from IRRI

%Farquhar model parameters
Gr = 49.6558492521714;% CO2 compensation point (µmol mol-1)
Rd = 1.33898789673117; % Day respiration (umol m-2 s-1)
Gm = 0.545813281251705; % Mesophyll conductance (umol m-2 s-1)
Kc_air = 606.946980366782;% Michaelis Menten of Rubisco carboxylation in 21% O2 (umol mol-1)

O=210;%mbar

% Specify range of CAs 
%Ca=(100:50:2000)';
Ca=(100:50:1200)';

% Preallocate vectors to store outputs including Ci calculated from CA, and
% constants for Rubisco limited assimilation, RuBP regeneration limited
% assimilation and TPU limited assimilation
[Ci, ...
    Ac_b, Ac_c, Aj_b, ...
    Aj_c, Ap_b, Ap_c, ...
    Net_Ac, Net_Aj, Net_Ap, ...
    Cc_Ac, Cc_Aj, Cc_Ap, ...
    Gross_Ac, Gross_Aj, Gross_Ap] = deal(zeros(numel(Ca),1));

for i = 1:numel(Ca)
% Calculate Ci from CA and various assimilation rates
    Ci(i)= 0.7*Ca(i);
    Ac_b(i) = Vcmax_m-Rd+(Ci(i)+Kc_air) * Gm;
    Ac_c(i) = ((Ci(i)-Gr) * Vcmax_m-(Ci(i)+Kc_air) * Rd) * Gm;
    Aj_b(i) = (J/4) - Rd + (Ci(i)+ 2 * Gr) * Gm;
    Aj_c(i) = ((Ci(i) - Gr) * J/4 - (Ci(i) + 2 * Gr) * Rd) * Gm;
% Calculate net assimilation rates under Rubisco limitation
    Net_Ac(i) = ((Ac_b(i) - sqrt(Ac_b(i)^2 - 4 * Ac_c(i))) / 2);
    Cc_Ac(i) = (Ci(i) - Net_Ac(i)/Gm); % Use Net_Ac to calculate Cc
    Gross_Ac(i) = Net_Ac(i) + Rd; 
% Calculate net assimilation rates under RuBP limitation
    Net_Aj(i) = ((Aj_b(i) - sqrt(Aj_b(i)^2 - 4 * Aj_c(i))) / 2);
    Cc_Aj(i) = (Ci(i) - Net_Aj(i)/Gm); % Use Net_Aj to calculate Cc
    Gross_Aj(i) = Net_Aj(i) + Rd; 
% Calculate net assimilation rates under TPU limitation
    Net_Ap(i) = 3 * TPU-Rd;
    Cc_Ap(i) = (Ci(i) - Net_Ap(i)/Gm); % Use Net_Ap to calculate Cc
    Gross_Ap(i) = Net_Ap(i) + Rd;
end

%% Plotting Cc vs A (combined graph)
% Use intersection points to get ranges to plot net A for each limitation
figure
for j = 1:numel(Net_Ac)
    if     Net_Ac(j) < Net_Aj(j)
           p = scatter(Cc_Ac(j),Net_Ac(j),'filled',MarkerFaceColor="#D95319");
           %p.Color = "#D95319";
    elseif Net_Aj(j) < Net_Ap(j)
           hold on
           q = scatter(Cc_Aj(j),Net_Aj(j),'filled',MarkerFaceColor="#0072BD");
           %q.Color = "#0072BD";
    else   
           hold on
           r=scatter(Cc_Ap(j),Net_Ap(j),'filled',MarkerFaceColor="#EDB120");
           %r.Color = "#EDB120";
    end
    hold on
end
xlabel('C_c (μmol mol^{-1})');
ylabel('Net Assimilation Rate (umol m^{-2} s^{-1})');
%lgd = legend([p,q,r],'Rubisco','RuBP regeneration','TPU','Location','southoutside');

% Insert textbox listing parameter values
dim = [0.62 0.5 0.12 0.2];
str = {'\it{V_{cmax}}\rm{} = 134.38 ± 5.40','\it{J}\rm{} = 189.09 ± 5.54','\it{TPU}\rm{} = 12.93 ± 0.32','\it{R_d}\rm{} = 1.34 ± 0.28','\it{g_m}\rm{} = 5.40 ± 0.28'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','left');

% Import data for curves
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Limitations_combined"),'-dpdf');

% Insert additional textbox containing Vcmax and Jmax into plot 
% dim = [0.68, 0.22, 0.1, 0.1];
% str = {'V_{cmax} = 105.15','J_{max}   = 158.69','TPU   = 10.13'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');