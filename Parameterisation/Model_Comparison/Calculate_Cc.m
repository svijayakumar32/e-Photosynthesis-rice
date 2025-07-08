%% Calculating Cc values equivalent to Ca for FvCB model fit

%% Specify average photosynthetic parameters
% Farquhar model parameters
Vcmax_m=134.380530667391;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
J=189.091184437385;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
TPU=12.9252127787927;% Avg Estimated using msuRACiFit (Gregory et al 2021) with rice Gamma Star at 28.9C
Gr = 49.6558492521714;%µmol mol-1
Rd = 1.33898789673117; %umol m-2 s-1
Gm = 0.545813281251705; %umol m-2 s-1
Kc_air = 606.946980366782;%umol mol-1

% Create a vector of CA whose elements increment by 10 (starting from 10).
CA=[0:1:2000]';
% To calculate the assimilation rates, we need three equations associated with limitation phases

%% Assimilation rates for Farquhar model divided by limitation processes (Ac, Aj, Ap) 
Farq_Matrix=zeros(2001,8);
k = 1; % Initialize row index
    % Farquhar model A calculation 
    for i = 1:2001
        Air_CO2 = CA(i);
        Ci = Air_CO2 * 0.7; % intercellular CO2 
        %if (i <= 30) % Rubisco- limited assimilation (Ac)
        %if (CA(i) <= 300) % Rubisco- limited assimilation (Ac)
        Ac_b = Vcmax_m-Rd+(Ci+Kc_air)*Gm;
        Ac_c = ((Ci-Gr)*Vcmax_m-(Ci+Kc_air)*Rd)*Gm;
        Net_Ac = ((Ac_b - sqrt(Ac_b^2 - 4*Ac_c)) / 2); % Net A expressed as a function of Cc does not include Rd
        Cc_Ac = (Ci - Net_Ac/Gm); % Use Net_Ac to calculate Cc
        Gross_Ac = Net_Ac + Rd; 
        %elseif (i >= 31) && (i <= 100) % RuBP limited assimilation (Aj)
        %elseif (CA(i) >= 301) && (CA(i) <= 1140) % RuBP limited assimilation (Aj)
        Aj_b = (J/4) - Rd + (Ci + 2*Gr) * Gm;
        Aj_c = ((Ci - Gr) * J/4 - (Ci + 2*Gr) * Rd) * Gm;
        Net_Aj = ((Aj_b - sqrt(Aj_b^2 - 4*Aj_c)) / 2); % Net A expressed as a function of Cc does not include Rd
        Cc_Aj = (Ci - Net_Aj/Gm); % Use Net_Aj to calculate Cc
        Gross_Aj = Net_Aj + Rd; 
        %else % TPU limited assimilation (Ap)
        Ap_b = (3*TPU)-Rd + (Ci- (1+(3*0))* Gr) *Gm;
        Ap_c = ((Ci-Gr)*3*TPU - (Ci-(1+(3*0))*Gr)*Rd)*Gm;
        Net_Ap = ((Ap_b - sqrt(Ap_b^2 - 4*Ap_c)) / 2); % Net A expressed as a function of Cc does not include Rd
        Cc_Ap = (Ci - Net_Ap/Gm); % Use Net_Ac to calculate Cc
        Gross_Ap = Net_Ap + Rd;
        % Compile results in matrix - Include all assimilation rates and Cc
        % values corresponding to CA and Ci
        Farq_Matrix(k, 1) = CA(i);
        Farq_Matrix(k, 2) = Ci;
        Farq_Matrix(k, 3) = Gross_Ac;
        Farq_Matrix(k, 4) = Gross_Aj;
        Farq_Matrix(k, 5) = Gross_Ap;
        Farq_Matrix(k, 6) = Cc_Ac; 
        Farq_Matrix(k, 7) = Cc_Aj; 
        Farq_Matrix(k, 8) = Cc_Ap; 
        k = k + 1; % Move to the next row
    end

%% Plotting Cc vs A 
figure
plot(Farq_Matrix(:,6),Farq_Matrix(:,3),'Color',"#A2142F",'LineWidth',1.5) %plot Cc_Ac against Gross_Ac in red
hold on
plot(Farq_Matrix(:,7),Farq_Matrix(:,4),'Color',"#0072BD",'LineWidth',1.5) %plot Cc_Aj against Gross_Aj in blue
hold on
plot(Farq_Matrix(:,8),Farq_Matrix(:,5),'Color',"#EDB120",'LineWidth',1.5) %plot Cc_Ap against Gross_Ap in yellow
legend('\it A_{c}','\it A_{j}','\it A_{p}')
xlabel('C_c (μmol mol^{-1})')
ylabel('Gross Assimilation Rate (μmol m^{-2} s^{-1})')
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Limitations"),'-dpdf');
% %% Plotting Cc vs A (combined graph)
% % Find intersection points to get ranges to plot for each limitation
% figure
% if Ac(i) < Aj(i)
%     p=plot(Cc(1:10),Ac(1:10),'-o');
%     %p=plot(Ci,Ac,'-o');
%     p.Color = "#D95319";
%     hold on
%     elseif Aj(i)<Ap(i)
%     q=plot(Cc(10:20),Aj(10:20),'-o');
%     %q=plot(Ci,Aj(),'-o');
%     q.Color = "#0072BD";
%     hold on
%     else
%     r=plot(Cc(20:33),Ap(20:33),'-o');
%     %r=plot(Ci,At,'-o');
%     r.Color = "#EDB120";
% end
% xlabel('C_c (ppm)');
% ylabel('A (umol m^-2 s^-1)');
% legend('Rubisco','RuBP regeneration','TPU','Location','NorthWest')
% set(gcf, 'PaperOrientation', 'landscape');
% print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Limitations_combined"),'-dpdf');
%% Plotting Ci vs A 
figure
plot(Farq_Matrix(:,2),Farq_Matrix(:,3),'Color',"#A2142F",'LineWidth',1.5) %plot Ci_Ac against Gross_Ac in red
hold on
plot(Farq_Matrix(:,2),Farq_Matrix(:,4),'Color',"#0072BD",'LineWidth',1.5) %plot Ci_Aj against Gross_Aj in blue
hold on
plot(Farq_Matrix(:,2),Farq_Matrix(:,5),'Color',"#EDB120",'LineWidth',1.5) %plot Ci_Ap against Gross_Ap in yellow
legend('\it A_{c}','\it A_{j}','\it A_{p}')
xlabel('C_i (μmol mol^{-1})')
ylabel('Gross Assimilation Rate (μmol m^{-2} s^{-1})')
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Ci_vs_Limitations"),'-dpdf');
%% Plotting CA vs A (Check intersection of phases)
figure
plot(Farq_Matrix(:,1),Farq_Matrix(:,3),'Color',"#A2142F",'LineWidth',1.5)
hold on
plot(Farq_Matrix(:,1),Farq_Matrix(:,4),'Color',"#0072BD",'LineWidth',1.5)
hold on
plot(Farq_Matrix(:,1),Farq_Matrix(:,5),'Color',"#EDB120",'LineWidth',1.5)
legend('\it A_{c}','\it A_{j}','\it A_{p}')
xlabel('C_a (μmol mol^{-1})')
ylabel('Gross Assimilation Rate (μmol m^{-2} s^{-1})')
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Ca_vs_Limitations"),'-dpdf');