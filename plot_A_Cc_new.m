%% Plotting Cc vs A for FvCB and e-Photosynthesis models (combined graph)
% Use intersection points to get ranges to plot gross A for each limitation
load('Farq_ePhoto_results_new2.mat') % result of running Farq_ePhoto_comparison but with new alpha values
figure
for j = 1:numel(Net_Ac)
    if     Net_Ac(j) < Net_Aj(j)   
           p = scatter(Cc_Ac(j),Gross_Ac(j),20,'filled', MarkerFaceColor="#D95319");
           %p.Color = "#D95319";
    elseif Net_Aj(j) < Net_Ap(j)
           hold on
           q = scatter(Cc_Aj(j),Gross_Aj(j),20,'filled', MarkerFaceColor="#0072BD");
           %q.Color = "#0072BD";
    else   
           hold on
           r=scatter(Cc_Ap(j),Gross_Ap(j),20,'filled', MarkerFaceColor="#EDB120");
           %r.Color = "#EDB120";
    end
end

hold on

% Plot Cc (x) against Gross Assimilation (y) for alpha values=1 (black dashed lines)
nonopt_plot = plot(Cc_ePhoto, Gross_A_nonopt, ...
            '-x', ...
            'Color', 'k');
            %'MarkerFaceColor', 'k', ... % k=black markers and line
            %'MarkerEdgeColor', 'k', ...
            %'MarkerSize', 6, ...
            %'DisplayName', 'ePhotosynthesis (α_Enzymes = 1, α_Rubisco = 1)');

% Plot Cc (x) against Gross Assimilation (y) for aRubisco = 1.32, aEnzymes = 1.12  (red dashed lines)
opt_plot = plot(Cc_ePhoto, Gross_A_opt, ...
            '-x', ...
            'Color', 'r');
            % 'MarkerFaceColor', 'r', ... % k=black markers and line
            % 'MarkerEdgeColor', 'r', ...
            % 'MarkerSize', 6, ...
            % 'DisplayName', 'ePhotosynthesis (α_Enzymes = 1.12, α_Rubisco = 1.32)');

hold on 

xticks(0:100:800)
xlim([0 800]) % Plot from Cc = Gamma Star
ylim([0 40]);
xlabel('C_c (μmol mol^{-1})');
ylabel('Gross Assimilation Rate (μmol m^{-2} s^{-1})');
lgd = legend([p,q,r,nonopt_plot,opt_plot],...
    'Rubisco','RuBP regeneration','TPU',...
    'α_{Enzymes} = 1, α_{Rubisco} = 1','α_{Enzymes} = 1.12, α_{Rubisco} = 1.32',...
    'Location','southeast');
% Add legend, including FvCB
%legend('show', 'Location', 'east');

% Insert textbox listing parameter values
%dim = [0.61 0.5 0.12 0.2];
%str = {'\it{V_{cmax}}\rm{} = 134.38 ± 5.40','\it{J}\rm{} = 189.09 ± 5.54','\it{TPU}\rm{} = 12.93 ± 0.32','\it{R_d}\rm{} = 1.34 ± 0.28','\it{g_m}\rm{} = 5.40 ± 0.28'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','left');

hold off 

% Export plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Farq_ePhoto_new2"),'-dpdf','-bestfit');