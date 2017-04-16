% Make goodness of fit Supplementary figure


num        = 1;
wid        = 17.6; % total width
hts        = [6];
cols       = {2};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 3.5, 2, [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time

for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        % Look at goodness of fit statistics per subject:
        hold on
        plot([0, 1], [0, 1], '--k')
        plot(simp_baseR2, baseR2, 'o')
        plot(simp_baseR2(simp_baseP<.05), baseR2(simp_baseP<.05), 'o', 'markerEdgeColor', 'k', 'markerFaceColor', 'b')
        ylabel('R-squared with nuisance')
        xlabel('R-squared base model')
        title('Baseline pupil goodness of fit')
        set(gca, 'box', 'off')
        
        
    elseif xx==2
        hold on
        plot([0, 1], [0, 1], '--k')
        plot(simp_respR2, respR2, 'o')
        plot(simp_respR2(simp_respP<.05), respR2(simp_respP<.05), 'o', 'markerEdgeColor', 'k', 'markerFaceColor', 'b')
        ylabel('R-squared with nuisance')
        xlabel('R-squared base model')
        title('Peak response goodness of fit')
        set(gca, 'box', 'off')
        
        
    elseif xx==3 

        hold on
        bar(nanmean(rel_baseline_AIC))
        plot([1, 1; 2, 2]', [nanmean(rel_baseline_AIC)+ nanstd(rel_baseline_AIC)./sqrt(length(rel_baseline_AIC)); ...
            nanmean(rel_baseline_AIC)- nanstd(rel_baseline_AIC)./sqrt(length(rel_baseline_AIC))], 'k', 'lineWidth', 6)
        ylim([-1, 1].*400)
        ylabel('realtive AIC ')
        set(gca, 'xtick', [1, 2], 'xticklabels', xTickLabs, 'box', 'off')
    
    elseif xx==4
        
        hold on
        bar(nanmean(rel_response_AIC))
        plot([1, 1; 2, 2]', [nanmean(rel_response_AIC)+ nanstd(rel_response_AIC)./sqrt(length(rel_response_AIC)); ...
            nanmean(rel_response_AIC)- nanstd(rel_response_AIC)./sqrt(length(rel_response_AIC))], 'k', 'lineWidth', 6)
        ylim([-1, 1].*850)
        ylabel('relative AIC')
        set(gca, 'xtick', [1, 2], 'xticklabels', xTickLabs)

    end
    
    setPLOT_panelLabel(gca, xx);
end


disp(sprintf('Mean baseline AIC:  \n \n null: %g  \n base: %g  \n full: %g', nanmean(null_base_AIC), nanmean(simp_base_AIC), nanmean(full_base_AIC)))
disp(sprintf('\n \n \nMean response AIC:  \n \n null: %g  \n base: %g  \n full: %g', nanmean(null_resp_AIC), nanmean(simp_resp_AIC), nanmean(full_resp_AIC)))


kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure S2', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
%plot2svg('berlinFig4.svg', gcf)
%saveas(gcf,  'berlinFig4.fig', 'fig')
%saveas(gcf,  'audPaperFigS2.eps', 'epsc2')
%close(gcf)