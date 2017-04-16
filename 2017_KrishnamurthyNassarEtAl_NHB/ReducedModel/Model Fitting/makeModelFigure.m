% Make figure of model driven prior weight behavior



num        = 1;
wid        = 17; % total width
hts        = [5];
cols       = {3};
[axs,fig_] = getPLOT_axes(num, wid, hts, cols, 2, 2, [], '');
set(axs,'Units','normalized');
% draw in each panel, one at a time



highNoiseColor=cbColors(2,:);
lowNoiseColor=cbColors(3,:);


lw=2;
ms=8
pc=[.5 .5 .5];

for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    if xx==1
        highNoise=allNoise==20;
        makeTACfig(allStableTerm(highNoise), allTAC(highNoise)+1, 7, true, highNoiseColor, ms, -1,'k')
        makeTACfig(allStableTerm(~highNoise), allTAC(~highNoise)+1, 7, true, lowNoiseColor, ms, -1,'k')
        ylabel('Stable probability')
        xlabel('TAC')
        
        
        
        set(gca, 'box', 'off')
        setPLOT_panelLabel(gca, xx);
    elseif xx==2

        highNoise=allNoise==20;
        [a, p1]=makeTACfig(allUncPriorWeightTerm(highNoise), allTAC(highNoise)+1, 7, true, highNoiseColor, ms, -1,'k')
        [b, p2]=makeTACfig(allUncPriorWeightTerm(~highNoise), allTAC(~highNoise)+1, 7, true, lowNoiseColor, ms, -1,'k')
        ylabel('Uncertainty based prior weight')
        xlabel('TAC');
        ylim([0, .6]);
        set(gca, 'box', 'off')
        setPLOT_panelLabel(gca, xx);

    elseif xx==3
        priorWeightMod_SEM_hi=nanstd(condPriorWeightMod(:,:,2))./sqrt(size(medPercErrHi, 1))
        priorWeightMod_SEM_lo=nanstd(condPriorWeightMod(:,:,1))./sqrt(size(medPercErrLo, 1))
        
        errBars_hi=[nanmean(condPriorWeightMod(:,:,2))+priorWeightMod_SEM_hi; nanmean(condPriorWeightMod(:,:,2))-priorWeightMod_SEM_hi];
        errBars_lo=[nanmean(condPriorWeightMod(:,:,1))+priorWeightMod_SEM_lo; nanmean(condPriorWeightMod(:,:,1))-priorWeightMod_SEM_lo];
        
        hold on
        plot(barInds,errBars_hi, '-k', 'lineWidth', 1)
        plot(nanmean(condPriorWeightMod(:,:,2)), 'o', 'markerSize', ms, 'markerFaceColor', highNoiseColor, 'markerEdgeColor', 'k', 'lineWidth', 1)
        plot(barInds,errBars_lo, '-k', 'lineWidth', 1)
        plot(nanmean(condPriorWeightMod(:,:,1)), 'o', 'markerSize', ms, 'markerFaceColor', lowNoiseColor, 'markerEdgeColor', 'k', 'lineWidth', 1)
        ylabel('Prior weight')
        xlabel('Sounds after change-point')
        ylim(pw_lims)
        xlim([.5, 6.5])
        set(gca, 'box', 'off')
        setPLOT_panelLabel(gca, xx);
    end
end



kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al', 'position', [0.85 0.95 0.15 0.05], 'EdgeColor', 'none')
saveas(gcf,  'audTaskModFig.fig', 'fig')
saveas(gcf,  'audTaskModFig.eps', 'epsc2')
close(gcf)