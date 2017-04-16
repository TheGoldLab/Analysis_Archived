% makeFig8 --- this is now supplementary figure 3

% Use the J. Neurosci figure template
fig_num = 'S3'; % figure number
total_width = 17; % total width
col_height = [7.5 6]; % height of each column
cols_per_row = {0.7, [0.35 0.35]}; % columns in each row
psh = 2.5; % panel separation height
psw = 2.5; % panel separation width


num        = 8;
wid        = 18; % total width
hts        = [5];
cols       = {3};
% [axs,fig_] = getPLOT_axes(num, wid, hts, cols, 1.5, 1.5, [], 'Krishnamurthy et al');
% set(axs,'Units','normalized');
[axs,fig_] = getPLOT_axes(fig_num, total_width, col_height, cols_per_row,...
    psh, psw, [], '');
set(axs,'Units','normalized');

% draw in each panel, one at a time
cpColInd=2;
expColInd=3;
noiseColInd=4;
ms=8;


for xx = 1:length(axs)
    axes(axs(xx)); hold on; cla(gca)
    
    if xx==1
       % panel A:  perceptual bias, sorted by pupil predicted bias
        
        hold on        
        plot([-100, 100], [0, 0], '-','Color',[0.5 0.5 0.5])
        plot([-100, 100].*.2, [-100, 100].*.2, '-','Color',[0.5 0.5 0.5])
        [a, ~, p R points_high]=snow_makeBinnedPlot(selPredErr(hiBin),selPercErr(hiBin), .2, 'prediction error', 'perceptual error', 10, [0, 0, 0] , true, false)
        [a, ~, p R points_low]=snow_makeBinnedPlot(selPredErr(~loBin & ~hiBin),selPercErr(~loBin & ~hiBin), .2, 'prediction error', 'perceptual error', 10, [.3 .3 .3], true, false)
        [a, ~, p R points_low]=snow_makeBinnedPlot(selPredErr(loBin),selPercErr(loBin), .2, 'Prediction error', 'Perceptual error', 10, [.6, .6, .6], true, false)
        ff=legend([points_high, points_low], 'High pupil-predicted bias', 'Low pupil-predicted bias')
        ylim([-15, 15])
        xlim([-60, 60])
        set(ff, 'box', 'off', 'location', 'northwest')
        set(gca, 'box', 'off')
        grid on
        
        
    elseif xx==2
        % Panel B: Coefficients for pupil terms in pupil predicted bias
        % model
        hold on
        plot([0, length(Bw)+1], [0, 0], '-','Color',[0.5 0.5 0.5])
        Scale=max(max(abs(BINTw((size(baseTerms, 2)+1):end,:))));
        labels1={'S', 'R', 'W', 'S', 'R', 'W'};
        for i = (size(baseTerms, 2)+1):length(Bw)
            plot([i,i], BINTw(i,:), '-', 'color', 'k', 'lineWidth', 2)
            plot([i], Bw(i), 'o', 'markerFaceColor', [.7, .7, .7], 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', ms)
            
        end
        
        if useResTimepoints
            xlim([size(baseTerms, 2)+.5, size(baseTerms, 2)+3.5])
            labels1={'within', 'between', 'int.'}
            
        else
            plot([7.5, 7.5], [-1, 1].*.1, '-k', 'lineWidth', 1)
            xlim([4, 11])
        end
        
        ylabel({'Effect of pupil response', 'on perceptual bias'})
        ylim([-Scale, Scale])
        set(gca, 'xtick', (size(baseTerms, 2)+1):length(Bw), 'xticklabels', labels1);
        set(gca, 'box', 'off')
        
        
        % rows: within, between, interaction
        % columns: beta, t-stat, pvalue:
        
        rows={'within', 'between', 'interaction'}
        for i = 1:3
        disp(sprintf('Coefficient statistics for: %s pupil term: \n beta : %s \n t stat : %s \n p-value : %s', rows{i}, num2str(coeffStats(i,1)), num2str(coeffStats(i,2)), num2str(coeffStats(i,3))))
        end
        
        
    elseif xx==3
        % Panel C: AIC advantage over base model for pupil models
        % Create stars for likelihood ratio test results:
        %        LR_isSig=[pValue1, pValue2, pValue3]<.05;
        %        LR_isReallySig=[pValue1, pValue2, pValue3]<.0001;;
        
        
        % Print statistics for inclusion in text:
        rows={'base', 'fixed', 'random'}
        for i = 1:3
            disp(sprintf('Likelihood ratio statistics for %s model: \n test statistic : %s \n criterion value : %s \n p-value : %s', rows{i}, num2str(LikeRatio_stats(i,1)), num2str(LikeRatio_stats(i,2)), num2str(LikeRatio_stats(i,3))))
        end
        
        
        Scale=5+max(abs(AIC1(2:2:6)-AIC1(1:2:5)));
        a=bar( (AIC1(2:2:6)-AIC1(1:2:5)).*-1 );
        
        set(a, 'FaceColor', [.7 .7 .7])
        ylim([-Scale, Scale])
        ylabel({'AIC improvement', 'to base model'});
        
        
        allMods=1:3;
        
        plot(allMods(LR_isReallySig), 5+max(abs(AIC1(2:2:6)-AIC1(1:2:5))), '*k', 'lineWidth', 1,  'markerSize', 10, 'markerEdgeColor', 'k')
        
        labels={'NE', 'FE', 'RE'}
        %labels={'base', 'fixed\newlinemodel\newlineeffects', 'random\newlinemodel\newlineeffects'};
        
        set(gca, 'xtick', 1:3, 'xticklabel', labels);        
        set(gca, 'box', 'off');
        
    end
    
    
    setPLOT_panelLabel(gca, xx,0.14,1.05);
end


kk=annotation('textbox');
set(kk, 'string', 'Krishnamurthy et al, figure S3', 'position', ...
    [0.75 0.85 0.15 0.05], 'EdgeColor', 'none','fontsize',14)


figName = '~/Dropbox/auditoryTaskManuscript/Figure_set_for_writeup/figure_S3.pdf';
saveas(gcf,  'audPaper_figS3.fig', 'fig')
saveas(gcf,  'audPaper_figS3.eps', 'epsc2')

print(gcf, figName,'-dpdf', '-fillpage', '-painters')

%close(gcf)