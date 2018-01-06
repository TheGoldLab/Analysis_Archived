function fig_ = figBIAS_wkByMonk(num)
% function fig_ = figBIAS_wkByMonk(num)
%

if nargin < 1 || isempty(num)
    num = 6;
end

% get monk names, n monks, session end
[monks,monkn,mse] = getBIAS_monks;

% units should be in inches, from wysifig
wid  = 6.0; % total width
ht   = 1;
cols = {2,2,2,2};
[axs,fig_] = getBIAS_axes(num, wid, ht, cols);
set(axs, 'Units', 'normalized');

%% get the fig4 data, psdata...
wkdat = FS_loadProjectFile('2008_Bias', 'figBIAS_wkFits');

%%FOUR ROWS -- summary by monkey
%
% plot avg raw wks
%%%
wki   = [2 3];
tmax  = 20;
wkst  = 0.01;
wkmv  = 0.2;
wkgs  = (-wkmv:wkst:wkmv)';
wkgn  = length(wkgs);
wkg   = zeros(wkgn, tmax);
wkn   = length(wki);
%nsb   = 3;
%co = {'r' 'g' 'b'};

% for each monk...
for mm = 1:monkn

    % for each wk...
    for ww = 1:wkn

        % wkdat = wk x vins x fits x session
        % plot wk density, average
        axes(axs((mm-1)*2+ww)); cla reset; hold on;
        wkr    = permute(wkdat{mm,1}(1:tmax,wki(ww),1,1:mse(mm)),[4 1 2 3]);
        for ll = 1:tmax
            [N,X] = hist(wkr(:,ll), wkgs);
            wkg(:,ll) = N';
        end
        pcolor(1:tmax,wkgs,max(wkg(:))-wkg)
        shading('interp');
        colormap(hot); % bone
        set(gca,'XTick', [], 'YTick', []);
        
        plot([0 tmax], [0 0], 'k:');
        wkf = permute(wkdat{mm,1}(1:tmax,wki(ww),1,1:mse(mm)),[4 1 2 3]);
        plot(prctile(wkf,25), '-', 'Color', 'k', 'LineWidth', 1);
        plot(nanmedian(wkf),  '-', 'Color', 'k', 'LineWidth', 2);
        plot(prctile(wkf,75), '-', 'Color', 'k', 'LineWidth', 1);
        
        plot(nanmedian(permute(wkdat{mm,1}(1:tmax,wki(ww),2,1:mse(mm)),[4 1 2 3])), ...
            'g--', 'LineWidth', 2)
        %         sbins = round(linspace(1,mse(mm), nsb+1));
        %         for ss = 1:nsb
        %             plot(nanmedian(permute(wkdat{mm,1}(1:tmax,wki(ww),2,sbins(ss):sbins(ss+1)),[4 1 2 3])), ...
        %                 '--', 'Color', co{ss}, 'LineWidth', 2)
        %         end
        
        axis([1 tmax -.06 .06])
        
        if ww == 1
            ylabel('Filter coefficient');
        end
        if mm == 4
            xlabel('Trial lag');
        elseif mm == 1
            if ww == 1
                title('Correct trials');
            else
                title('Error trials');
            end
        end
    end
end


return



if 1
    % rs vs session
    axes(axs((mm-1)*3+3)); cla reset; hold on;
    xs = (1:size(wkdat{mm,5},1))';
    plot(xs([1 end]), [0 0], 'k:');
    Lup = psdat{mm}(:,1)<0.05;
    Llo = psdat{mm}(:,1)>0.95;

    % plot error bars
    es = [wkdat{mm,5}(:,2) wkdat{mm,5}(:,3)]';
    ys = wkdat{mm,5}(:,1);
    plot([xs xs]', es, 'k-');
    plot(xs,ys,'k.');

    plot([xs(Lup) xs(Lup)]', es(:,Lup), 'r-');
    plot(xs(Lup), ys(Lup), 'r.');
    
    plot([xs(Llo) xs(Llo)]', es(:,Llo), 'b-');
    plot(xs(Llo), ys(Llo), 'b.');
    set(gca, 'XTick', [50 100 150 200]);
    axis([1 xs(end) -0.5 0.5]);
    
    % weighted regression
    Lnn = isfinite(ys);
    A   = [ones(sum(Lnn),1) xs(Lnn)];
    [b,bint,h,p] = regressW(ys(Lnn), diff(es(:,Lnn))', A);
    plot(xs(Lnn), A*b, 'k--', 'LineWidth', 2);
    disp([b(2) p(2) sum(Lup) sum(isfinite(psdat{mm}(:,1)))])
    ylabel('r');
    if mm == 4
        xlabel('Session');
    end
end

return

%%%%
%% ROW 4: corr vs err
%%%%
wk1 = 1;
wk2 = 500;
amm = [-.2 .2];
for mm = 1:monkn

    nses = size(figdat{mm,1}, 4);
    
    wdat = nans(nses,4);
    for ss = 1:nses
        wdat(ss,:) = [ ...
            sum(figdat{mm,1}(1:wk1,2,2,ss)), sum(figdat{mm,1}(1:wk1,3,2,ss)), ...            
            sum(figdat{mm,1}(1:wk2,2,2,ss)), sum(figdat{mm,1}(1:wk2,3,2,ss))];
    end

    % wk summary
    axes(axs(mm+8)); cla reset; hold on;
%    plot(wdat(:,3), wdat(:,4), 'r.');
    plot(wdat(:,1), wdat(:,2), 'k.', 'MarkerSize', 9, 'Color', 'k');
    plot(amm, [0 0], 'k-');
    plot([0 0], amm, 'k-');
    plot(amm, amm, 'k:');
    axis([amm amm]);
end
axes(axs(9));
ylabel('Correct');
xlabel('Error');
