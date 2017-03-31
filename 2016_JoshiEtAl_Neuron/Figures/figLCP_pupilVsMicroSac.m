function figLCP_pupilVsMicroSac(num, session)
% function figLCP_pupilVsMicroSac(num, session)
%  siteData{1}: trialsxcols matrix, cols are:
%   1 ... fix start time wrt fixation on
%   2 ... fix end time wrt fix start time (fix duration)
%   3 ... reported correcdt
%   4 ... beep on time, wrt to fix start time
%   5 ... trial begin time, wrt to fix start time
%   6 ... trial end time, wrt to fix start time
%   7 ... trial wrt time (cpu clock)
%   8 ... LFP index corresponding to fix start time (coded above)
%
%  siteData{2}: Analog:
%   dim1: trial
%   dim2: sample
%   dim3: 1=x, 2=y, 3=z-pupil, 4=corrected z-pupil, 5=slope
%   [remember first sample is considered time=0 (i.e., wrt fix start time)]
%     = eyedat(Lgood,:,:);
%
%  siteData{3}: spikes, re-coded wrt fix start time
%       first is multi-unit
%       rest are single-unit
%
%  siteData{4}: LFP, re-coded wrt fix start time
%
% SiteData{5}: pupil events... rows are events, columns are:
%   1. trial number
%   2. start time of event (wrt fix start time)
%   3. end time of event (wrt fix start time)
%   4. magnitude at start of event (raw z-score)
%   5. magnitude at end of event (raw z-score)
%   6. magnitude at start of event (corrected z-score)
%   7. magnitude at end of event (corrected z-score)
%   8. time of subsequent max slope
%   9. magnitude of subsequent max slope (corrected z/sample)


%% Collect event stats: duration, baseline, size, vmax
monkeys = {'Oz' 'Cicero'};
sites   = {'LC' 'IC' 'SC'};
for mm = 1:2
    for ss = 1:length(sites)
        [base_dir, fnames] = getLCP_cleanDataDir(monkeys{mm}, sites{ss});
        for ff = 1:length(fnames)
            disp(sprintf('%s: %s, file %d/%d', ...
                monkeys{mm}, sites{mm}, ff, length(fnames)))
            load(fullfile(base_dir, fnames{ff}));
            Lnb = ~isfinite(siteData{1}(:,4));
            Fnb = find(Lnb);
            for tt = 1:length(Fnb)
                Ltt = isfinite(siteData{2}(Fnb(tt),:,4));
                tax = (1:sum(Ltt));
                Xin = siteData{2}(Fnb(tt),Ltt,1)';
                Yin = siteData{2}(Fnb(tt),Ltt,2)';
                Lg  = isfinite(Xin) & isfinite(Yin);
                if any(~Lg)
                    Xin(~Lg) = linterp(tax(Lg), Xin(Lg), tax(~Lg));
                    Yin(~Lg) = linterp(tax(Lg), Yin(Lg), tax(~Lg));
                end
                
                 
                sacs = findMicroSaccades(Xin, Yin, true);
                r = input('next');
            end
        end
        
                %  sac(1:num,1)   onset of saccade
                %  sac(1:num,2)   end of saccade
                %  sac(1:num,3)   peak velocity of saccade (vpeak)
                %  sac(1:num,4)   horizontal component     (dx)
                %  sac(1:num,5)   vertical component       (dy)
                %  sac(1:num,6)   horizontal amplitude     (dX)
                %  sac(1:num,7)   vertical amplitude       (dY)
                [sac, radius] = microsacc_plugin(...
                    XinC, YinC,[0 0; diff(XinC) diff(YinC)], ...
                    0.01,5,1,1);
                
                subplot(2,1,1); cla reset; hold on
                plot(Xin)
                plot(Yin,'r')
                plot(XinC, 'c')
                plot(YinC,'m')
                for xx = 1:size(sac,1)
                    plot(sac(xx,1).*[1 1], [-1 1], 'k-');
                end
                
                
                
                subplot(2,1,2); cla reset; hold on                
                vel  = [0; sqrt(diff(XinC).^2 + diff(YinC).^2)*1000];
                vels = vel - median(vel);
                plot(vels);
                r = input('next');
            end
                
                fc = 25;
                fs = 1000;
                H0 = 2*fc/fs;
                M  = 28;
                Ms = -M:M;
                HM = nans(length(Ms),1);
                for mm = 0:M
                    HM(abs(Ms)==mm) = sin(H0*pi*mm)./(pi*mm);
                end
                HM(Ms==0) = H0;
                
                [sac, radius] = microsacc_plugin( ...
                    [Xin, Yin],vecvel([Xin, Yin],1000),VFAC,MINDUR,MSDX,MSDY)
                sacs = findMicroSaccades(Xin, Yin, ...
                    1000, 10, true);
                plot(siteData{2}(Fnb(tt),Ltt,4),'g', 'LineWidth', 2)
                r = input('next');
            end
                
                cla reset; hold on;
                Ltt = isfinite(siteData{2}(Fnb(tt),:,4));
                plot(siteData{2}(Fnb(tt),Ltt,1),'r')
                plot(siteData{2}(Fnb(tt),Ltt,2),'g')
                plot(siteData{2}(Fnb(tt),Ltt,3),'b')
                plot(siteData{2}(Fnb(tt),Ltt,4),'k')
                r = input('next');
            end
            %             rdat = nans(nnb,2,2); % R,P per trial for x/y eye vs PD
            %             for tt = 1:nnb
            %                 for rr = 1:2
            %                     Lg = isfinite(siteData{2}(Fnb(tt),:,rr)) & ...
            %                         isfinite(siteData{2}(Fnb(tt),:,4));
            %                     [R,P] = corr(siteData{2}(Fnb(tt),Lg,rr)', ...
            %                         siteData{2}(Fnb(tt),Lg,4)', 'type', 'Spearman');
            %                     rdat(tt,:,rr) = [R P];
            %                 end
            %             end
            %             vals{mm,2} = cat(1, vals{mm,2}, rdat);
        end
    end
end

%% Third and fourth panels are histograms of PD "event" durations
xs = 0:20:2000;
for mm = 1:2
    subplot(4,2,4+mm); cla reset; hold on;    
    Nd = hist(vals{mm}(vals{mm}(:,3)>0,1), xs);
    H(1) = bar(xs,Nd);
    Nc = hist(vals{mm}(vals{mm}(:,3)<0,1), xs);
    H(2) = bar(xs,-Nc);
    set(H,'FaceColor','k') 
    xlim([0 1500])
    xlabel('Duration (ms)')
    ylabel('Count')
    disp(sprintf('Monkey %s: median d=%.2f, median c=%.2f', ...
        monkeys{mm}, nanmean(vals{mm}(vals{mm}(:,3)>0,1)), ...
        nanmean(vals{mm}(vals{mm}(:,3)<0,1))))
end


%% Fifth and sixth panels are density maps of event size vs baseline
for mm = 1:2
    subplot(4,2,6+mm); cla reset; hold on;
    hist3(vals{mm}(:,2:3),[200 200]);
    set(gca, 'View', [0 90])
    set(gcf,'renderer','opengl');
    set(get(gca,'child'), ...
        'LineStyle', 'none', 'FaceColor','interp','CDataMode','auto');
    axis([-2.5 2.5 -2.5 2.5])
    xlabel('Baseline PD (z-score)');
    ylabel('Phasic PD (z-score)');
end
             
%% JUNK BELOW
%  density maps of vmax vs event size
% for mm = 1:2
%     subplot(4,2,6+mm); cla reset; hold on;
%     hist3(vals{mm}(:,2:3),[200 200]);
%     set(gca, 'View', [0 90])
%     set(gcf,'renderer','opengl');
%     set(get(gca,'child'), ...
%         'LineStyle', 'none', 'FaceColor','interp','CDataMode','auto');
%     %axis([-2.5 2.5 -2.5 2.5])
% end
% 
%             % 2: periodogram per trial, session
%             wrt0  = siteData{1}(1,7);
%             adat = nans(nnb, size(siteData{2},2)+1);
%             for tt = 1:nnb
%                 Ltt = isfinite(siteData{2}(Fnb(tt),:,4));
%                 [A,LAGS] = xcorr(siteData{2}(Fnb(tt),Ltt,4), 'coeff');
%                 adat(tt,1:sum(Ltt)) = A(LAGS>=0);
%             end
%             
%                 subplot(2,1,1); cla reset; hold on;
%                 plot(siteData{2}(Fnb(tt),Ltt,4));
%                 subplot(2,1,2); cla reset; hold on;
%                 plot(XX);
%                 r = input('next')
%             end
%                 
%             
%             tax   = (0:size(siteData{2},2)-1)';
%             ptdat = [];
%             for tt = find(Lnb)'
%                 Ltt = isfinite(siteData{2}(tt,:,4));
%                 ptdat = cat(1, ptdat, ...
%                     [((siteData{1}(tt,7)-wrt0)+tax(Ltt))./1000, ...
%                     siteData{2}(tt,Ltt,4)']);
%             end
%             %[Pf,Ff,ALPHAf] = lomb(ptdat(:,2), ptdat(:,1));
%             [P,F,ALPHA] = fastlomb(ptdat(:,2), ptdat(:,1));
%             vals{mm,2} = cat(1, vals{mm,2}, interp1(F,P,xi,'spline'));
%         end
%     end
% end
