function getLCP_uStimVsSTA
% function getLCP_uStimVsSTA
%
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
%   dim3: 1 = x, 2 = y, 3 = pupil, 4 = corrected z-pupil, 5 = pupil slope
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

%% collect ustim filenames/data
% site, {monkey, type}, example_data (monkey type file)
ustim_sites = { ...
    'LC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'IC', {'Cicero', {'short'}; 'Oz', {'short' 'long'}}, [2 1 6]; ...
    'SC', {'Cicero', {'short'}; 'Oz', {'short'}},        [2 1 4]};
nus    = size(ustim_sites,1);
% load data from file
udat = FS_loadProjectFile('2013_LCPupil', 'uStim');
tmin = -1000;
tmax = 1000;
pax  = (tmin:tmax)';
npax = length(pax);
Lmp  = pax>0;
usdat = cell(nus, 2); % filename/[mm tt ff peak value]
for ss = 1:nus
    Lg = isfinite(udat{ss,2}(:,4)) & isfinite(udat{ss,2}(:,5));
    for mm = 1:size(ustim_sites{ss,2},1)
        Lm = Lg & udat{ss,2}(:,1)==mm;
        for tt = 1:size(ustim_sites{ss,2}{mm,2},2)
            Lmt = Lm & udat{ss,2}(:,2)==tt;
            % get directory
            cdir = fullfile(dirnames('local'), 'Data', 'Projects', ...
                '2013_LCPupil', 'Data', 'uStim', ...
                ustim_sites{ss,2}{mm,1}, ustim_sites{ss,1}, ustim_sites{ss,2}{mm,2}{tt}, 'clean');
            files = dir(fullfile(cdir, '*.mat'));
            for ff = 1:size(files,1)
                Lf = Lmt & udat{ss,2}(:,3)==ff;                
                if sum(Lf) > 4
                    % get PD at mean peak, plot index vs mean value
                    meanp = nanmean(udat{ss,1}(Lf,:,2));
                    mvs   = meanp(Lmp);
                    maxi  = find(mvs==max(mvs),1);
                    pis   = (-10:10)+maxi-tmin+1;
                    pis   = pis(pis>=1&pis<=npax);
                    fdat  = nans(sum(Lf), 1);
                    Ff    = find(Lf);
                    for rr = 1:sum(Lf)
                        fdat(rr) = nanmean(udat{ss,1}(Ff(rr),pis,2)-...
                            udat{ss,2}(Ff(rr),4));
                    end
                    usdat{ss,1} = cat(1, usdat{ss,1}, {files(ff).name});
                    usdat{ss,2} = cat(1, usdat{ss,2}, [mm tt ff maxi nanmean(fdat) nan nan nan]);
                end
            end
        end
    end
end


%% collect associated STA data
%% pick filenames by hand
fnames = {{ ...
    'Cicero_012214_15_5_FixnBeep.mat'; ...
    'Cicero_012314_d15_090_FixnBeep.mat'; ...
    'Cicero_022714_d18_469_FixnBeep.mat'; ...
    'Cicero_033114_d15_431_FixnBeep.mat'; ...
    'Cicero_041514_d16_120_FixnBeep.mat'; ...
    'Oz_060713_d_11_7_Fixn.mat'; ...
    ''; ...
    ''; ...
    'Oz_061813_d11_715_Fixn.mat'; ...
    ''; ...
    ''; ...
    ''}, ...
    { ...
    'Cicero_030314_d7_681_FixnBeep.mat'; ...
    'Cicero_031214_d10_705_FixnBeep.mat'; ...
    ''; ...
    ''; ...
    ''; ...
    'Cicero_032814_d10_099_FixnBeep.mat'; ...
    'Oz_060513_d8_868_Fixn.mat'; ...
    'Oz_061013_d8_167_Fixn.mat'; ...
    ''; ...
    'Oz_061213_d8_5_Fixn.mat'; ...
    'Oz_061213_d8_5_Fixn.mat'; ...
    'Oz_061313_d8_6_Fixn.mat'; ...
    'Oz_061713_d6_660_Fixn.mat'; ...
    'Oz_061713_d6_660_Fixn.mat'; ...
    'Oz_061313_d8_6_Fixn.mat'; ...
    'Oz_061713_d6_660_Fixn.mat'; ...
    'Oz_061913_d7_3.mat'; ...
    'Oz_062113_d11_216_FixnEnd.mat'}, ...
    { ...
    'Cicero_022414_d9_249_FixnBeep.mat'; ...
    ''; ...
    ''; ...
    'Cicero_032614_d7_148_FixnBeep.mat'; ...
    'Cicero_032614_d7_376_FixnBeep.mat'; ...
    ''; ...
    'Cicero_040114_d3_675_FixnBeep.mat'; ...
    'Cicero_040114_d3_675_FixnBeep.mat'; ...
    ''; ...
    ''}};

xax = (-1000:2000)';
Lx  = xax>0&xax<700; 
for ss = 1:nus    
    % get mean STA
    for ff = 1:size(fnames{ss},1)
        if ~isempty(fnames{ss}{ff})
            [base_dir, ~] = getLCP_cleanDataDir(...
                fnames{ss}{ff}(1:strfind(fnames{ss}{ff},'_0')-1), ustim_sites{ss,1});
            load(fullfile(base_dir, fnames{ss}{ff}));
            disp([ss ff])
            
            %% STA
            %  sdat is:pupil measure
            %   1: per unit (multi first, if appropriate)
            %   2: NUM_BINS
            %   3: 1-mean, 2-sem, 3-shuffled mean, 4=ps
%             [sdat,~] = getLCP_stpds(siteData, 'slope', 'multi', 0);
%             if size(sdat,1)>0
%                 usdat{ss,2}(ff,6:7) = [min(sdat(1,Lx,1)) max(sdat(1,Lx,1))];
%             end
%             
            %% Tonic correlation
            % Use only no-beep, no-stim trials
            Fnb = find(~isfinite(siteData{1}(:,4)) & ...
                ~isfinite(siteData{1}(:,9)));
            num_tr     = length(Fnb);
            tr_ends    = siteData{1}(Fnb,6)-200;
            sdat       = nans(num_tr, 1);
            pdat       = nans(num_tr, 1);
            if size(siteData{3},2)==2
                uu = 2;
            else
                uu = 1;
            end
            
            for tt = 1:num_tr
                
                % get pupil average(s)
                pi = max(500, find(isfinite(siteData{2}(Fnb(tt),:,4)),1));
                pe = find(isfinite(siteData{2}(Fnb(tt),:,4)),1,'last');
                
                if pe > pi + 1000
                    
                    % mean PD during trial
                    pdat(tt) = mean(siteData{2}(Fnb(tt),pi:pe,4));
                    
                    sp = siteData{3}{Fnb(tt), uu};
                    sdat(tt,1) = ...
                        sum(sp>=200&sp<=tr_ends(tt))./...
                        (tr_ends(tt)-200).*1000;
                end
            end
            
            % compute correlations
            Lu = isfinite(pdat) & isfinite(sdat);

            % correlation of residuals after linear
            % regressions
            tax = siteData{1}(Fnb(Lu),7);
            tax = tax - tax(1);
            X = [ones(sum(Lu),1) tax];
            % X = [ones(sum(Lu),1) (1:sum(Lu))'];
            % [Bp,BINTp,Rp] = regress(pdat(Lu), X);
            Bp = robustfit(X(:,2), pdat(Lu));
            Rp = pdat(Lu)-X*Bp;
            Bs = robustfit(X(:,2), sdat(Lu));
            Rs = sdat(Lu)-X*Bs;
            usdat{ss,2}(ff,8) = corr(Rp, Rs, 'type', 'Spearman');
        end
    end
end

for ss = 1:nus
    subplot(3,1,ss); cla reset; hold on;
%     Lf = isfinite(usdat{ss,2}(:,5)) & isfinite(usdat{ss,2}(:,6));
%     xs = usdat{ss,2}(Lf,7)-usdat{ss,2}(Lf,6);
    Lf = isfinite(usdat{ss,2}(:,8));
    xs = usdat{ss,2}(Lf,8);
    ys = usdat{ss,2}(Lf,5);
    plot(xs, ys, 'k.')
    [R,P] = corr(xs,ys,'type','Spearman');
    title(sprintf('r=%.2f, p=%.2f', R, P))
%    axis([0 1000 0 0.02]) 
end
    
