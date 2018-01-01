function Fig_S1
%% Supplementary Figure 1
%
% Performance trade-offs
%
% Figure-generating code for
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% Set up figure
gr      = 0.6.*ones(3,1);
wid     = 11.6; % total width
ht      = 4;
cols    = {2};
[axs,~] = getPLOT_axes(1, wid, ht, cols, 1.5, 1.5, [], 'Glaze et al', true);
set(axs,'Units','normalized','FontSize', 12);

%% Where to find the data
[file_list, data_dir,raw_data_dir] = getDataInfo;

%% Loop through the subjects to collect data for analysis
num_subjects = length(file_list);
all_data     = [];
MU_DIST      = 30;

for ss = 1:num_subjects
    
    % get the data file
    load(fullfile(raw_data_dir, file_list{ss}))
    
    % loop through the sessions for this subject
    num_sessions = session_ind - 1;
    for ii = 1:num_sessions
       
       % get the appropriate data
       eval(['data=data' num2str(ii) ';'])
       
       % useful stuff
       midpt = mean(data.muall(:,1));
       musgn = sign(data.muall(data.muinds,1)-midpt);
       musgn(data.signaled==0) = 0;
       musgn = [0;musgn(1:end-1)];
       if data.H(1)>1
          data.H = 1./data.H;
       end
       data.H = round(data.H,2);
              
       % Collect:
       %   1. Likelihood of x given left
       %   2. Likelihood of x given right
       %   3. Choice (0/1)
       %   4. logJ0 = log(H/(1-H))
       %   5. Signaled mean from previous trial (if given)
       %   6. Generative mean
       all_data = cat(1, all_data, [ ...
          normpdf(data.X(:,1)-midpt, MU_DIST/2, data.sigma), ...
          normpdf(data.X(:,1)-midpt,-MU_DIST/2, data.sigma), ...
          double(data.pred==2), ...
          log(data.H(:,1)./(1-data.H(:,1))), ...
          musgn, ...
          (data.muinds-1.5).*2]);
    end
end

%% A: predicted performance, using real and swapped choice variability
%
% Get the regress fits, columns are:
%   1. J_default
%   2. H_slope
%   3. noise in the DV
%   4. lapse rate
load(fullfile(data_dir, 'adaptivityModelFits.mat'), 'fits');
rfits = fits;
nfits = size(fits,1);
fdat  = nans(nfits,2,2);
[Y,I] = sort(fits(:,2));
Lc1   = all_data(:,6)==1;
Lj    = abs(all_data(:,4))>2; % sort by H near/far from 0.5

% simulate, assuming no lapses
for ff = 1:nfits
    
    J_default = fits(I(ff),1);
    H_slope   = fits(I(ff),2);
    eta       = fits(I([ff nfits-ff+1]),3);
    
    % simulate performance with real, eta-swapped fits
    for xx = 1:2
        
        % L = Log-posterior-odds
        J       = J_default + H_slope*all_data(:,4);
        Hvc     = 1./(1+exp(-J));
        [q1,q2] = dsprt_mixed_Hvc2(all_data(:,1), ...
            all_data(:,2), Hvc, all_data(:,5));
        L       = log(q1)-log(q2);
        
        % LogLR from predicted and actual choices
        choicehat = 1./(1+exp(-L/eta(xx)));
        fdat(ff,xx,1) = mean([choicehat(Lc1); 1-choicehat(~Lc1)]).*100;
        
        % separated by H
        if xx==1
            fdat(ff,1,2) = mean([choicehat(Lc1& Lj); 1-choicehat(~Lc1& Lj)]).*100;
            fdat(ff,2,2) = mean([choicehat(Lc1&~Lj); 1-choicehat(~Lc1&~Lj)]).*100;
        end
    end
end

% stats for m_H, H_default
hd = fits(I,1);
mh = fits(I,2);
ys = fdat(:,1,1);
[R1,P1] = corr(hd,ys,'type','Spearman');
[R2,P2] = corr(mh,ys,'type','Spearman');
disp(sprintf('Hdefault vs performance: R=%.4f, P=%.4f',R1,P1))
disp(sprintf('mH       vs performance: R=%.4f, P=%.4f',R2,P2))

%% Panel A
axes(axs(1)); cla reset; hold on;
xs  = Y;
y1s = fdat(:,1,1);
y2s = fdat(:,2,1);
plot(xs, y1s, 'ko');
plot(xs, y2s, 'o', 'Color', gr, 'MarkerFaceColor', gr);
hs=lsline;
set(hs, 'LineWidth', 1.5)
xlabel('Adaptivity');
ylabel('Simulated peformance (% correct)');
axis([-0.1 1.1 90 97]);

% make linear model, for stats
Xv = cat(2, [xs; xs], [xs; zeros(size(xs))], [ones(size(xs)); zeros(size(xs))]);
Yv = [y1s; y2s];
LM = fitlm(Xv, Yv);

% just reg'lar
LM = fitlm(xs,y1s)

% just swapped
LM = fitlm(xs,y2s)

%% Panel B
axes(axs(2)); cla reset; hold on;
plot(Y, fdat(:,1,2), 'ks', 'MarkerSize', 7);
plot(Y, fdat(:,2,2), 'k+', 'MarkerSize', 7);
Lg = isfinite(Y) & isfinite(fdat(:,1,2)) & isfinite(fdat(:,2,2));
[R1,P1] = corr(Y(Lg), fdat(Lg,1,2), 'type', 'Spearman');
[R2,P2] = corr(Y(Lg), fdat(Lg,2,2), 'type', 'Spearman');
title(sprintf('R1=%.4f,P1=%.2g, R2=%.4f,P2=%.2g', R1,P1,R2,P2))
axis([0 1.1 90 97])


