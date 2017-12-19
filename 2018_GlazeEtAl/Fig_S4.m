function Fig_S4
%% Supplementary Figure 4
%
% Additional choice variability in the Sampling Model
%
% Figure-generating code for:
%  Glaze CM, Filipowicz ALS, Kable JW, Balasubramanian V, and Gold JI
%  "A bias-variance trade-off governs individual differences in on-line 
%     learning in an unpredictable environment"
%

%% Set up the figure
wid     = 6; % total width
ht      = 6;
cols    = {1};
[axs,~] = getPLOT_axes(3, wid, ht, cols, 2, .5, [], 'Glaze et al', true);
set(axs,'Units','normalized');

% assume fixed H0, range of r0 and M
H0 = 0.5;
r0 = exp(linspace(0.7,7,5));
nM = 2.^(0:5);
N = 1000;

nr0 = length(r0);
nnM = length(nM);
sdat = nans(nr0, nnM);

for rr = 1:nr0
    for mm = 1:nnM
        
        % get H, J samples
        H_samples = betarnd(r0(rr)*H0, r0(rr)*(1-H0), N, nM(mm));        
        J_samples = log(H_samples./(1-H_samples));
        
        % find nearest to 0
        sdat(rr,mm) = sum((min(abs(J_samples),[],2).^2))./(N-1);
    end
end

axes(axs(1)); cla reset; hold on;
set(gca, 'FontSize', 12);
for rr = 1:nr0
    plot(nM, sdat(rr,:), 'k-', 'LineWidth', (nr0-rr+1)/2)
end
xlabel('Number of hazard-rate hypotheses');
ylabel('Variance of sampled log[H/(1?H)]')
legs = cell(nr0,1);
for rr = 1:nr0
    legs{rr} = num2str(log(r0(rr)));
end
legend(legs{:});
axis([0 32 0 3.2])
