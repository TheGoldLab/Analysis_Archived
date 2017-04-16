% Script for generating Figure 1, panels d-f, for
% Krishnamurthy, Nassar et al, 2017
%
% D. LOW RELIABILITY, HIGH RELEVANCE
% E. LOW RELEVANCE
% F. HIGH RELIABILITY, HIGH RELEVANCE
%
% Written by jigold

% x-axis
xax = -180:0.1:180;

% define likelihoods
L_SIGMA = 12;       % width 
L_MUS   = [9 67 39];  % means

% define prior means by hand
P_MUS       = [18 20 30];

% RELIABILITY: samples & generative sigma define "relative uncertainty"
P_SAMPLES   = [1 10 5];
G_SIGMA     = 8; % use this as generative width to make figs look nice

% RELEVANCE: define change-point probability by hand
CPP = [0.0001 0.9 0.0001];

% plot params
lw = 4; % line width
axlims = [-20 100];

for ii = 1:3

    % get panel
    subplot(3,1,ii); cla reset; hold on;
    
    % get/plot likelihood
    likelihood = normpdf(xax, L_MUS(ii), L_SIGMA);
    likelihood = likelihood./sum(likelihood);
    plot(xax, likelihood, 'b-', 'LineWidth', lw);

    % get/plot prior
    prior = normpdf(xax, P_MUS(ii), sqrt(G_SIGMA.^2 + (G_SIGMA/sqrt(P_SAMPLES(ii))).^2));
    prior = prior./sum(prior);
    plot(xax, prior, 'r-', 'LineWidth', lw);

    % posterior for realsies
    posterior = (1-CPP(ii))*likelihood.*prior + CPP(ii)*likelihood;
    posterior = posterior./sum(posterior);
    plot(xax, posterior, 'm-', 'LineWidth', lw);
    xlim(axlims);
    set(gca, 'FontSize', 14, 'yTick', [], 'XTick', []);
    ylabel('Probability');
end
xlabel('Simulated azimuthal location')

