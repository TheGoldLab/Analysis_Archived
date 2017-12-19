clear
clc
close all

logKspc = linspace(-10,log(.1),20);
logr0spc = linspace(log(.5),log(2000),50);
partMspc = [5 10 20 100 1000];
partMspc = [1000];

[logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc] = get_all_data;

alphamat = zeros(numel(partMspc),numel(logr0spc),numel(logKspc));
alphamatb = zeros(numel(logr0spc),numel(logKspc));

for logr0ind = 1:numel(logr0spc)
    for logKind = 1:numel(logKspc)
        params = [logr0spc(logr0ind) .5 logKspc(logKind)];
        alpha = get_mod_regression_fast(params,'bayesultimate',logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc);
        alphamatb(logr0ind,logKind) = mean(alpha,'omitnan');
       
        logr0ind
        logKind
        
%         for partMind = 1:numel(partMspc)
%             
%             params = [logr0spc(logr0ind) .5 logKspc(logKind) partMspc(partMind)];
%             
%             alpha = get_mod_regression_fast(params,'partfilt',logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc);
%             
%             alphamat(partMind,logr0ind,logKind) = mean(alpha,'omitnan');
%             
%         end
    end
end

