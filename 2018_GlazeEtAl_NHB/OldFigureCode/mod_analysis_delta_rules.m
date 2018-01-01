clear
clc
close all

logammaspc = linspace(-10,log(.1),20);

alphavc_norm = zeros(length(logammaspc),1);
alphavc_leakint = alphavc_norm;
[logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc] = get_all_data;

for gammaind = 1:numel(logammaspc)
    params = [logammaspc(gammaind) .5];
    alpha = get_mod_regression_fast(params,'delta_rule_norm',logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc);
    alphavc_norm(gammaind) = mean(alpha,'omitnan');
    
    alpha = get_mod_regression_fast(params,'delta_rule_leakint',logJarr,xarr,l1arr,l2arr,musgnarr,sigmavc,subjvc);
    alphavc_leakint(gammaind) = mean(alpha,'omitnan');

    gammaind

end

figure; plot(logammaspc,alphavc_norm,'-','linewidth',2); hold on; plot(logammaspc,alphavc_leakint,'-','linewidth',2)
set(gca,'box','off','ticklength',[0,0],'fontsize',12)
legend({'delta rule on H, normative inference','delta rule on H, leaky integration'},'fontsize',12)

