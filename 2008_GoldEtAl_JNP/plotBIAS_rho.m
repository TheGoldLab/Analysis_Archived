function plotBIAS_rho(xs, ys, es, ps, ax, axlm)

if nargin < 6 || isempty(axlm)
    axlm = [-.7 .7];
elseif isscalar(axlm)
    axlm = [-axlm axlm];
end

axes(ax); cla reset; hold on;
Lp = ps<0.05 & ys>0;
Ln = ps<0.05 & ys<0;
Lf = isfinite(ys);
g1 = 0.8*ones(1,3);
g2 = 0.99*ones(1,3);

if ~isempty(xs)

    %% scatter
    plot([xs xs]', [ys-es ys+es]', '-', 'Color', g1);
    plot(xs(~Ln&~Lp), ys(~Ln&~Lp), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 2);
    plot([0 5], [0 0], 'k:');
    plot(xs(Lp), ys(Lp), 'ko', 'MarkerFaceColor', g1, 'MarkerSize', 2);
    plot(xs(Ln), ys(Ln), 'ko', 'MarkerFaceColor', g1, 'MarkerSize', 2);
    A        = [ones(size(xs)) xs];
    [b,bint] = regressW(ys(Lf), es(Lf), A(Lf,:));
    xx = (0:0.1:5)';
    plot(xx, [ones(size(xx)) xx]*b, 'k--');
    disp([nanmedian(ys) signtest(ys(Lf)) sum(Lp) sum(Ln) sum(Lf) b(2) bint(2,1) bint(2,2)])
    ylim(axlm);

    plot([0 2], nanmedian(ys).*[1 1], 'r--');
else
    
    %% Histogram
    xax = (-1:.1:1)';
    bar(xax, hist(ys, xax),1,'k');

    h2=bar(xax, hist(ys(Lp), xax),1);
    set(h2,'EdgeColor', g1, 'FaceColor', g1);

    h3=bar(xax, hist(ys(Ln), xax),1);
    set(h3,'EdgeColor', g1, 'FaceColor', g1);

    disp([nanmedian(ys) signtest(ys(Lf)) sum(Lp) sum(Ln) sum(Lf)])
    xlim(axlm);
end

