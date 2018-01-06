function plotBIAS_bHist(bs, es, xs, ax, xlm)

axes(ax); cla reset; hold on;

Lsp = bs-1.96*es>0;
Lsn = bs+1.96*es<0;
disp([sum(Lsp) sum(Lsn) sum(isfinite(bs))])

gr  = 0.99*ones(1,3);

if nargin < 3 || isempty(xs)
    
    %% Histogram
    bax = (-4:.1:4)';
    bar(bax, hist(bs, bax),1,'k');

    h2=bar(bax, hist(bs(Lsp), bax),1,'r');
    set(h2,'EdgeColor','r');

    h3=bar(bax, hist(bs(Lsn), bax),1,'b');
    set(h3,'EdgeColor','b');

    disp(sprintf('n-:%d,n+:%d,n:%d', sum(bs+1.96*es<0), sum(bs-1.96*es>0), length(bs)))
    title(sprintf('%.3f [%.3f] (p=%.3f)', nanmedian(bs), iqr(bs(isfinite(bs))), signtest(bs(isfinite(bs)))))

    if nargin > 4 && ~isempty(xlm)
        xlim([-xlm xlm]);
    else
        xlim([bax(1) bax(end)]);
    end
    
else
    
    %% scatter
    xs(xs>1.5) = 1.5;
    plot([0 5], [0 0], 'k:');
    plot([xs xs]', [bs-es bs+es]', '-', 'Color', 0.8*ones(1,3));
    plot(xs(~Lsp&~Lsn), bs(~Lsp&~Lsn), 'kx', 'MarkerSize', 3);
    plot(xs(Lsp), bs(Lsp), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
    plot(xs(Lsn), bs(Lsn), 'ko', 'MarkerFaceColor', gr,  'MarkerSize', 3);
    
    A         = [ones(size(xs)) xs];
    Lg        = isfinite(xs) & isfinite(bs) & isfinite(es);
    [b, bint] = regressW(bs(Lg), es(Lg), A(Lg,:));
    xax       = (0:.1:5)';
    plot(xax, [ones(size(xax)) xax]*b, 'k--');
%    title(sprintf('b=%.2f [%.2f %.2f]', b(2), bint(2,1), bint(2,2)))
%    axis([xax(1) xax(end) -5 5]);
end
