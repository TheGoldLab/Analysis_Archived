function fig_ = figBIAS_LIPROC(num)
%function fig_ = figBIAS_LIPROC(num)
%

if nargin < 1 || isempty(num)
    num = 14;
end

monkc = {'Cyrus' 'Cy' 3; 'ZsaZsa' 'ZZ' 4};
monkn = 2;

%% Set up Fig
% units should be in inches, from wysifig
wid  = 4; % total width
hts  = 0.8;
cols = {2,2,2};
[axs,fig_] = getBIAS_axes(num, wid, hts, cols);
set(axs,'Units','Normalized');

% get the data
rocdat = FS_loadProjectFile('2008_Bias', 'figBIAS_LIPROC');

if isempty(rocdat)

    rocdat = cell(size(monkc,1), 1);
    ddir   = '/Users/jigold/Desktop/data/MTLIP/';

    for mm = 1:monkn

        a        = getML_txt([ddir monkc{mm,2} 'TRain_LIP.txt']);
        fn       = a.data{strcmp(a.name,'dat_fn')};
        usable   = find(a.data{strcmp(a.name,'usable')}==1);
        uid      = a.data{strcmp(a.name,'uid')};

        % rocdat{mm,1} = nans(length(usable), 27, 7, 3, 2);
        rocdat{mm,1} = nans(length(usable), 27, 6);

        for ss = 1:length(usable) %ss=100 %length(ses):-1:1

            disp(fn{usable(ss)})
            openFIRA([ddir monkc{mm,1} '/' fn{usable(ss)}]);

            % get spike id, slopes
            sid   = getFIRA_spikeByID(uid(usable(ss)));
            rocdat{mm,1}(ss,:,:,:,:) = getBIAS_LIPdata6(sid);
        end
    end

    FS_saveProjectFile('2008_Bias', 'figBIAS_LIPROC', rocdat);
end

tax = mean( [(0:50:1300)', (200:50:1500)'], 2);
Lt  = tax < 900;
xs  = tax(Lt);
for mm = 1:2
    ses   = (1:size(rocdat{mm,1},1))';
    sbins = prctile(ses, linspace(0,100,4));
    for bb = 1:3
        Lses = ses>=sbins(bb) & ses<=sbins(bb+1);
        tdat = nans(sum(Lt), 3, 3);
        for tt = 1:sum(Lt)
            for pp = 1:2
                b1 = rocdat{mm}(Lses,tt,(pp-1)*2+1);
                b2 = rocdat{mm}(Lses,tt,(pp-1)*2+2);
                Lg = isfinite(b1) & isfinite(b2);
                tdat(tt,:,pp) = [ ...
                    nanmean(b2-b1), nanse(b2-b1), signtest(b1(Lg),b2(Lg))];
            end
            tdat(tt,1,3) = [nanmean(rocdat{mm}(Lses,tt,6))];
        end


        axes(axs((bb-1)*2+mm)); cla reset; hold on;
        %        plot(xs, tdat(:,1,3), '--', 'Color', 0.5*ones(3,1));
        %        axis([xs(1) xs(end) 0 1])
        plot(xs([1 end]), [0 0], 'k:');
        for pp = 1:2
            ys = tdat(:,1,pp);
            es = tdat(:,2,pp);
            plot([xs xs]', [ys-es ys+es]', '-', 'Color', 0.7*ones(3,1));
            plot(xs, ys, 'k-', 'LineWidth', 3-pp);
            Lp = tdat(:,3,pp) < 0.05;
            plot(xs(Lp), ys(Lp), 'k*');
        end
        axis([xs(1) xs(end) -0.1 0.17]);
    end
end
