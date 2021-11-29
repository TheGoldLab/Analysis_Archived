function [th thd] = getML_RLSimulatedThForMTFine(Monk, N, recompute)
% simulate fine discrimination thresholds for model MT neurons

[h] = dirnames;
fn  = ['getML_RLSimulatedThForMTFine_' Monk '_' int2str(N) '.mat'];
savepath = [h, filesep, 'code', filesep, 'matlab', filesep, 'current', ...
    filesep, '01MTLIP', filesep, 'reports', filesep, '3-ReinforcementLearning', ...
    filesep, 'Preparation', filesep, 'mat', filesep, fn];



if recompute
    %%%
    %   Get statistics of MT responses for each neuron
    %   (0% response, pref gain, null gain, and fano factor)
    %%%
    bb   = 100;
    be   = 1000;
    bw   = 100;
    bs   = 25;
    bins = [zeros(size([bb+bw:bs:be]')) [bb+bw:bs:be]'];

    [pp, np, fp  ] = getML_SigNoise2([Monk 'PRe_MT.txt'], 1, [0 1], 1, 0);
    [p, n, f]      = getML_SigNoise2([Monk 'TRain_MT.txt'], 1, [0 1], 1, 0);

    % combine pre- and during training data
    PREF = [pp p];
    NULL = [np n];
    FANO = [fp f];

    COHS = [3.2 6.4 12.8 25.6 51.2 99.9]'/100;
    th   = nans(length(PREF),1);
    thd  = nans(length(PREF),1);
    for i = 1:length(PREF)
        if ~isnan(PREF(2,i)) & ~isnan(NULL(2,i)) & ~isnan(FANO(i)) & PREF(2,i)>0.5
            p    = nans(length(COHS),1);
            for j = 1:length(COHS)
                bl = nanmean([PREF(1,i) NULL(1,i)]);
                
                
                % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)        
                atune = normpdf(45, 10, 40);
                atune = atune/max(max(atune));          
                acoh  = (PREF(2,i)-NULL(2,i)).*atune+NULL(2,i);              
                m     = bl+COHS(j)*acoh;
                sd    = sqrt(FANO(i)*abs(m));
                rp    = m+sd*randn(N,1);

                % acoh will be a gaussian function, with peak at pref-slope, and baseline equal to null-slope (which can be a negative number)        
                atune = normpdf(45, -10, 40);
                atune = atune/max(max(atune));          
                acoh  = (PREF(2,i)-NULL(2,i)).*atune+NULL(2,i);              
                m     = bl+COHS(j)*acoh;
                sd    = sqrt(FANO(i)*abs(m));
                rn    = m+sd*randn(N,1);

                p(j) = rocN(rp,rn,100);
            end

            [fits, sems] = ctPsych_fit(@quick2,COHS, [p N*ones(size(p))]);
            th(i)  = fits(1);
            thd(i) = sems(1);
        end
    end

    save(savepath, 'th', 'thd')
    
    
    if 0
        [fitsp, semsp] = getML_neurometricROCT([Monk 'PRe_MT.txt'], bins, 1, [0 1], 0, 0);
        [fits, sems]   = getML_neurometricROCT([Monk 'TRain_MT.txt'], bins, 1, [0 1], 0, 0);
        th_ = [fitsp(1,:) fits(1,:)];
        
        cla
        hold on
        plot(th, th_, '.k')
        plot([0.1 1], [0.1 1], 'r')
        hold off
        set(gca, 'xscale', 'log', 'yscale', 'log')
    end
else
    load(savepath, 'th', 'thd')
end















