function [rocs, ns, sems, p] = getML_MTChoiceProbability(monk, recompute)
% get choice probability of MT neurons

% TRain (choice probability)
fn     = [monk 'TRain_MT.txt'];
a      = getML_txt(fn);
fname  = a.data{strcmp(a.name,'dat_fn')};
ddir   = a.data{strcmp(a.name,'ddir')};
uid    = a.data{strcmp(a.name,'uid')};
usable = a.data{strcmp(a.name,'usable')};
train  = a.data{strcmp(a.name,'train')};
ses    = a.data{strcmp(a.name,'session')};
tkflg  = a.data{strcmp(a.name,'taskflag')};


[hdir, ldir, cdir, tdir] = dirnames;
savepath = [tdir '/getML_MTChoiceProbability_' monk '.mat'];


if recompute
    rocs = nans(length(fname), 7);
    ns   = nans(length(fname), 7);
    sems = nans(length(fname), 7);
    p    = nans(length(fname), 7);
    coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9]';

    global FIRA
    for i = 1:length(fname)
        if usable(i)==1 & train(i)==1
            fprintf('%d: %s\n', i, fname{i})
            openFIRA(fname{i})

            trials      = find(~isnan(FIRA.ecodes.data(:,1)));

            %% get selection arrays
            % for dot_dir
            [Ld, Ud]    = selectFIRA_trialsByUniqueID('dot_dir', trials);
            Ld          = Ld(:,[find(Ud==round(ddir(i))) find(Ud==round(mod(ddir(i)+180,360)))]);

            % for choice
            [Lch, Uch]  = selectFIRA_trialsByUniqueID('choice', trials);
            Lch         = Lch(:,[find(Uch==round(ddir(i))) find(Uch==round(mod(ddir(i)+180,360)))]);

            % for dot coherence
            [Lc, Uc]    = selectFIRA_trialsByUniqueID('dot_coh', trials);

            % for correct/incorrect (only include crt and incrt trials, ignore
            % broken fix & no choice trials)
            Lct         = FIRA.ecodes.data(trials, getFIRA_ecodeColumnByName('correct'));
            Lct         = Lct>=0;

            % for task
            [Ltk, Utk]  = selectFIRA_trialsByUniqueID('task', trials);
            if tkflg(i)<=2
                Ltk     = Ltk(:,Utk==3);
            else
                Ltk     = Ltk(:,Utk==6);
            end

            for j = 1:length(Uc)
                % choose only trials when dots moved to the pref direction, except for 0% coherence
                % i.e marginal choice probability
                if Uc(j)~=0
                    Lc(:,j) = Lc(:,j) & Lct & Ltk & Ld(:,1);
                else
                    Lc(:,j) = Lc(:,j) & Lct & Ltk;
                end
            end

            spikei      = getFIRA_spikeByID(uid(i));

            if strcmp(monk, 'Cy') 
                minn = 2;
                bt          = getFIRA_ecodeTimesByName('dot_on',  300, trials);
                et          = getFIRA_ecodeTimesByName('dot_off', 0, trials);           
            else
                minn = 5;
                bt          = getFIRA_ecodeTimesByName('dot_on',  0, trials);
                et          = getFIRA_ecodeTimesByName('dot_off', 0, trials);            
            end
                
            % get grid of ROCs
            [rocs_, ns_, sems_] = getFIRA_neurometricROCT(trials, Lch, Lc(:,ismember(Uc,coh)),...
                spikei, bt, et, [], minn);

            % put roc in the right place, because sometimes some coh are not
            % used in training
            c = [0 3.2 6.4 12.8 25.6 51.2 99.9];
            Uc(~ismember(Uc,c))=[]; %remove strange coh
            for j = 1:length(Uc)
                Ic = find(c==Uc(j));
                rocs(i,Ic) = rocs_(j);
                ns(i,Ic)   = ns_(j);
                sems(i,Ic) = sems_(j);
                p(i,Ic)    = nan;
                
                %p(i,Ic)    = p_(j);
            end
        end
    end
    
    save(savepath, 'rocs', 'ns', 'sems', 'p')
else

    load(savepath)
end

% % plot result and see
% figure
% hold on
% for i = 1:length(roc)
%     plot(ses(i), roc(i,1), 'o')
% end
% hold off
% ylim([0.3 0.8])
%
% % plot result and see
% figure
% hold on
% for i = 1:length(roc)
%     plot(ses(i), nanmean(roc(i,1)), 'o')
% end
% hold off
% ylim([0.3 0.8])
%
%
%
%
% boxplot(roc(:,1:7)')
% %xlim([min(ses), max(ses)])
% ylim([0.3 0.8])
%
%


