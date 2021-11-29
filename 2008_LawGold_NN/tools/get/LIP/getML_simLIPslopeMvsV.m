function [xsim, ysim] = getML_simLIPslopeMvsV(Monk, N, recompute)

savepath = ['/work/mirror_jeff/code/matlab/tmp-mat/getML_simLIPslopeMvsV_' Monk '.mat'];

if recompute
    [R m v] = getML_LIPSlopeDist(Monk, 0);

    % get simulated m/v for the mixture model
    coh  = [0 3.2 6.4 12.8 25.6 51.2 99.9];
    xsim = nans(N,2*length(coh));
    ysim = nans(N,2*length(coh));
    
    for i = 1:7
        fprintf(sprintf('%d-pref\n',i))
        %%% get sim for pref
        %
        % find session with the lowest and highest 10%
        ind = 1:length(m);
        x   = m(:,i);
        y   = v(:,i);

        if strcmp(Monk, 'Cy')
            ilo = ind(x<=myprctile(x,10));
            ihi = ind(x>=myprctile(x,90));
        else % remove some outliers for ZsaZsa
            ilo = ind(x>=myprctile(x,3) & x<=myprctile(x,13));
            ihi = ind(x>=myprctile(x,87) & x<=myprctile(x,97));
        end
        
        % do simulation
        lenlo = length(ilo);
        lenhi = length(ilo);
        M     = 20;
        for j = 1:N
            if ~mod(j,round(N/10))
                fprintf(sprintf('%d\n',j))
            end
            p         = round(M*rand);   % portion drawn from LO group
            sesLo     = ilo(unidrnd(lenlo,1,p));
            sesHi     = ihi(unidrnd(lenhi,1,M-p));
            L         = ismember(R(:,1),sesLo) & R(:,3)==coh(i) & R(:,4)==1 & R(:,5)==1;
            spLo      = R(L,6);
            L         = ismember(R(:,1),sesHi) & R(:,3)==coh(i) & R(:,4)==1 & R(:,5)==1;
            spHi      = R(L,6);
            xsim(j,i) = nanmean([spLo; spHi]);
            ysim(j,i) = nanvar([spLo; spHi]);
        end
        [xsim(:,i) I] = sort(xsim(:,i));
        ysim(:,i)     = ysim(I,i);
   
        %%% get sim for null
        %
        % find session with the lowest and highest 10%
        fprintf(sprintf('%d-null\n',i))
        
        ind = 1:length(m);
        x   = m(:,i+7);
        y   = v(:,i+7);
        
        if strcmp(Monk, 'Cy')
            ilo = ind(x<=myprctile(x,10));
            ihi = ind(x>=myprctile(x,90));
        else % remove some outliers for ZsaZsa
            ilo = ind(x>=myprctile(x,3) & x<=myprctile(x,13));
            ihi = ind(x>=myprctile(x,87) & x<=myprctile(x,97));
        end
        
        % do simulation
        lenlo = length(ilo);
        lenhi = length(ilo);
        M     = 20;
        for j = 1:N
            if ~mod(j,round(N/10))
                fprintf(sprintf('%d\n',j))
            end
            p           = round(M*rand);   % portion drawn from LO group
            sesLo       = ilo(unidrnd(lenlo,1,p));
            sesHi       = ihi(unidrnd(lenhi,1,M-p));
            L           = ismember(R(:,1),sesLo) & R(:,3)==coh(i) & R(:,4)==1 & R(:,5)==1;
            spLo        = R(L,6);
            L           = ismember(R(:,1),sesHi) & R(:,3)==coh(i) & R(:,4)==1 & R(:,5)==1;
            spHi        = R(L,6);
            xsim(j,i+7) = nanmean([spLo; spHi]);
            ysim(j,i+7) = nanvar([spLo; spHi]);
        end        
        [xsim(:,i+7) I] = sort(xsim(:,i+7));
        ysim(:,i+7)     = ysim(I,i+7);
        
   
    end  
    save(savepath, 'xsim', 'ysim')

else
    load(savepath)
end
