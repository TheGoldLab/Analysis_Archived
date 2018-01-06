function [bs_,rocs_]   = getBIAS_LIPdata2(sid)
% function [bs_,rocs_] = getBIAS_LIPdata2(sid)

global FIRA

% get tin, wk
tin   = getBIAS_tinFromFIRA(sid);
wk    = getBIAS_wkFromFIRA(tin);

times = getFIRA_ecodesByName({'trg_on', 'dot_on', 'dot_off', 'fp_off', 'sac_lat'});
cohs  = getFIRA_ecodesByName('dot_coh');
cors  = getFIRA_ecodesByName('correct');
stin  = sign(cos(tin*pi/180));
chcs  = stin.*sign(cos(getFIRA_ecodesByName('choice')*pi/180));
dirs  = stin.*sign(cos(getFIRA_ecodesByName('dot_dir')*pi/180));
Lgood = getFIRA_ecodesByName('task')==6 & cors>=0;
Lcor  = cohs==0 | cors==1;
Ltin  = chcs==1;

% fit to PMF, with wk
ctdat = [cohs./100 (times(:,3)-times(:,2))./1000 dirs chcs cors];
fits  = ctPsych_fit(@ddExp5fz_F, ctdat(Lgood, 1:4), ctdat(Lgood, 5), ...
    [], [], [], [], wk(Lgood));
wk    = wk.*fits(4) + fits(5);

% save mean abs wk
bs_ = [nanmean(abs(wk(Lgood))) nanse(abs(wk(Lgood))) nanmean(wk(Lgood))];

% get data
bins = [(0:50:1400)', (100:50:1500)'];
sb   = size(bins,1);
rs   = nans(size(Lgood,1), 4+sb);
[r,zs,rs(Lgood,1)]        = getFIRA_rateByBin(Lgood, sid, times(Lgood,1)-300, times(Lgood,1));
[r,zs,rs(Lgood,2)]        = getFIRA_rateByBin(Lgood, sid, times(Lgood,1),     times(Lgood,2));
[r,zs,rs(Lgood,2+(1:sb))] = getFIRA_rateByBin(Lgood, sid, times(Lgood,2),     times(Lgood,3), [], bins);
[r,zs,rs(Lgood,end-1)]    = getFIRA_rateByBin(Lgood, sid, times(Lgood,4)-300, times(Lgood,4));
[r,zs,rs(Lgood,end)]      = getFIRA_rateByBin(Lgood, sid, times(Lgood,4)+times(Lgood,5)-200, ...
    times(Lgood,4)+times(Lgood,5));
    
% pres -- roc tin vs tout, sep by bias
rocs_ = nans(1,8,size(rs,2));
for rr = 1:size(rs,2)
    rocs_(1,:,rr) = [ ...
        rocN(rs(Lgood&Lcor&Ltin,rr),      rs(Lgood&Lcor&~Ltin,rr)), ...
        rocN(rs(Lgood&Lcor&Ltin&wk>0,rr), rs(Lgood&Lcor&~Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&Lcor&Ltin&wk<0,rr), rs(Lgood&Lcor&~Ltin&wk>0,rr)), ...
        rocN(rs(Lgood&Ltin,rr),           rs(Lgood&~Ltin,rr)), ...
        rocN(rs(Lgood&Ltin&wk>0,rr),      rs(Lgood&~Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&Ltin&wk<0,rr),      rs(Lgood&~Ltin&wk>0,rr)), ...
        rocN(rs(Lgood&Ltin&wk>0,rr),      rs(Lgood&Ltin&wk<0,rr)), ...
        rocN(rs(Lgood&~Ltin&wk>0,rr),     rs(Lgood&~Ltin&wk<0,rr))];
end
