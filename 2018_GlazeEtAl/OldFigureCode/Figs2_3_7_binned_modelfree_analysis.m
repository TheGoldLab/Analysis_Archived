clear
clc
close all

% MAIN ROUTINE FOR GENERATING MODEL-FREE ANALYSIS, BINNED BY EITHER
% INDIVIDUAL ADAPTABILITY MEASURE OR PRIOR PRECISION

dirnm = pwd;
data_dir = [dirnm(1:strfind(dirnm,'Dropbox')+19) filesep 'fit_data' filesep];

% THRESHOLD FOR STRENGTH OF INITIAL EVIDENCE FOR ONE SIDE (ABS(LLR) >=
% THRESHOLD FOLLOWED BY PARTICULAR DATA SEQUENCES)
initialcutoff = 2.5;


% THRESHOLD FOR STRENGTH OF SUBSEQUENT EVIDENCE FOR SWITCHED SIDE (ABS(LLR) <=
% THRESHOLD)
evcutoff = 2.5;

% GET BASIC STATS AND FOCUS ON SUBJECTS WHO COMPLETED >1 SESSION AND OMIT
% THE ONE WE EXPERIMENTED WITH WHO STARTED OFF WITH UNSIGNALED DATA AT
% BEGINNING
[substrct,sessionNvc,sigmamedvc,Hlowvc,Hhighvc,signaledmedvc,signaledfirstvc] = get_subjstats;
indsval = (sessionNvc>1 & signaledfirstvc>0);

% GET ALL DATA AND PARALLEL VECTOR OF LEARNING PARAMETERS CONCATENATED
% ACROSS SUBJECTS
[correctvc,RTvc,subjvc,Hvc,rvc,musgnvc,evc,initialevc,blockindvc,isprevc,sessionindsvc,paramat,alphavc,datacorrectvc,xvc,choicevc,xprevc] = get_modelfree_data('data');

RTmed = grpmedian(RTvc,subjvc);
RTdevc = RTvc;

subjU = unique(subjvc);

for subjind = 1:numel(subjU)
   indsval = subjvc==subjU(subjind) & RTvc<3000;
   RTdevc(indsval) = RTvc(indsval) - median(RTvc(indsval));
end

RTdevc = RTvc;


% COMPUTE PR(SWITCH) AS FUNCTION OF TRIALS POST CHANGE OVER DIFFERENT
% HAZARD RATE REGIMES, WHEN INITIAL EVIDENCE FOR PRE-CHANGE SIDE WAS STRONG
Hbins = [0 .2 .8 1];
rbins = [1 2 3 4 inf];

adaptvc = alphavc;
adaptbins = [-inf .1 .3 inf];
% 
% adaptvc = paramat(:,end-2);
% adaptbins = [-inf 2 4 inf];

pcorrectmat = zeros(numel(adaptbins)-1,numel(Hbins)-1,numel(rbins))+nan;
pcorrectmat(:,:,1)=0;
RTmat = pcorrectmat;
pcorrectsemat = pcorrectmat;
runlenmat = RTmat;
Nmat = RTmat;

trialindsval = musgnvc==0 & initialevc>initialcutoff & blockindvc>0 & evc<evcutoff;

for adaptind = 1:numel(adaptbins)-1
    adaptindsval = adaptvc>=adaptbins(adaptind) & adaptvc<adaptbins(adaptind+1);
    for Hind = 1:numel(Hbins)-1
        Hindsval = Hvc>=Hbins(Hind) & Hvc<Hbins(Hind+1);
        RTmed0 = median(RTdevc(trialindsval & adaptindsval & Hindsval));
        
        ispreinds = blockindvc>0 & adaptindsval & Hindsval & isprevc & evc>initialcutoff & sessionindsvc>0;
        pcorrectmat(adaptind,Hind,1) = 1 - mean(correctvc(ispreinds));
        RTmat(adaptind,Hind,1) = median(RTdevc(ispreinds & RTvc<3000));
        Nmat(adaptind,Hind,1) = sum(ispreinds);
        
        for rind = 1:numel(rbins)-1
            rindsval = rvc>=rbins(rind) & rvc<rbins(rind+1);
            indsval = trialindsval & adaptindsval & Hindsval & rindsval;
            
            if sum(indsval)>5
                pcorrectmat(adaptind,Hind,rind+1) = mean(correctvc(indsval));
                pcorrectsemat(adaptind,Hind,rind+1) = std(bootstrp(200,@mean,correctvc(indsval)));
                RTmat(adaptind,Hind,rind+1) = median(RTdevc(indsval & RTvc<3000  & sessionindsvc>0));
            end
            
            Nmat(adaptind,Hind,rind+1) = sum(indsval);
        end
    end
end


figure
for i = 1:numel(adaptbins)-1
    subplot(1,numel(adaptbins)-1,i)
    errorbar(repmat([1:numel(rbins)],3,1)',squeeze(pcorrectmat(i,:,:))',squeeze(pcorrectsemat(i,:,:))','-','linewidth',1.5)
    
    ylim([-.05 1.05])
    xlim([.5 numel(rbins)+.5])
    set(gca,'box','off','ticklength',[0 0],'fontsize',12,'xtick',1:numel(rbins),'xticklabel',{'pre','1','2','3','4+'})

end

set(gcf,'Position',[440   682   560   116])

xbins = [-50 -15:3:15 50];


% COMPUTE PR(SWITCH) AS FUNCTION OF EVIDENCE FOR SWITCHED SIDE OVER
% DIFFERENT HAZARD RATE REGIMES. HERE WE ARE GOING TO BIN HAZARD RATE MORE
% NARROWLY BECAUSE THE PURPOSE IS TO SHOW CHOICE FUNCTION SLOPE (NOISE)
% DIFFERENCES AND WE DO NOT WANT TO BIAS THE QUALITATIVE ANALYSIS BY
% LUMPING DIFFERENT HAZARD RATES TOGETHER INTO THE SAME REGIME (THERE MAY
% STILL BE ISSUES AND IF A REVIEWER ASKS WE CAN JUST SHOW SAME).

Hbins = [0 .07 .15 .35 .55 .8 .92 1];

pcorrectmat2 = zeros(numel(adaptbins)-1,numel(Hbins)-1,numel(xbins)-1)+nan;
pcorrectsemat2 = pcorrectmat2;
Nmat2 = pcorrectmat2;
xvals = zeros(numel(xbins)-1,1);
Nvals = xvals;

trialindsval = blockindvc>100 & (musgnvc~=0 | abs(xprevc)>15);
choicevc(xprevc>0) = 1-choicevc(xprevc>0);
xvc = -1 * xvc .* sign(xprevc);

for adaptind = 1:numel(adaptbins)-1
    adaptindsval = adaptvc>=adaptbins(adaptind) & adaptvc<adaptbins(adaptind+1);
    for Hind = 1:numel(Hbins)-1
        Hindsval = Hvc>=Hbins(Hind) & Hvc<Hbins(Hind+1);

        for xind = 1:numel(xbins)-1
            xindsval = xvc>=xbins(xind) & xvc<xbins(xind+1);
            indsval = trialindsval & adaptindsval & Hindsval & xindsval;
            
            if sum(indsval)>0
                pcorrectmat2(adaptind,Hind,xind) = mean(choicevc(indsval));
                pcorrectsemat2(adaptind,Hind,xind) = std(bootstrp(200,@mean,choicevc(indsval)));
                xvals(xind) = xvals(xind) + sum(xvc(indsval));
                Nvals(xind) = Nvals(xind) + sum(indsval);
            end
            
            Nmat2(adaptind,Hind,xind) = sum(indsval);
        end
    end
end

xvals = xvals./Nvals;
xvals = ceil(xvals);

pcorrectmat2 = pcorrectmat2(:,[1 4 7],:);
pcorrectsemat2 = pcorrectsemat2(:,[1 4 7],:);

figure
for i = 1:numel(adaptbins)-1
    subplot(1,numel(adaptbins)-1,i)
    errorbar(repmat(xvals',3,1)',squeeze(pcorrectmat2(i,:,:))',squeeze(pcorrectsemat2(i,:,:))','-','linewidth',1.5)
    
    ylim([-.05 1.05])
    xlim([-30 30])
    set(gca,'box','off','ticklength',[0 0],'fontsize',12)

end

set(gcf,'Position',[440   682   560   116])
