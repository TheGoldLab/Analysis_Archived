clear
clc

cd /Users/lab/Documents/MATLAB/

variance = repmat({'high';'low'},5,1);
block_len = [repmat({'short'},4,1);repmat({'long'},6,1)];
type = repmat({'signaled';'signaled';'unsignaled';'unsignaled'},2,1);
type = [type;{'mixed';'mixed'}];

varvc = repmat([10;sqrt(150)],5,1);
taulen = [5*ones(4,1);2*ones(6,1)];

dirs = {'SIGNALED','UNSIGNALED','SEMISIGNALED'};

condN = zeros(length(type),1);

subjall = cell(100,1);
flindall = 0;

for dirind = 1:length(dirs)
    dirnm = ['./data_2AFC/' dirs{dirind} '/'];
    fls = dir(dirnm);
    
    for flind = 3:length(fls)
        flnm = fls(flind).name;
        
        if strcmp(flnm(end-2:end),'mat')
            load([dirnm flnm]);
            
            if data.N>200
                flindall = flindall + 1;
                subj = flnm(5:6);
                subjall{flindall} = subj;
                
                switch dirind
                    case 1
                        if data.sigma==sqrt(150)
                            if length(data.tau)>2
                                condN(1) = condN(1)+1;
                            else
                                condN(5) = condN(5)+1;
                            end
                        elseif data.sigma==10
                            if length(data.tau)>2
                                condN(2) = condN(2)+1;
                            else
                                condN(6) = condN(6)+1;
                            end
                        end
                        
                        
                    case 2
                        if data.sigma==sqrt(150)
                            if length(data.tau)>2
                                condN(3) = condN(3)+1;
                            else
                                condN(7) = condN(7)+1;
                            end
                        elseif data.sigma==10
                            if length(data.tau)>2
                                condN(4) = condN(4)+1;
                            else
                                condN(8) = condN(8)+1;
                            end
                        end
                        
                        
                        
                    case 3
                        if data.sigma==sqrt(150)
                            condN(9) = condN(9)+1;
                        elseif data.sigma==10
                            condN(10) = condN(10)+1;
                        end
                        
                end
                
            end
            
            clear data
            
        end
        
    end
    
end

subjall = subjall(1:flindall);

subjids = unique(subjall);

save subj_info.mat subjids