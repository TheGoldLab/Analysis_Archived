function pHmat= calc_pH_partfilt(Hsampmat,Hspc)

[partM,N] = size(Hsampmat);
H_M = numel(Hspc);

pHmat = zeros(H_M,N);

for n = 1:N
    [N,xout] = hist(Hsampmat(:,n),Hspc);
    pHmat(:,n) = N;
end
