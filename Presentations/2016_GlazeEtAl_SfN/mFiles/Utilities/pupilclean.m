function [cleaned_data,dropout] = pupilclean(ldil, threshold)
% PUPIL CLEAN    Clean pupil data with  blinks.
% Argument 1 = your pupil data vector
% Argument 2 = a z-score used to eliminate unusual points in data. Low
% z-score to exclude much data, high z-score to exclude only the most
% unusual data.
% By M. Kabir, reach out for help at kabir.naim@gmail.com

 
    orientation = size(ldil);
    
    if orientation(2) > orientation(1)
        ldil = ldil';
    end
%     
%     negs = ldil == -1; %find where validity code is bad
% 
%     %Find places where the derivative of the dilation is suspiciously large
%     deriv = [0; abs(diff(ldil))];
%     threshold = mean(deriv) + z*std(deriv); %set threshold above which things are disincluded
    bad_idx = ldil < threshold;
%     bad_idx = bad_idx + negs; bad_idx = bad_idx > 0;

    idx = find(~bad_idx); %get good data indices

    %Interpolate areas where there is bad data, using non-bad data as the
    %interpolation function. idx designates the indices of values looked at in
    %the 'function' ldil(idx), in order to fill in query points provided by 
    %find(negs) 
    ldil(bad_idx) = interp1(idx, ldil(idx), find(bad_idx));
    
    dropout = mean(bad_idx);
    
    %Deal with any NaNs from interpolation, substitute with signal mean
    nans = isnan(ldil);
    averaged = mean(ldil(~nans));
    ldil(nans) = ones(1, length(ldil(nans)))*averaged;
    
    cleaned_data = ldil;
end
    

