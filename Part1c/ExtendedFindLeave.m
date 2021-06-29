function [leave] = ExtendedFindLeave(IBmatrix, as, xb, phase, artificial)
% Given the entering column as and the vector xb of basic variables
% findleave finds a leaving column of basis matrix Bmatrix
% It returns 0 if no column can be found (i.e unbounded)
% leave = p indicates the pth column of the Bmatrix leaves.

most_important_vector = IBmatrix*as;
% If theres some artificials and Phase II remove them even if ratio is 0.
if any(artificial == 1) && (phase == 2)
    AVindex = find(artificial == 1);
    for i = AVindex
        if most_important_vector(i) ~= 0
            leave = AVindex(i);
            return 
        end
    end
end
          
% Check for unboundness
if all(most_important_vector <= 0)
    leave = 0;
    return
end

% Check if divinding by zero or negative number if so set to NaN
most_important_vector(most_important_vector <= 0) = NaN;
% Calculate Ratio
ratio = xb./most_important_vector;
% Return the minimum of the ratio
[~, leave] = min(ratio);
end
