function [leave] = findleave(Bmatrix, as, xb)
% Given: Entering column as; Vector xb of basic variables; Basis Matrix
% findleave finds a leaving column of basis matrix Bmatrix
% It returns 0 if no column can be found (i.e unbounded)
% leave = p indicates the pth column of the Bmatrix leaves.

most_important_vector = Bmatrix\as;

% Check for unboundness
if all(most_important_vector <= 0)
    leave = 0;
    return
end

% Check if dividing by zero or negative number if so set to NaN
most_important_vector(most_important_vector <= 0) = NaN;
% Calculate Ratio
ratio = xb./most_important_vector;
% Return the minimum of the ratio
[~, leave] = min(ratio);
end
