function [as, cs, s] = findenter(A, pie, c)
% Given the complete m by n matrix Amatrix,
% the complete cost vector c with n components,
% the vector pie with m components
% findenter finds the index of the entering variable and its column
% It returns the column as, its cost coefficient cs, and its column index s
% Returns s=0 if no entering variable can be found (i.e. optimial)
% This will happen when minimum reduced cost > tolerance
% where tolerance = -1.0E-6

tolerance = -1.0E-6;
% Calculate the reduced costs vector
rc = c' - pie'*A;

% Check optimiality criteria
if all(rc > tolerance)
    % If basic matrix is optimal set the following
    s = 0;
    as = 0;
    cs = 0;
else
    % If basic matrix is not optimal set entering variable to lowest
    % index which has a negative reduced cost. Bland's Rule.
    s = find(rc <= tolerance, 1);
    as = A(:, s);
    cs = c(s);
end 
end

