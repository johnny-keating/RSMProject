function [z, x, pie, indices, exitflag] = PartCsimplexHelper(A, b, c, m,...
    n, Bmatrix, indices, phase)
% Solves min cx s.t Ax=b, x>=0
% starting with basic variables listed in vector indices
% and basis matrix Bmatrix
% exitflag is 0 if solved successfully, -1 if problem is unbounded
% returns optimal vector x and its value z, along with pi, and indices of
% basic variables. If the problem is unbounded simplex returns the last BFS
% found.

exitflag = 0;
% Create IBmatrix, the inverse of the B matrix
IBmatrix = inv(Bmatrix); 
% Find values of the variables in the solution 
x = IBmatrix*b;
% Find the optimal value 
z = c(indices).'*x;
% Calculate pi
pie = (c(indices)'*IBmatrix).';
% Return the column as, cost coefficent cs and the column index s
[as, cs, s] = findenter(A, pie, c, indices);

% Exits if no entering variable can be found. 
while s ~= 0
    % Return the index of leaving variable in the basis matrix
    artificial = indices > n;
    [leave] = ExtendedFindLeave(IBmatrix, as, x, phase, artificial);
    % If the problem is unbounded. 
    if leave == 0
        exitflag = -1; 
        return 
    end
    
    % Update the indices, Bmatrix and create cb
    [IBmatrix, indices, cb] = updateGJ(IBmatrix, indices, c(indices),...
        cs, as, s, leave);
    % Same operation as above
    x = IBmatrix*b;
    z = c(indices)'*x;
    pie = (c(indices)'*IBmatrix)';
    [as, cs, s] = findenter(A, pie, c, indices);
end
end

