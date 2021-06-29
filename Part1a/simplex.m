function [z, x, pie, indices, exitflag] = simplex(A, b, c, m, n,...
    Bmatrix, indices)
% Solves min cx s.t Ax=b, x>=0
% starting with basic variables listed in vector indices
% and basis matrix Bmatrix
% exitflag is 0 if solved successfully, -1 if problem is unbounded
% returns optimal vector x and its value z, along with pi, and indices of
% basic variables. If the problem is unbounded simplex returns the last BFS
% found.
exitflag = 0; 
% Find values of the variables in the solution 
x = Bmatrix\b;
% Find the optimal value 
z = c(indices)'*x;
% Calculate pi
pie = (c(indices)'/Bmatrix)';
% Return the column as, cost coefficent cs and the column index s
[as, cs, s] = findenter(A, pie, c);

% Exits if no entering variable can be found. 
while s ~= 0
    % Return the index of leaving variable in the basis matrix
    [leave] = findleave(Bmatrix, as, x);
    
    % If the problem is unbounded. 
    if leave == 0
        temp = zeros(n, 1);
        temp(indices) = x;
        x = temp;
        exitflag = -1; 
        return 
    end
    
    % Update the indices, Bmatrix and create cb
    [Bmatrix, indices, cb] = update(Bmatrix, indices, c(indices), cs, s,...
        as, leave);
    % Same operation's as above
    x = Bmatrix\b;
    z = c(indices)'*x;
    pie = (c(indices)'/Bmatrix)';
    [as, cs, s] = findenter(A, pie, c);
end

% Set x to the full variable array
temp = zeros(n, 1);
temp(indices) = x;
x = temp;
end

