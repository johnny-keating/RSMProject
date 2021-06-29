function [z, x, pie, indices, exitflag] = rsm(A, b, c, m, n)
% Solves min cx s.t. Ax=b, x>=0
% exitflag is 0 if solved successfully, 1 if infeasible, -1 if unbounded
% Performs a Phase I procedure starting with an all artificial basis
% and then calls function simplex with modified leaving variable 
% criterion in Phase II. 

tolerance = 1.0E-6;
% Phase 1:
% Construct inputs for Phase I
phase = 1;
indices = (n + 1):(n + m);
cP1 = zeros((n + m), 1);
cP1(indices) = 1;
Bmatrix = eye(m);
AP1 = [A, Bmatrix];

% Input phase 1 structure into simplex method 
[z, x, pie, indices, exitflag] = PartCsimplexHelper(AP1, b, ...
    cP1, m, (n+m), Bmatrix, indices, phase);

% Check if problem is infeasible 
if z > tolerance
    exitflag = 1; 
else
    % Phase 2:
    phase = 2;
    AP2 = [A, AP1(:, indices(indices > n))];
    Bmatrix = AP1(:, indices);
    indices = [indices(indices <= n), (n + 1):n + ...
        length(indices(indices > n))];
    cP2 = [c; zeros(length(AP2) - length(A), 1)];
    [z, x, pie, indices, exitflag] = PartCsimplexHelper(AP2, b, cP2, m, ...
        n, Bmatrix, indices, phase);
end
temp = zeros(n, 1);
temp(indices <= n) = x(indices <= n);
x = temp;
end