function [IBmatrix, indices, cb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave)
% updateGJ updates the inverse matrix by using the pivot operation
% Parameters:
% IBmatrix: The inverse Bmatrix to update 
% Indices: Identifies the basic variables in order of the Bmatrix
% cb: Column vector of basic costs in the order of B columns 
% as: The entering column
% s: The index of the entering variable
% Leave: (p) Column index of the Bmatrix that leaves. 

y = IBmatrix*as;
% Divide the leaving row by the the corresponding y value
IBmatrix(leave, :) = IBmatrix(leave, :)/y(leave);
% y value will become 1 in the leaving row
y(leave) = 1;

% Pivot
for i = 1:length(as)
    if i ~= leave
        IBmatrix(i, :) = IBmatrix(i, :) - IBmatrix(leave, :)*y(i);
    end 
end

% Update and reorder the indices
indices(leave) = s;
[indices, id] = sort(indices);
% Reorder the IBmatrix
IBmatrix = IBmatrix(id, :);
% Update and reorder the cost vector
cb(leave) = cs; 
cb = cb(id);
end

