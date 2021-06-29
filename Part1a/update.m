function [Bmatrix, newIndices, newcb] = update(Bmatrix, indices, cb, cs,...
    s, as, leave)
% Bmatrix is current my m basis matrix
% indices is a column vector current identifiers for basic varaibles in
% order of B columns 
% cb is a column vector of basic costs in the order of B columns 
% as is the entering column
% s is the index of the entering variable
% leave is the column (p) of the basis matrix that must leave (not its
% variable index t)
% updates replaces column leave of Bmatrix with as to give newBmatrix
% replaces row leave of indices with enter to give newIndices
% replaces row leave of cb with cs to give new_cb

% Update and sort the indices, to help with Bland's Rule.
indices(leave) = s;
[newIndices, id] = sort(indices);
% Update and reorder the Bmatrix
Bmatrix(:, leave) = as;
Bmatrix = Bmatrix(:, id);
% Update and reorder the cb matrix
cb(leave) = cs;
newcb = cb(id);
end

