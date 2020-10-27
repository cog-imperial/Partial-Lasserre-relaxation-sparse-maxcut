function clique = gen_cliques(W)
%%
% Taken from SparsePOP version 3.01 by Hayato Waki, Sunyoung Kim, Masakazu Kojima, Masakazu Muramatsu, Hiroshi Sugimoto and Makoto Yamashita.
% This function is taken from lines 97 to 135 of the file genClique.m contained in SparsePOP301/subPrograms/Mfiles
% Availability: https://sparsepop.sourceforge.io/
% Refer to the license in the file SparsePOP301/readme.txt for usage of this function.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generates maximal cliques of a chordal extension of W.
nDim = size(W,1);
W_aux = spones(W);
rmat = W_aux + 5*nDim*speye(nDim);
% Step 2
% Chordal Extension by Cholesky decomposition
% minimum degree ordering
I = symamd(rmat);
%% cholesky decomposition
[R,p] = chol(rmat(I,I));
if (p > 0)
    error('Correlative sparsity matrix is not positive definite.');
end

% Step3
% Finding the maximal cliques
%

% put 1 for nonzero element ofR
Cliques = spones(R);
[value,orig_idx] = sort(I);
remainIdx = 1;
for i=2:nDim
    checkSet = Cliques(i,i:nDim);
    one = find(checkSet);
    noOfone = length(one);
    cliqueResult = Cliques(1:i-1,i:nDim)*checkSet';
    yesno = find(cliqueResult == noOfone);
    %
    % Remove the set included into other set.
    %
    if ~any(yesno)
        remainIdx = [remainIdx;i];
    end
end
clique = Cliques(remainIdx,orig_idx);

end
