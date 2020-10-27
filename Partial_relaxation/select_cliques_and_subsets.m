function [] = select_cliques_and_subsets(problem, r, p, H)
%%
% This function calculates maximal cliques (and/or subsets of the maximal cliques)for the maxcut POP corresponding
% to the graph contained in the file problem (see Campos, Misener and Parpas, "Partial Lasserre relaxation for sparse Max-Cut", 2020, for more details).
% The function returns 2 text files: 
%          1. clique_problem_r_p_H.txt: File with the maximal 
%             cliques (and subsets if p > 0): each row has 1 in position i if the 
%             ith node is part of the clique or zero ow.
%          2. cliqueaux_problem_r_p_H.txt: File with 1 if the
%             set in position i in the previous file needs to be included as
%             1st order constraint, or 2 if needs to be included as 2nd order.
% Inputs:
% problem: name of the instance 
% r: max size of the maximal clique (or subset) to be added to the
%        relaxation as a second order constraints (r parameter in paper).
% p: the number of subset to include for maximal cliques larger than
%          r (zero if no subsets) as second order constraints (p parameter in paper).
% H: heuristics to follow when including subsets of maximal cliques.
%    (see Campos, Misener and Parpas, "A partial Lasserre relaxation for sparse Max-Cut", 2020, for more details). 
%    H1: select subsets randomly
%    H2: select subsets with maximal absolute value of the weights in the graph.
%    H3: select subsets that repeat the most in the entire set of maximal cliques.
%    H4: select subsets that repeat the least in the entire set of maximal cliques.
%    H5: combination of H2 and H5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BSD 3-Clause License

% Copyright (c) 2020, Juan S. Campos
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    r = 5;
    p = 0;
    rep = "false";
    use_norm = "false";
elseif nargin < 3
    p = 0;
    rep = "false";
    use_norm = "false";
elseif nargin < 4
    rep = "false";
    use_norm = "false";
else
    if strcmp("H1",H)
        rep = "false";
        use_norm = "false";
    elseif strcmp("H2",H)
        rep = "false";
        use_norm = "true";
    elseif strcmp("H3",H)
        rep = "max";
        use_norm = "false";
    elseif strcmp("H4",H)
        rep = "min";
        use_norm = "false";
    elseif strcmp("H5",H)
        rep = "min";
        use_norm = "true";
    end
end


use_rand = 0;
if (strcmp("false",rep) && strcmp("false",use_norm))
    use_rand = 1;
end

if strcmp(rep,"min")
    rep_aux = 0;
elseif strcmp(rep,"max")
    rep_aux = 1;
else
    rep_aux = 2;
end

if strcmp(use_norm,"false")
    use_norm_aux = 0;
elseif  strcmp(use_norm,"true")
    use_norm_aux = 1;
else
    error("the norm parameter is not true or false");
end




max_attempts = 20; % Maximal subsets to consider for maximal cliques

str = sprintf("graphs/%s.txt",problem);
W_aux = dlmread(str);
for i=2:size(W_aux,1)
    W(W_aux(i,1), W_aux(i,2)) = W_aux(i,3);
    W(W_aux(i,2), W_aux(i,1)) = W_aux(i,3);
end
n = size(W,1);

clique_aux = [];


tic
clique = gen_cliques(W);
[lll,sss] = sort(sum(clique,2));
clique = clique(sss,:);
%clique(sum(clique,2) < 2,:) = [];
tot_c = size(clique,1);
fprintf('The sparsity of W is: %f \n',full(nnz(W)/(n^2-n)));
fprintf('The mean of the cliques is: %f \n',full(mean(sum(clique,2))));
fprintf('The median of the cliques is: %f \n',full(median(sum(clique,2))));
fprintf('The max cliques is: %f \n',full(max(sum(clique,2))));
fprintf('The min cliques is: %f \n',full(min(sum(clique,2))));

count_add = 1;


if r == 0
    clique_aux = clique;
    first_second_order = ones(size(clique,1),1);
    
elseif p == 0
    idx_2 = find(sum(clique,2) <= r);
    clique_aux = clique;
    first_second_order = ones(size(clique_aux,1),1);
    first_second_order(idx_2) = 2;
    
else
    
    for i =1:size(clique,1)
        
        c = find(clique(i,:));
        s_c = size(c,2);
        
        
        if s_c<= r
            clique_aux = [clique_aux;clique(i,:)];
            first_second_order(count_add,1) = 2;
            count_add = count_add +1;
            continue;
        else
            clique_aux = [clique_aux;clique(i,:)];
            first_second_order(count_add,1) = 1;
            count_add = count_add +1;
        end
        add = 0;
        
        % total_aux = nchoosek(s_c,r);
        total_aux = max_attempts;
        use_perm = 0;
        max_attempts_aux = max_attempts;
        if total_aux < max_attempts
            use_perm = 1;
            max_attempts_aux = total_aux;
            g = nchoosek(c,r);
        end
        subclique_aux = [];
        info_subclique.norm = [];
        info_subclique.count_rep = [];
        info_subclique.max_eig = [];
        counter = 0;
        rng('default')
        idx_2 = find(first_second_order == 2);
        for j = 1:max_attempts_aux
            
            if use_perm == 0
                c1 = randsample(c, r);
            else
                c1 = g(j,:);
            end
            aux_0_1 = zeros(1,n);
            aux_0_1(c1) = 1;
            
            if size(unique([clique_aux(idx_2,:);aux_0_1],'rows'),1) == size(clique_aux(idx_2,:),1)
                continue;
            end
            if size(unique([subclique_aux;aux_0_1],'rows'),1) == size(subclique_aux,1)
                continue;
            end
            
            subclique_aux = [subclique_aux;aux_0_1];
            counter = counter + 1;
            if use_norm_aux == 1
                info_subclique.norm(counter,1) = norm(W(c1,c1));
            end
           
            if rep_aux < 2
                info_subclique.count_rep(counter,1) = 0;
                for k=1:tot_c
                    c2 = find(clique(k,:));
                    if k == i || size(c2,2)<=r
                        continue;
                    else
                        if sum(ismember(c2,c1)) == r
                            info_subclique.count_rep(counter,1) = info_subclique.count_rep(counter,1) + 1;
                        end
                    end
                end
            end
        end
        size_clique_aux = size(subclique_aux,1);
        p_aux = min(size_clique_aux, p);
        if use_rand == 1
            clique_aux = [clique_aux;subclique_aux(1:p_aux,:)];       
        elseif use_norm_aux == 1 && rep_aux == 2
            [aa,bb] = sort(info_subclique.norm,'descend');
            clique_aux = [clique_aux;subclique_aux(bb(1:p_aux),:)];
        elseif use_norm_aux == 0 && rep_aux == 1
            [aa,bb] = sort(info_subclique.count_rep,'descend');
            clique_aux = [clique_aux;subclique_aux(bb(1:p_aux),:)];
        elseif use_norm_aux == 0 && rep_aux == 0
            [aa,bb] = sort(info_subclique.count_rep,'ascend');
            clique_aux = [clique_aux;subclique_aux(bb(1:p_aux),:)];
        elseif use_norm_aux == 1 && rep_aux == 0
            idx_aux = info_subclique.count_rep == 0;
            if sum(idx_aux) == 0
                idx_aux = logical([1:1:size(subclique_aux,1)]);
            end
            sub_aux = subclique_aux(idx_aux,:);
            info_aux = [info_subclique.norm(idx_aux)];
            [aa,bb] = sort(info_aux,'descend');
            clique_aux = [clique_aux;sub_aux(bb(1:min(p,sum(idx_aux))),:)];
            p_aux = min(p,sum(idx_aux));
        end
        first_second_order(count_add:count_add + p_aux-1,1) = 2;
        count_add = count_add + p_aux;
    end
end

toc

name_aux = sprintf("clique_%s_%d_%d_%s.txt",problem,r,...
    p, H);
fileID = fopen(name_aux,'w');
fprintf(fileID,'%d,%d \n',size(clique_aux,1),n);
fclose(fileID);
dlmwrite(name_aux,full(clique_aux),'-append');
name_aux = sprintf("clique_aux_%s_%d_%d_%s.txt",problem,r,...
    p, H);
fileID = fopen(name_aux,'w');
fprintf(fileID,'%d,%d \n',size(clique_aux,1),1);
fclose(fileID);
dlmwrite(name_aux,full(first_second_order),'-append');

end