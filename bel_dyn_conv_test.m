% BEL_DYN_CONV_TEST tests whether or not a system has converged to a steady
% state according to some prescribed criteria.
% Called by: bel_dyn_onesimulation.
% Inputs: D (dimensionality of opinion space); N (number of agents); v
% (opinions at current time); tol (tolerance for convergence).
function [N_groups,N_individuals_per_group] = bel_dyn_conv_test(D,N,v,tol)
% N_groups is how many groups have formed
% N_individuals_per_group is a vector of how many agents there are per group
%%
groups = eye(N);
% In the end, the non-zero elements of the jth row will all be 1,
% and will indicate the indices of agents who are in the same group as j.
%
unaccountedfor = 1:N; % agents which have not found a group
%%
j = 0; % index for groups
jj = 0; % index for v; basically the same thing as j except we never do jj=jj-1
NN = N; % we never do NN=NN-1
while j < N
    j = j+1;
    jj = jj+1;
    if ~isempty( find(unaccountedfor==jj,1) ) % if jj is still unaccounted for
        for k = 1:NN
            if ( k ~= jj ) && ( ~isempty( find(unaccountedfor==k,1) ) ) 
                % if k~=jj and k is still unaccounted for
                %
                alignment = 0;
                for m = 1:D
                    alignment = alignment + (v(k,m) - v(jj,m)).^2;
                end % now we have the alignment between agents jj&k
                if alignment < tol % then we consider jj&k to be flocking
                    groups(j,k) = 1;
                    unaccountedfor(unaccountedfor==k) = []; % delete k from 'unaccountedfor'
                end
            end
        end
        if length(find(groups(j,:)>0)) > 1 % if this row doesn't look like 0 0 0 ... 0 0 1 0 0 ... 0 0 0
            unaccountedfor(unaccountedfor==jj) = []; % delete jj from 'unaccountedfor'
        end
    else % if jj is already in a group with someone else
        groups(j,:) = [];
        j = j-1;
        N = N-1;
    end
    if isempty(unaccountedfor) % i.e. if every individual has been accounted for
        if j < N
            groups(j+1:end,:) = []; % delete the remaining rows of groups
        end
        break % exit the big for-loop
    end
end
%%
% Now j is the number of rows still remaining in groups,
% so it is exactly N_groups
N_groups = j; % this is the number of groups that have formed
N_individuals_per_group = ones(N_groups,1);
for i = 1:N_groups
    N_individuals_per_group(i) = sum(groups(i,:));
end