function [v_new,acceleration,Q_new] = bel_dyn_onestep(D,N,rho,mu,beta,alpha,v,i) 
% In a vector, the elements are arranged in the following manner: 
% the first D elements correspond to the 1st individual, then 
% elements D+1 to 2D correspond to the 2nd individual, and so on. 
% In a matrix, each row corresponds to a different individual.
%%
alignment = zeros(N,N); % pairwise alignments
% MAKE THE (j,k) ELEMENT OF THIS sum_{tau=t-mu+1}^t of ||v_k(tau)-v_j(tau)||^2
%
affinity = zeros(N,N); % affinity matrix (a_jk)
%
ca = zeros(N,N); % c_jk * a_jk
cav = zeros(N,N*D); % c_jk * a_jk * (v_k - v_j)
%
acceleration = zeros(N,D); 
% accelerations of the N individuals, each having D dimensions
%
Q_new = zeros(N,1); % number of agents that j listens to
%%
for j = 1:N % for each of the N individuals
    %% calculate rho_r
    r = sqrt(sum( v(j, (min(i,mu)-1)*D+1 : (min(i,mu)-1)*D+D).^2 )); % distance of agent j to origin
    rho_r = rho + (1-rho)*(1-exp(-alpha*r));
    %%
    for k = 1:N % compute the distances between j and each of the others
        temp = zeros(1,mu); 
        % row vector with space to input the opinion distance for each timestep once it is calculated
        % defining it like this should mean we can use the same code for both i<=mu and i>mu
        for m = 1:length(temp)   % for each of the last i/mu timesteps worth of opinion blocks
            temptemp = 0;
            for n =  (m-1)*D+1:(m-1)*D+D  % for each of the D opinions in each block
                temptemp = temptemp + (v(k,n) - v(j,n)).^2; 
                % gives the euclidean distance between j and k for each of the timesteps
            end     
            temp(m) = temptemp;    
            % for each (j,k) temp should now be a row vector storing the last 
            % i/mu euclidean opinion distances between agents j and k
        end 
        alignment(j,k) = sum(temp);    % sums the last i/mu opinion distances for each (j,k)
    end
    %%
    for k = 1:N
        zeroth_index_of_k = (k-1)*D;
        %
        if beta > 0
            affinity(j,k) = 1 / ((1 + alignment(j,k) ).^beta);
        else
            affinity(j,k) = 1;
        end
        %
        if affinity(j,k) > rho_r
            ca(j,k) = affinity(j,k);
            Q_new(j) = Q_new(j) + 1;
        else
            ca(j,k) = 0;
        end
        %
        for m = 1:D % run through the D dimensions  
            if i <= mu
                cav(j,zeroth_index_of_k + m) = ca(j,k) * (v(k,(i-1)*D+m) - v(j,(i-1)*D+m));
            else
                cav(j,zeroth_index_of_k + m) = ca(j,k) * (v(k,end-D+m) - v(j,end-D+m));
            end
        end
    end
    %%
    cavsum = zeros(1,D);  
    % these get reassigned to zeros every time we start a new j
    %
    for k = 1:N
        zeroth_index_of_k = (k-1)*D;
        for m = 1:D % run through the D dimensions
            cavsum(m) = cavsum(m) + cav(j,zeroth_index_of_k + m);
            % cavsum(m) is the sum of affinity values of agents with j that  
            % are above the threshold of rho, in the mth dimension; 
            %
            acceleration(j,m) = cavsum(m) / ( Q_new(j) ); % acceleration of individual j
        end
    end
end
%%
if i <= mu
    v_new = v(:,(i-1)*D+1:(i-1)*D+D) + acceleration;
else
    v_new = v(:,(mu-1)*D+1: (mu-1)*D+D) + acceleration;
end
