% BEL_DYN_ONESIMULATION runs one simulation.
% Calls: bel_dyn_onestep; bel_dyn_conv_test.
% Called by: bel_dyn_main.
%% 
convergence_steps_threshold = 100;
alignment_tol = 1e-6;
% This is what pairwise alignment, i.e. ||vk-vj||^2, is compared to 
% in order to determine whether j&k have reached consensus.
%% Initialise opinions and variables
tempstruct = matfile(sprintf('bel_dyn_init_D=%d_N=%d.mat',D,N));
v0 = tempstruct.v0;
if (s >= 0) && (s < 1000) && (init_configs_used(s+1) == 0)
    v0 = v0(:,s*D+1:(s+1)*D); % read v0 from pre-generated initial configurations
    init_configs_used(s+1) = 1;
elseif s >= 1000
    return
else
    % we would get here if a simulation with the 
    % (s+1)th init config has been done and failed
    failed_init_configs = [failed_init_configs,v0(:,s*D+1:(s+1)*D)];                    
    failed_init_indices = [failed_init_indices;s]; % s between 0 and 999
    return
end
clearvars tempstruct
fprintf('D = %d, N = %d, rho = %g, mu = %d, beta = %g, alpha = %g, simulation %d.\n',...
    D,N,rho,mu,beta,alpha,s+1)
%
timer = tic; % start timer
%
N_groups_old = N; % how many groups have formed
N_individuals_per_group_old = ones(N,1); % number of individuals in each group
t_conv = Inf; % t at which convergence occurs
%
v = [v0, zeros(N,T*D)];% to store all the v; D columns for each time
Q = zeros(N,T); % to store the number of agents each is listenting to at each timestep
v0_mean = zeros(1,D); % to store the initial mean opinion
vfinal_mean = v0_mean; % to store the final mean opinion
for m = 1:D
    v0_mean(m) = sum(v0(:,m))/N;
end
%% Run simulation from t=0 to t=T
steps_with_fixed_groups = 0;
for i = 1:T
    % Advance by one step:
    if i <= mu
        [v(:,i*D+1:i*D+D),acceleration,Q(:,i)] = ...
            bel_dyn_onestep(D,N,rho,mu,beta,alpha,v(:,1:(mu-1)*D+D),i);
    else
        [v(:,i*D+1:i*D+D),acceleration,Q(:,i)] = ...
            bel_dyn_onestep(D,N,rho,mu,beta,alpha,v(:,(i-mu)*D+1:(i-1)*D+D),i);
    end
    % Test for convergence:
    [N_groups, N_individuals_per_group] = bel_dyn_conv_test(D,N,v(:,(i-1)*D+1:(i-1)*D+D),alignment_tol);
    if (N_groups ~= N_groups_old) || ...
            ( (N_groups == N_groups_old) && ...
            (sum(N_individuals_per_group == N_individuals_per_group_old) < N_groups) ) || ...
            (max(max(abs(acceleration)) ) > alignment_tol)
        % If the number of groups has changed, 
        % or if the number of groups has not changed but the grouping has changed,
        % or if any component of any acceleration is larger than some
        % threshold.   
        %
        steps_with_fixed_groups = 0;
    else
        steps_with_fixed_groups = steps_with_fixed_groups + 1;
        % This is the number of steps that have elapsed during which
        % grouping has remained fixed.
        %
    end
    N_groups_old = N_groups;
    N_individuals_per_group_old = N_individuals_per_group;
    % If convergence reached, set t_conv break out of for-loop:
    if (steps_with_fixed_groups >= convergence_steps_threshold) ...
            && ( max(max(abs( v(:,i*D+1:i*D+D) - v(:,(i-convergence_steps_threshold)*D+1:(i-convergence_steps_threshold)*D+D) ))) ...
                <= alignment_tol ) 
        % the second condition: no opinion has moved by more than
        % alignment_tol in any dimension for at least
        % convergence_steps_threshold steps
        t_conv = i - convergence_steps_threshold;
        break
    end
end
%
if steps_with_fixed_groups < convergence_steps_threshold
    t_conv = Inf;
elseif i == convergence_steps_threshold % when no changes in grouping ever occurs, since t=0
    t_conv = 0; 
end
v(:,(i+1)*D+1:end) = [];
for m = 1:D
    vfinal_mean(m) = sum(v(:,end-D+m))/N;
end
vdist = sqrt((vfinal_mean - v0_mean) * transpose(vfinal_mean - v0_mean)); % distance between v0_mean and vfinal_mean
timesteps = i;
computation_time = toc(timer);
%% Save data
if t_conv < Inf
    betastring = sprintf('%0.4f',beta);
    betastring(betastring=='.') = [];
    alphastring = sprintf('%0.4f',alpha);
    alphastring(alphastring=='.') = [];
    %
    rhostring = sprintf('%0.4f',rho);
    rhostring(rhostring=='.') = [];
    %
    if alpha == 0
        save(fullfile(sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s',D,N,rhostring,mu,betastring),...
            (sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_%05d.mat',D,N,rhostring,mu,betastring,s))),...
                 'computation_time','N_groups','N_individuals_per_group','t_conv','timesteps','v',...
                 'v0_mean','vfinal_mean','vdist'); 
    else
        save(fullfile(sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s',D,N,rhostring,mu,betastring,alphastring),...
            (sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s_%05d.mat',D,N,rhostring,mu,betastring,alphastring,s))),...
                 'computation_time','N_groups','N_individuals_per_group','t_conv','timesteps','v',...
                 'v0_mean','vfinal_mean','vdist'); 
    end
else
    if alpha == 0
        save(fullfile(sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s',D,N,rhostring,mu,betastring),...
            (sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_%05d.mat',D,N,rhostring,mu,betastring,s))),...
                 'computation_time','t_conv','timesteps','v','v0_mean'); 
    else
        save(fullfile(sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s',D,N,rhostring,mu,betastring,alphastring),...
            (sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s_%05d.mat',D,N,rhostring,mu,betastring,alphastring,s))),...
                 'computation_time','N_groups','N_individuals_per_group','t_conv','timesteps','v',...
                 'computation_time','t_conv','timesteps','v','v0_mean'); 
    end
    s = s-1;
    fails = fails + 1;
end