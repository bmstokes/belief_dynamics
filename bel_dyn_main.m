% BEL_DYN_MAIN performs a series of opinion dynamics simulations according 
% to the model Stokes et al, in "Extremism and segregation emerge naturally
% through collective opinion dynamics in a novel agent-based model".
% Calls: bel_dyn_onesimulation.
%% 
fid = fopen('parameters.txt'); % read parameters
params_cell = cell(11,1);
params = zeros(10,1);
tline = fgetl(fid);
i = 1;
while ischar(tline)
    params_cell{i,1}=tline;
    i = i+1;
    tline = fgetl(fid);
end
fclose(fid);
for i = 1:10
    p = char(params_cell{i,1});
    p(1:find(p=='='))=[];
    params(i) = str2double(p);
end
D = params(1);
N = params(2);
rho_min = params(3);
rho_max = params(4);
del_rho = params(5);
mu = params(6);
beta = params(7); % always 0.5 in the paper, but can be anything >= 0
alpha = params(8);
T = params(9);
S = params(10);
S_input = S;
%
clearvars fid params_cell params tline i p ans
%% Main
betastring = sprintf('%0.4f',beta);
betastring(betastring=='.') = [];
alphastring = sprintf('%0.4f',alpha);
alphastring(alphastring=='.') = [];
for rho = rho_min : del_rho : rho_max 
    init_configs_used = zeros(1000,1); % to indicate which of the pre-generated initial configs have been used
    rhostring = sprintf('%0.4f',rho);
    rhostring(rhostring=='.') = [];
    if alpha == 0
        tempstr = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s',D,N,rhostring,mu,betastring);
    else
        tempstr = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s',D,N,rhostring,mu,betastring,alphastring);
    end
    if ~exist(tempstr,'dir') 
        mkdir (tempstr);
        fails = 0;    
        failed_init_configs = []; % to store init configs which fail to lead to steady state
        failed_init_indices = []; % to store indices of init configs which fail to lead to steady state
        s = 0;
    else
        if alpha == 0
            tempstr = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s/D=%d_N=%d_rho=%s_mu=%d_beta=%s.mat',...
                D,N,rhostring,mu,betastring,D,N,rhostring,mu,betastring);
        else
            tempstr = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s/D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s.mat',...
                D,N,rhostring,mu,betastring,alphastring,D,N,rhostring,mu,betastring,alphastring);
        end
        tempstruct = matfile(tempstr);
        fails = tempstruct.fails;
        failed_init_configs = tempstruct.failed_init_configs; 
        failed_init_indices = tempstruct.failed_init_indices;
        S =  S_input + tempstruct.S;  
        s = tempstruct.S; 
    end
    while s < S
        bel_dyn_onesimulation;
        s = s + 1;
    end
    % save global mat file:
    if alpha == 0
        tempdir = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s',D,N,rhostring,mu,betastring);
        tempfn = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s.mat',D,N,rhostring,mu,betastring);
    else
        tempdir = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s',D,N,rhostring,mu,betastring,alphastring);
        tempfn = sprintf('D=%d_N=%d_rho=%s_mu=%d_beta=%s_alpha=%s.mat',D,N,rhostring,mu,betastring,alphastring);
    end
    save(fullfile(tempdir,tempfn),...
        'alignment_tol','beta','alpha','D','mu','N','rho','S','fails',...
        'failed_init_configs','failed_init_indices','convergence_steps_threshold','T'); 
    clear tempdir tempfn
end
%%
clearvars acceleration i m N_groups_old N_individuals_per_group_old Q rho rhostring 
clearvars steps_with_fixed_groups N_groups N_individuals_per_group timesteps
clearvars t_conv v v0_mean vdist vfinal_mean tempstruct tempstr computation_time timer alphastring betastring