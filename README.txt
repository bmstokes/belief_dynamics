v1.0.0

To perform opinion dynamics simulations according to the model in 
"Extremism and segregation emerge naturally through collective opinion dynamics in a novel agent-based model", 
do:

1. Set parameters.txt, where:
D = dimensionality of opinion space
N = number of agents
rho_min = minimum value of the threshold parameter
rho_max = maximum value of the threshold parameter
del_rho = increment in the threshold parameter
mu = memory capacity
beta = decay of influence (the power to which the 
       (1+sum(w*(v_j-v_i))) factor is raised in the expression for pairwise affinity; 
       equals 0.5 in the paper, but can take any non-negative value)
alpha = reinforcement rate
T = maximum number of time-steps to run each simulation for
S = the number of different initial conditions to use per parameter set (D,N,rho,mu,beta,alpha)

2. Run bel_dyn_main.m, which performs S*(1+(rho_max-rho_min)/del_rho) simulations.

The .mat files, e.g. bel_dyn_init_D=1_N=100.mat, each contains 1000 pre-generated random initial conditions, 
to be used by the simulations. 

Further information can be accessed in MATLAB by calling, e.g., "help bel_dyn_main".