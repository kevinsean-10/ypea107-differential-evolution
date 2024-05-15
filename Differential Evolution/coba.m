clear;clc;close all;

dim = 2;
boundaries = repmat([-10, 10], dim, 1);
population_size = 250;
max_iter = 250;
scaling_factor = 0.8;
crossover_prob = 0.2;
seed = 0;

rng(seed)


dec = DE_class(boundaries,max_iter,population_size, ...
                scaling_factor,crossover_prob,seed);

clc;
[p1,b1] = dec.generate_points(population_size,boundaries,seed);
[population,BestSol] = dec.DE(p1,b1,boundaries,max_iter,scaling_factor,crossover_prob,true)
