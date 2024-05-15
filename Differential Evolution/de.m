%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Differential Evolution (DE) in MATLAB (URL: https://yarpiz.com/231/ypea107-differential-evolution), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

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

%% Problem Definition

% CostFunction = @(x) Sphere(x);    % Cost Function

dim = 2;            % Number of Decision Variables

VarSize = [1 dim];   % Decision Variables Matrix Size

VarMin = -10;          % Lower Bound of Decision Variables
VarMax = 10;          % Upper Bound of Decision Variables

%% DE Parameters

MaxIt = 250;      % Maximum Number of Iterations

nPop = 250;        % Population Size

beta_min = 0.2;   % Lower Bound of Scaling Factor
beta_max = 0.8;   % Upper Bound of Scaling Factor

pCR = 0.2;        % Crossover Probability

%% Initialization

% empty_individual.Position = [];
% empty_individual.Cost = [];
% 
% BestSol.Cost = inf;
% 
% pop = repmat(empty_individual, nPop, 1);
% 
% for i = 1:nPop
% 
%     pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
% 
%     pop(i).Cost = dec.objective_function(pop(i).Position);
% 
%     if pop(i).Cost<BestSol.Cost
%         BestSol = pop(i);
%     end
% 
% end
% 
BestCost = zeros(MaxIt, 1);

[pop,BestSol] = dec.generate_points(nPop,boundaries,seed);

%% DE Main Loop

for it = 1:MaxIt
    
    for i = 1:nPop
        
        x = pop(i).Position;
        
        A = randperm(nPop);
        
        A(A == i) = [];
        
        a = A(1);
        b = A(2);
        c = A(3);
        
        % Mutation
        %beta = unifrnd(beta_min, beta_max);
        % beta = unifrnd(beta_min, beta_max, VarSize);
        beta = repmat(scaling_factor,VarSize);
        y = pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y = max(y, boundaries(:,1)');
        y = min(y, boundaries(:,2)');
		
        % Crossover
        z = zeros(size(x));
        j0 = randi([1 numel(x)]);
        for j = 1:numel(x)
            if j == j0 || rand <= crossover_prob
                z(j) = y(j);
            else
                z(j) = x(j);
            end
        end
        
        NewSol.Position = z;
        NewSol.Cost = dec.objective_function(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i) = NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol = pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it) = BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
% BestSol
%% Show Results

figure;
%plot(BestCost);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
pause(5)
close;
%%
clc;

% Accessing all the 'Position' fields using arrayfun
populationcell = arrayfun(@(x) x.Position, pop, 'UniformOutput', false);

population = cell2mat(populationcell)


