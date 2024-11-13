function [bestSolution, bestFitness, iterCurve] = DISH(fitnessFunc, dim, ub, lb, NP_init, maxEvals)
    % 初始化参数
    NP = NP_init;  % 初始种群大小
    G = 0;  % 代数
    maxGen = maxEvals;  % 最大代数
    pmin = 2 / NP_init;  % p 最小值
    pmax = 0.2;  % p 最大值
    A = [];  % 初始化外部存档
    k = 1;  % Memory index counter
    archiveSize = round(1.5 * NP);  % 最大存档大小
    iterCurve = zeros(1, maxGen);  % 保存每代的最优适应度

    NP_final = NP*0.2;
    H=5;
    
    % 初始化种群和适应度
    population = initialization(NP, ub, lb, dim);
    fitness = zeros(1, NP);
    for i = 1:NP
        fitness(i) = feval(fitnessFunc, population(i, :));
    end
    
    [bestFitness, idx] = min(fitness);
    bestSolution = population(idx, :);
    
    % 初始化F和CR的记忆
    memoryCR = 0.8 * ones(H, 1);  % CR记忆
    memoryF = 0.5 * ones(H, 1);   % F记忆
    
    P_new = [];  % 新种群
    G_max = maxGen;  % 最大代数
    SF = []; SCR = [];  % 用于存储成功的 F 和 CR
    
    % 主循环
    while G < maxGen
        G = G + 1;
        FES = G * NP;  % 当前的评估次数
        CR = zeros(1, NP);
        F = zeros(1, NP);
        
        for i = 1:NP
            % 随机选择记忆中的 F 和 CR
            r = randi(H);
            F(i) = cauchyrnd(memoryF(r), 0.1);  % Cauchy分布
            F(i) = max(0.1, min(F(i), 1));  % 限制 F 在 [0.1, 1]
            
            CR(i) = normrnd(memoryCR(r), 0.1);  % 正态分布
            CR(i) = max(0, min(CR(i), 1));  % 确保 CR 在 [0, 1] 之间
            
            % 动态选择 pbest
            p = pmin + (pmax - pmin) * (FES / maxEvals);
            pbestIndex = selectPBest(fitness, p);
            
            % 选择随机数 r1, r2
            r1 = selectRandomIndex(NP, i);
            r2 = selectRandomIndex(NP, i, r1);
            
            % 变异操作 (Eq. 15)
            mutant = population(i, :) + F(i) * (population(pbestIndex, :) - population(i, :)) + ...
                     F(i) * (population(r1, :) - population(r2, :));
            mutant = ensureBounds(mutant, ub, lb);  % 确保变异体在边界内
            
            % 二进制交叉 (Eq. 6)
            crossover = crossoverBin(population(i, :), mutant, CR(i));
            
            % 计算新个体的适应度
            newFitness = feval(fitnessFunc, crossover);
            
            % 选择操作
            if newFitness < fitness(i)
                P_new = [P_new; crossover];
                SF = [SF; F(i)];
                SCR = [SCR; CR(i)];
                
                % 更新存档
                if size(A, 1) < archiveSize
                    A = [A; crossover];
                else
                    A(randi(size(A, 1)), :) = crossover;
                end
                
                fitness(i) = newFitness;  % 更新适应度
            else
                P_new = [P_new; population(i, :)];
            end
        end
        
        % 更新记忆
        if ~isempty(SF) && ~isempty(SCR)
            memoryF(k) = weightedLehmerMean(SF);
            memoryCR(k) = weightedLehmerMean(SCR);
            k = k + 1;
            if k > H
                k = 1;  % 环形更新记忆
            end
        end
        
        % 种群大小缩减 (Eq. 13)
        NP_new = round(NP_init - (FES / maxEvals) * (NP_init - NP_final));
        if NP_new < NP
            % 根据适应度排序，删除最差的个体
            [~, sortedIndices] = sort(fitness);
            worstIndices = sortedIndices(NP_new+1:end);
            population(worstIndices, :) = [];
            fitness(worstIndices) = [];
            NP = NP_new;
        end
        
        % 更新种群
        population = P_new;
        P_new = [];
        
        % 记录每代的最优适应度
        [bestFitness, idx] = min(fitness);
        bestSolution = population(idx, :);
        iterCurve(G) = bestFitness;
        
        % 打印当前迭代状态（可选）
        % disp(['Generation ', num2str(G), ' Best Fitness: ', num2str(bestFitness)]);
    end
end

% 初始化种群
function pop = initialization(NP, ub, lb, dim)
    pop = lb + (ub - lb) .* rand(NP, dim);
end

% 选择 pbest
function pbestIndex = selectPBest(fitness, p)
    [~, sortedIndices] = sort(fitness);
    numPBest = max(round(p * numel(fitness)), 1);
    pbestIndex = sortedIndices(randi(numPBest));
end

% 确保变异体在边界内
function mutant = ensureBounds(mutant, ub, lb)
    mutant = max(lb, min(ub, mutant));
end

% 二进制交叉
function crossover = crossoverBin(target, mutant, CR)
    dim = numel(target);
    crossover = target;
    jrand = randi(dim);
    for j = 1:dim
        if rand <= CR || j == jrand
            crossover(j) = mutant(j);
        end
    end
end

% 选择随机数r，确保不同于i, r1, r2
function r = selectRandomIndex(NP, idx1, idx2)
    r = randi(NP);
end

% 生成柯西分布随机数
function r = cauchyrnd(mu, sigma)
    r = mu + sigma * tan(pi * (rand - 0.5));
end

% 计算加权的Lehmer均值
function meanWL = weightedLehmerMean(S)
    w = abs(S - mean(S));  % 权重
    meanWL = sum(w .* S.^2) / sum(w .* S);
end
