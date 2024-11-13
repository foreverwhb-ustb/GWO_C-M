function [bestSolution, bestFitness, iterCurve] = JSO(fitnessFunc, dim, ub, lb, NP, maxEvals)
    % 初始化参数
    H = 5;  % Memory size for F and CR
    memoryCR = 0.8 * ones(H, 1);  % CR memory
    memoryF = 0.5 * ones(H, 1);   % F memory
    archive = [];  % Archive
    archiveSize = round(1.5 * NP);  % Archive max size
    % maxGen = floor(maxEvals / NP);  % Maximum generations
    maxGen = maxEvals;
    pmin = 0.1;  % 最小 p 值
    pmax = 0.2;  % 最大 p 值

    % 初始化种群
    population = initialization(NP, ub, lb, dim);
    fitness = zeros(1, NP);
    for i = 1:NP
        fitness(i) = feval(fitnessFunc, population(i, :));
    end

    % 确定初始最优解
    [bestFitness, idx] = min(fitness);
    bestSolution = population(idx, :);
    iterCurve = zeros(1, maxGen);  % 保存每代的最优适应度

    k = 1;  % Memory index counter
    g = 1;  % Generation counter
    nfes = NP;  % 初始化评估次数

    % 主迭代循环
    while nfes < maxEvals
        SCR = [];  % 保存成功的CR
        SF = [];   % 保存成功的F
        for i = 1:NP
            % 选择随机数 r1，r2，r3
            r1 = selectRandomIndex(NP, i);
            r2 = selectRandomIndex(NP, i, r1);
            r3 = selectRandomIndex(NP, i, r1, r2);
            % r1=randi(NP);
            % r2=randi(NP);
            % r3=randi(NP);

            % Eq.(1): pbest 的动态更新
            p = (pmax - pmin) * nfes / maxEvals + pmin;
            pbestIndex = selectPBest(fitness, p);  % 选择 pbest
            
            % 随机选择记忆中的 F 和 CR
            r = randi(H);
            F = cauchyrnd(memoryF(r), 0.1);  % Cauchy分布
            F = min(1, max(0, F));  % 确保 F 在 0 和 1 之间
            CR = normrnd(memoryCR(r), 0.1);  % 正态分布
            CR = min(1, max(0, CR));  % 确保 CR 在 0 和 1 之间
            
            % Eq.(3) 变异策略 current-to-pbest/1/bin
            mutant = population(i, :) + F * (population(pbestIndex, :) - population(i, :)) + ...
                     F * (population(r1, :) - population(r2, :));
            mutant = ensureBounds(mutant, ub, lb);  % 确保变异体在范围内

            % 二进制交叉
            crossover = crossoverBin(population(i, :), mutant, CR);

            % 计算新个体的适应度
            newFitness = feval(fitnessFunc, crossover);
            nfes = nfes + 1;

            % 选择：如果新个体比旧个体好，更新
            if newFitness < fitness(i)
                population(i, :) = crossover;
                fitness(i) = newFitness;
                SCR = [SCR; CR];  % 保存成功的CR
                SF = [SF; F];     % 保存成功的F
                
                % 如果新个体好于当前最优个体，更新全局最优
                if newFitness < bestFitness
                    bestFitness = newFitness;
                    bestSolution = crossover;
                end
            end
        end

        % 更新记忆中的F和CR
        if ~isempty(SCR)
            memoryF(k) = sum(SF.^2) / sum(SF);
            memoryCR(k) = sum(SCR) / numel(SCR);
            k = mod(k, H) + 1;
        end

        % 记录当前代的最优适应度
        iterCurve(g) = bestFitness;

        % 种群大小线性缩减（LPSR）
        if nfes > maxEvals * 0.9
            NP = round(NP * (maxEvals - nfes) / (0.1 * maxEvals));
            population = population(1:NP, :);
            fitness = fitness(1:NP);
        end

        g = g + 1;
    end
end


% 选择随机数r，确保不同于i, r1, r2
function r = selectRandomIndex(NP, idx1, idx2, idx3)
    r = randi(NP);
end

% 选择pbest
function pbestIndex = selectPBest(fitness, p)
    [~, sortedIndices] = sort(fitness);
    numPBest = max(round(p * numel(fitness)), 1);
    pbestIndex = sortedIndices(randi(numPBest));
end

% 确保变异体在边界内
function mutant = ensureBounds(mutant, ub, lb)
    mutant = max(lb, min(ub, mutant));
end

% 二进制交叉操作
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

% 生成柯西分布随机数
function r = cauchyrnd(mu, sigma)
    r = mu + sigma * tan(pi * (rand - 0.5));
end
