% LSHADE算法主体
function [bestSolution, bestFitness, iterCurve] = LSHADE(fitnessFunc, constraints, lower, upper, pop, maxIterations, dim)
    NP = pop;  % 种群大小
    F = 0.5;   % 缩放因子
    CR = 0.5;  % 交叉概率
    pbest_rate = 0.05;  % pbest选择概率
    memorySize = 10;  % 成功历史记忆大小

    population = initialization(lower, upper, dim, NP);  % 初始化种群
    fitness = zeros(NP, 1);  % 初始化适应度值
    conv = zeros(NP, 1);  % 初始化约束

    for i = 1:NP
        fitness(i) = feval(fitnessFunc, population(i, :));  % 计算初始种群的适应度值
        conv(i) = feval(constraints, population(i, :));  % 计算初始种群的约束违反程度
    end

    [~, sortedIndices] = sort(fitness);
    bestSolution = population(sortedIndices(1), :);  % 全局最优解
    bestFitness = fitness(sortedIndices(1));  % 全局最优适应度值
    iterCurve = zeros(1, maxIterations);  % 迭代曲线

    memoryCR = zeros(memorySize, 1);  % 成功历史CR
    memoryF = zeros(memorySize, 1);  % 成功历史F
    memoryIndex = 1;

    Ninit = NP;
    Nmin = 4;
    NFE = 0;
    max_NFE = 10000 * dim;

    for iter = 1:maxIterations
        % 线性减少种群大小
        NP = max(Nmin, Ninit - floor((Ninit-Nmin)*iter/maxIterations));
        popReduction(population, NP, dim);

        [F, CR] = updateParameters(memoryF, memoryCR, memoryIndex);  % 更新缩放因子和交叉概率

        for i = 1:NP
            pbestIndex = selectPBest(fitness, pbest_rate);  % 选择pbest

            r1 = randi([1, NP]);  % 选择邻居r1
            while r1 == i
                r1 = randi([1, NP]);
            end

            r2 = randi([1, NP]);  % 选择邻居r2
            while r2 == i || r2 == r1
                r2 = randi([1, NP]);
            end

            r3 = randi([1, NP]);  % 选择邻居r3
            while r3 == i || r3 == r1 || r3 == r2
                r3 = randi([1, NP]);
            end

            mutant = population(r1, :) + F * (population(pbestIndex, :) - population(r3, :));  % 生成新个体
            crossover = crossoverBin(population(i, :), mutant, CR);  % 交叉操作

            newFitness = feval(fitnessFunc, crossover);  % 计算新个体的适应度值
            newConv = feval(constraints, crossover);  % 计算新个体的约束违反程度

            if (conv(i) > newConv || (conv(i) == 0 && newConv == 0 && fitness(i) > newFitness))  % 更新个体和种群
                population(i, :) = crossover;
                fitness(i) = newFitness;
                conv(i) = newConv;

                memoryF(memoryIndex) = F;  % 更新成功历史CR和F
                memoryCR(memoryIndex) = CR;
                memoryIndex = mod(memoryIndex, memorySize) + 1;
            end
        end

        [~, sortedIndices] = sort(fitness);
        bestSolution = population(sortedIndices(1), :);  % 更新全局最优解
        bestFitness = fitness(sortedIndices(1));  % 更新全局最优适应度值
        iterCurve(1,iter) = bestFitness;  % 记录迭代曲线
    end
end

% 种群缩减函数
function popReduction(population, NP, dim)
    population = population(1:NP, :);
end

% pbest选择函数
function pbestIndex = selectPBest(fitness, pbest_rate)
    [~, sortedIndices] = sort(fitness);
    numPBest = round(pbest_rate * numel(fitness));
    pbestIndex = sortedIndices(1:numPBest);
end

% 缩放因子和交叉概率更新函数
function [F, CR] = updateParameters(memoryF, memoryCR, memoryIndex)
    F = mean(memoryF);
    CR = mean(memoryCR);
end

% 二进制交叉
function crossover = crossoverBin(target, mutant, CR)
    dim = numel(target);
    crossover = target;
    crossover(rand(dim, 1) < CR) = mutant(rand(dim, 1) < CR);
end

% 种群初始化函数
function population = initialization(lower, upper, dim, NP)
    population = lower + (upper - lower) .* rand(NP, dim);
end