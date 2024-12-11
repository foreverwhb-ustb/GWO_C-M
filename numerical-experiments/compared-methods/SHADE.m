% % SHADE算法主体
% function [bestSolution, bestFitness, iterCurve] = SHADE(fitnessFunc, dim, ub, lb, pop, maxIterations, constraintsFunc)
%     NP = pop;  % 种群大小
%     F = 0.5;   % 缩放因子
%     CR = 0.5;  % 交叉概率
%     pbest_rate = 0.05;  % pbest选择概率
%     memorySize = 10;  % 成功历史记忆大小
% 
%     population = initialization(NP, ub, lb, dim);  % 初始化种群
%     fitness = zeros(1, NP);  % 初始化适应度值
% 
%     % 计算初始种群的适应度并处理约束
%     for i = 1:NP
%         fitness(i) = handle_constraints(fitnessFunc(population(i, :)), constraintsFunc(population(i, :)));
%     end
% 
%     [~, sortedIndices] = sort(fitness);
%     bestSolution = population(sortedIndices(1), :);  % 全局最优解
%     bestFitness = fitness(sortedIndices(1));  % 全局最优适应度值
%     iterCurve = zeros(1, maxIterations);  % 迭代曲线
% 
%     memoryCR = zeros(memorySize, 1);  % 成功历史CR
%     memoryF = zeros(memorySize, 1);  % 成功历史F
%     memoryIndex = 1;
% 
%     for iter = 1:maxIterations
%         [F, CR] = updateParameters(memoryF, memoryCR, memoryIndex);  % 更新缩放因子和交叉概率
% 
%         for i = 1:NP
%             pbestIndex = selectPBest(fitness, pbest_rate);  % 选择pbest
% 
%             r1 = randi(NP);  % 选择邻居r1
%             while r1 == i
%                 r1 = randi(NP);
%             end
% 
%             r2 = randi(NP);  % 选择邻居r2
%             while r2 == i || r2 == r1
%                 r2 = randi(NP);
%             end
% 
%             r3 = randi(NP);  % 选择邻居r3
%             while r3 == i || r3 == r1 || r3 == r2
%                 r3 = randi(NP);
%             end
% 
%             % 差分变异操作
%             mutant = population(r1, :) + F * (population(r2, :) - population(r3, :));
%             mutant = BoundrayCheck(mutant, ub, lb, dim);  % 边界检查
% 
%             % 二进制交叉操作
%             crossover = crossoverBin(population(i, :), mutant, CR);
% 
%             % 计算新个体的适应度并处理约束
%             newFitness = handle_constraints(fitnessFunc(crossover), constraintsFunc(crossover));
% 
%             % 检查新个体是否优于当前个体，并更新个体和适应度
%             if newFitness < fitness(i)
%                 population(i, :) = crossover;
%                 fitness(i) = newFitness;
% 
%                 % 更新成功的缩放因子和交叉概率历史
%                 memoryF(memoryIndex) = F;
%                 memoryCR(memoryIndex) = CR;
%                 memoryIndex = mod(memoryIndex, memorySize) + 1;
%             end
% 
%             % 检查新个体是否是全局最优解，并更新全局最优解
%             if fitness(i) < bestFitness
%                 bestSolution = population(i, :);
%                 bestFitness = fitness(i);
%             end
%         end
% 
%         % 记录当前最优解的适应度值到迭代曲线中
%         iterCurve(iter) = bestFitness;
% 
%         % 显示当前迭代和最优适应度
%         disp(['Iteration: ', num2str(iter), ' Best Fitness: ', num2str(bestFitness)]);
%     end
% end
% 
% % pbest选择函数
% function pbestIndex = selectPBest(fitness, pbest_rate)
%     [~, sortedIndices] = sort(fitness);
%     numPBest = round(pbest_rate * numel(fitness));
%     pbestIndex = sortedIndices(1:numPBest);
% end
% 
% % 缩放因子和交叉概率更新函数
% function [F, CR] = updateParameters(memoryF, memoryCR, memoryIndex)
%     F = mean(memoryF);
%     CR = mean(memoryCR);
% end
% 
% % 二进制交叉操作
% function crossover = crossoverBin(target, mutant, CR)
%     dim = numel(target);
%     crossover = target;
%     crossover(rand(dim, 1) < CR) = mutant(rand(dim, 1) < CR);
% end
% 
% % 边界检查函数：确保解在上下界内
% function x = BoundrayCheck(x, ub, lb, dim)
%     x = min(max(x, lb), ub);  % 使用min和max确保解在边界内
% end
% 
% % % 约束处理函数：通过惩罚函数处理违反约束的解
% % function penalized_fitness = handle_constraints(fitness, constraints)
% %     penalty = 0;
% %     for i = 1:length(constraints)
% %         if constraints(i) > 0  % 如果违反约束
% %             penalty = penalty + constraints(i);  % 增加违反约束的惩罚
% %         end
% %     end
% %     penalized_fitness = fitness + 10000 * penalty;  % 添加惩罚系数
% % end




% SHADE算法主体
function [bestSolution, bestFitness, iterCurve] = SHADE(fitnessFunc, dim, ub, lb, pop, maxIterations)
    NP = pop;  % 种群大小
    F = 0.5;   % 缩放因子
    CR = 0.5;  % 交叉概率
    pbest_rate = 0.05;  % pbest选择概率
    memorySize = 10;  % 成功历史记忆大小

    population = initialization(NP, ub, lb, dim);  % 初始化种群
    fitness = zeros(1, NP);  % 初始化适应度值

    for i = 1:NP
        fitness(i) = feval(fitnessFunc, population(i, :));  % 计算初始种群的适应度值
    end

    [~, sortedIndices] = sort(fitness);
    bestSolution = population(sortedIndices(1), :);  % 全局最优解
    bestFitness = fitness(sortedIndices(1));  % 全局最优适应度值
    iterCurve = zeros(1, maxIterations);  % 迭代曲线

    memoryCR = zeros(memorySize, 1);  % 成功历史CR
    memoryF = zeros(memorySize, 1);  % 成功历史F
    memoryIndex = 1;

    for iter = 1:maxIterations
        [F, CR] = updateParameters(memoryF, memoryCR, memoryIndex);  % 更新缩放因子和交叉概率

        for i = 1:NP
            pbestIndex = selectPBest(fitness, pbest_rate);  % 选择pbest

            r1 = randi(NP);  % 选择邻居r1
            while r1 == i
                r1 = randi(NP);
            end

            r2 = randi(NP);  % 选择邻居r2
            while r2 == i || r2 == r1
                r2 = randi(NP);
            end

            r3 = randi(NP);  % 选择邻居r3
            while r3 == i || r3 == r1 || r3 == r2
                r3 = randi(NP);
            end

            mutant = population(r1, :) + F * (population(r2, :) - population(r3, :));  % 生成新个体
            crossover = crossoverBin(population(i, :), mutant, CR);  % 交叉操作

            newFitness = feval(fitnessFunc, crossover);  % 计算新个体的适应度值

            if newFitness < fitness(i)  % 更新个体和种群
                population(i, :) = crossover;
                fitness(i) = newFitness;

                memoryF(memoryIndex) = F;  % 更新成功历史CR和F
                memoryCR(memoryIndex) = CR;
                memoryIndex = mod(memoryIndex, memorySize) + 1;
            end
        end

        [sortfitness, sortedIndices] = sort(fitness);
        bestSolution = population(sortedIndices(1), :);  % 更新全局最优解
        bestFitness = fitness(sortedIndices(1));  % 更新全局最优适应度值
        % bestSolution = population(i, :);  % 更新全局最优解
        % bestFitness = fitness(i);  % 更新全局最优适应度值
        iterCurve(1,iter) = bestFitness;  % 记录迭代曲线

        % disp(['Iteration: ', num2str(iter), ' Best Fitness: ', num2str(bestFitness)]);
    end
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