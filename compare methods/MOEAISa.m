function [BestSol, BestFitness, IterCurve] = MOEAISa(N, Max_iter, dim, lb, ub, fobj)
    % 输入参数：
    % N - 种群大小
    % Max_iter - 最大迭代次数
    % dim - 维度
    % lb, ub - 解空间边界
    % fobj - 目标函数
    % dataset - 数据集

    
    % 示例 dataset，100 个样本，20 个特征，1 个分类标签
    dataset = rand(100, 20);  % 随机生成 100x20 的特征矩阵
    labels = randi([0, 1], 100, 1);  % 随机生成 0 或 1 的标签（例如二分类）
    dataset = [dataset, labels];  % 将特征和标签合并
    
    % 初始化参数
    Population = initialize_population(N, dim, lb, ub);  % 区间初始化
    Fitness = zeros(N, 1);
    W = calculate_feature_weights(dataset);  % 使用 ReliefF 计算特征权重
    
    % 评估初始种群
    for i = 1:N
        Fitness(i) = fobj(Population(i, :));
    end
    
    % 初始最优解
    [BestFitness, best_idx] = min(Fitness);
    BestSol = Population(best_idx, :);
    IterCurve = zeros(1, Max_iter);
    
    % 主循环
    for iter = 1:Max_iter
        % 选择父代
        Parents = select_parents(Population, Fitness, N);
        
        % 交叉操作（自适应交叉）
        Offspring = self_adaptive_crossover(Parents, W, dim);
        
        % 变异操作
        Offspring = mutation(Offspring, lb, ub);
        
        % 评估子代
        for i = 1:N
            NewFitness = fobj(Offspring(i, :));
            if NewFitness < Fitness(i)  % 更新
                Population(i, :) = Offspring(i, :);
                Fitness(i) = NewFitness;
                BestFitness = NewFitness;
                BestSol = Offspring(i, :);
                if NewFitness < BestFitness
                    BestFitness = NewFitness;
                    BestSol = Offspring(i, :);
                end
            end
        end
        
        % 记录每次迭代的最优解
        IterCurve(1,iter) = BestFitness;
    end
end

% 区间初始化函数
function Population = initialize_population(N, dim, lb, ub)
    Population = lb + (ub - lb) .* rand(N, dim);
end

% 交叉操作函数
function Offspring = self_adaptive_crossover(Parents, W, dim)
    % 使用父代和特征权重进行自适应交叉
    Offspring = zeros(size(Parents));
    for i = 1:size(Parents, 1)
        parent1 = Parents(i, :);
        parent2 = Parents(randi(size(Parents, 1)), :);
        similarity = sum(parent1 & parent2) / sum(parent1 | parent2);
        crossover_prob = similarity * 0.5 + 0.5;  % 自适应概率
        for d = 1:dim
            if rand() < crossover_prob
                Offspring(i, d) = parent1(d);
            else
                Offspring(i, d) = parent2(d);
            end
        end
    end
end

% 变异操作函数
function Offspring = mutation(Offspring, lb, ub)
    mutation_rate = 0.1;
    for i = 1:size(Offspring, 1)
        if rand() < mutation_rate
            Offspring(i, :) = lb + (ub - lb) .* rand(1, size(Offspring, 2));
        end
    end
end

% 特征权重计算函数（ReliefF）
function W = calculate_feature_weights(dataset)
    [~, W] = relieff(dataset(:,1:end-1), dataset(:,end), 20, 'method', 'classification');
end

% 选择父代函数
function Parents = select_parents(Population, Fitness, N)
    Parents = Population(tournament_selection(Fitness, N), :);
end

% 锦标赛选择函数
function idx = tournament_selection(Fitness, N)
    idx = zeros(N, 1);
    for i = 1:N
        competitors = randperm(N, 2);
        if Fitness(competitors(1)) < Fitness(competitors(2))
            idx(i) = competitors(1);
        else
            idx(i) = competitors(2);
        end
    end
end
