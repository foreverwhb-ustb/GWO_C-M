function [gBestScore, cg_curve] = SaPSO(N, Max_iteration, lb, ub, dim, fobj)

    % PSO Parameters
    Vmax = 6;
    noP = N;
    wMax = 0.9;
    wMin = 0.2;
    c1 = 2;
    c2 = 2;
    alpha = 0.1;  % 自适应权重更新
    pMax = 0.9;   % 候选解生成策略最大概率
    pMin = 0.1;   % 候选解生成策略最小概率

    % Initialize strategy probabilities
    strategies = {'standard', 'mutation', 'crossover'}; % 可以扩展更多策略
    prob = [1/3, 1/3, 1/3];  % 初始每个策略的选择概率相等

    % Initializations
    iter = Max_iteration;
    vel = zeros(noP, dim);
    pBestScore = zeros(noP, 1);
    pBest = zeros(noP, dim);
    gBest = zeros(1, dim);
    cg_curve = zeros(1, iter);

    % Random initialization for agents
    pos = initialization(noP, ub, lb, dim);

    for i = 1:noP
        pBestScore(i) = inf;
    end

    % Initialize gBestScore for a minimization problem
    gBestScore = inf;

    % Main loop
    for loop = 1:iter
        for i = 1:noP
            % Ensure particles stay within bounds
            Flag4ub = pos(i, :) > ub;
            Flag4lb = pos(i, :) < lb;
            pos(i, :) = (pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % Evaluate fitness for each particle
            fitness = fobj(pos(i, :));

            if (pBestScore(i) > fitness)
                pBestScore(i) = fitness;
                pBest(i, :) = pos(i, :);
            end
            if (gBestScore > fitness)
                gBestScore = fitness;
                gBest = pos(i, :);
            end
        end

        % Update inertia weight
        w = wMax - loop * ((wMax - wMin) / iter);

        % Update the velocity and position of particles based on strategies
        for i = 1:noP
            strategy = roulette_wheel_selection(prob);  % 根据策略选择
            switch strategy
                case 'standard'
                    % 使用标准PSO更新规则
                    for j = 1:dim
                        vel(i, j) = w * vel(i, j) + c1 * rand() * (pBest(i, j) - pos(i, j)) + ...
                            c2 * rand() * (gBest(j) - pos(i, j));
                        vel(i, j) = min(max(vel(i, j), -Vmax), Vmax);
                        pos(i, j) = pos(i, j) + vel(i, j);
                    end
                case 'mutation'
                    % 添加变异策略
                    mutation_factor = 0.1;  % 定义变异强度
                    pos(i, :) = pos(i, :) + mutation_factor * randn(1, dim);
                case 'crossover'
                    % 添加交叉策略
                    if rand < 0.5
                        for j = 1:dim
                            pos(i, j) = pBest(i, j) + rand() * (gBest(j) - pBest(i, j));
                        end
                    end
            end
        end

        % 自适应更新策略选择的概率
        prob = update_strategy_probabilities(prob, gBestScore, alpha, pMax, pMin);

        % 记录每次迭代的全局最优值
        cg_curve(loop) = gBestScore;

        % 显示当前代数和全局最优值
        % fprintf('Iteration %d, Best Fitness: %.4f\n', loop, gBestScore);
    end
end

%% 帮助函数

% 初始化粒子的位置
function pos = initialization(noP, ub, lb, dim)
    pos = lb + (ub - lb) .* rand(noP, dim);
end

% 轮盘赌选择策略
function strategy = roulette_wheel_selection(prob)
    cumulative_prob = cumsum(prob);
    r = rand();
    strategy_idx = find(r <= cumulative_prob, 1);
    strategies = {'standard', 'mutation', 'crossover'};
    strategy = strategies{strategy_idx};
end

% 更新策略的概率
function new_prob = update_strategy_probabilities(prob, gBestScore, alpha, pMax, pMin)
    improvement = max(0, rand() * (1 - gBestScore));  % 通过全局最优的改进程度更新
    new_prob = prob + alpha * improvement;
    new_prob = max(min(new_prob, pMax), pMin);  % 限制在 [pMin, pMax] 范围内
    new_prob = new_prob / sum(new_prob);  % 确保总和为1
end
