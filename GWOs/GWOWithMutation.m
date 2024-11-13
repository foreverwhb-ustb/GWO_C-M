function [Best_Pos, Best_Score, IterCurve] = GWOWithMutation(pop, dim, ub, lb, fobj, MaxIter, mutationProb, constraintsFunc)
% function [Best_Pos, Best_Score, IterCurve] = GWOWithMutation(pop, dim, ub, lb, fobj, MaxIter, mutationProb)
    % 初始化 Alpha, Beta, Delta 三个狼群
    Alpha_Pos = zeros(1, dim);
    Alpha_Score = inf;
    Beta_Pos = zeros(1, dim);
    Beta_Score = inf;
    Delta_Pos = zeros(1, dim);
    Delta_Score = inf;

    % 初始化种群
    x = initialization(pop, ub, lb, dim);
    Fitness = zeros(1, pop);
    for i = 1:pop
        Fitness(i) = handle_constraints(fobj(x(i, :)), constraintsFunc(x(i, :)));  % 计算每个个体的适应度值，处理约束
        % Fitness(i)=fobj(x(i,:)); %计算每个个体的适应度值
    end

    % 排序适应度，得到初始 Alpha, Beta, Delta
    [SortFitness, IndexSort] = sort(Fitness);
    Alpha_Pos = x(IndexSort(1), :);
    Alpha_Score = SortFitness(1);
    Beta_Pos = x(IndexSort(2), :);
    Beta_Score = SortFitness(2);
    Delta_Pos = x(IndexSort(3), :);
    Delta_Score = SortFitness(3);
    
    Group_Best_Pos = Alpha_Pos;
    Group_Best_Score = Alpha_Score;

    % 迭代开始
    for t = 1:MaxIter
        a = 2 - 2 * ((exp(t / MaxIter) - 1) / (exp(1) - 1));  % 非线性调整a值

        for i = 1:pop
            for j = 1:dim
                % Alpha 狼更新位置 X1
                r1 = rand; r2 = rand;
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_Alpha = abs(C1 * Alpha_Pos(j) - x(i,j));
                X1 = Alpha_Pos(j) - A1 * D_Alpha;

                % Beta 狼更新位置 X2
                r1 = rand; r2 = rand;
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_Beta = abs(C2 * Beta_Pos(j) - x(i,j));
                X2 = Beta_Pos(j) - A2 * D_Beta;

                % Delta 狼更新位置 X3
                r1 = rand; r2 = rand;
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_Delta = abs(C3 * Delta_Pos(j) - x(i,j));
                X3 = Delta_Pos(j) - A3 * D_Delta;

                % 更新位置
                x(i,j) = (X1 + X2 + X3) / 3;
            end

            % 变异操作
            if rand < mutationProb
                mpoint = randi([1, dim]);  % 随机选择一个变异点
                x(i, mpoint) = lb(mpoint) + rand * (ub(mpoint) - lb(mpoint));  % 进行变异操作，随机选择一个范围内的新值
            end        
        end

        % 边界检查
        x = BoundrayCheck(x, ub, lb, dim);

        % 更新适应度并替换 Alpha, Beta, Delta
        for i = 1:pop
            Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));  % 计算适应度，带有约束处理
            % Fitness(i)=fobj(x(i,:)); %计算每个个体的适应度值
            if Fitness(i) < Alpha_Score
                Delta_Pos = Beta_Pos; Delta_Score = Beta_Score;
                Beta_Pos = Alpha_Pos; Beta_Score = Alpha_Score;
                Alpha_Pos = x(i,:); Alpha_Score = Fitness(i);
            elseif Fitness(i) < Beta_Score
                Delta_Pos = Beta_Pos; Delta_Score = Beta_Score;
                Beta_Pos = x(i,:); Beta_Score = Fitness(i);
            elseif Fitness(i) < Delta_Score
                Delta_Pos = x(i,:); Delta_Score = Fitness(i);
            end
        end

        % 记录迭代过程中的最优解
        Group_Best_Pos = Alpha_Pos;
        Group_Best_Score = Alpha_Score;
        IterCurve(1, t) = Group_Best_Score;
    end

    % 返回最终的最优解
    Best_Pos = Group_Best_Pos;
    Best_Score = Group_Best_Score;
end


% 边界检查函数：确保解在上下界内
function x = BoundrayCheck(x, ub, lb, dim)
    x = min(max(x, lb), ub);  % 将解限制在上下界之间
end