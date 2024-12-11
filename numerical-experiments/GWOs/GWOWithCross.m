function [Best_Pos,Best_Score,IterCurve] = GWOWithCross(pop, dim, ub, lb, fobj, MaxIter, crossProb, constraintsFunc)
% function [Best_Pos,Best_Score,IterCurve] = GWOWithCross(pop, dim, ub, lb, fobj, MaxIter, crossProb)
    Alpha_Pos = zeros(1, dim);  % 初始化Alpha狼群
    Alpha_Score = inf;
    Beta_Pos = zeros(1, dim);   % 初始化Beta狼群
    Beta_Score = inf;
    Delta_Pos = zeros(1, dim);  % 初始化Delta狼群
    Delta_Score = inf;

    % 初始化种群
    x = initialization(pop, ub, lb, dim);  
    Fitness = zeros(1, pop);  % 初始化适应度函数
    for i = 1:pop
        Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));  % 计算适应度值（带约束处理）
        % Fitness(i)=fobj(x(i,:)); %计算每个个体的适应度值
    end
    
    [SortFitness, IndexSort] = sort(Fitness);  % 对适应度排序
    Alpha_Pos = x(IndexSort(1), :);
    Alpha_Score = SortFitness(1);
    Beta_Pos = x(IndexSort(2), :);
    Beta_Score = SortFitness(2);
    Delta_Pos = x(IndexSort(3), :);
    Delta_Score = SortFitness(3);
    Group_Best_Pos = Alpha_Pos;
    Group_Best_Score = Alpha_Score;

    % 主循环：迭代优化过程
    for t = 1:MaxIter
        a = 2 - 2 * ((exp(t / MaxIter) - 1) / (exp(1) - 1));  % 非线性调整a值
        
        for i = 1:pop
            for j = 1:dim
                % Alpha狼群更新位置X1
                r1 = rand; r2 = rand;
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_Alpha = abs(C1 * Alpha_Pos(j) - x(i,j));
                X1 = Alpha_Pos(j) - A1 * D_Alpha;

                % Beta狼群更新位置X2
                r1 = rand; r2 = rand;
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_Beta = abs(C2 * Beta_Pos(j) - x(i,j));
                X2 = Beta_Pos(j) - A2 * D_Beta;

                % Delta狼群更新位置X3
                r1 = rand; r2 = rand;
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_Delta = abs(C3 * Delta_Pos(j) - x(i,j));
                X3 = Delta_Pos(j) - A3 * D_Delta;

                x(i,j) = (X1 + X2 + X3) / 3;  % 更新位置
            end

            % 交叉操作
            if rand < crossProb
                tmpx = x;
                [px, py] = size(x);
                cgrey = randi([1, pop]);  % 随机选择交叉对象
                while cgrey == i  % 确保交叉对象不是自身
                    cgrey = randi([1, pop]);
                end
                cpoint = randi([1, dim-1]);  % 选择交叉点
                
                % 进行交叉
                x(i,:) = [tmpx(i, 1:cpoint), tmpx(cgrey, cpoint+1:end)];
                x(cgrey,:) = [tmpx(cgrey, 1:cpoint), tmpx(i, cpoint+1:end)];
            end
        end

        x = BoundrayCheck(x, ub, lb, dim);  % 边界检查

        % 更新适应度并替换 Alpha、Beta、Delta狼
        for i = 1:pop
            Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));  % 带约束的适应度处理
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
        
        Group_Best_Pos = Alpha_Pos;
        Group_Best_Score = Alpha_Score;
        IterCurve(1,t) = Group_Best_Score;
    end

    Best_Pos = Group_Best_Pos;
    Best_Score = Group_Best_Score;
end


% 边界检查函数
function x = BoundrayCheck(x, ub, lb, dim)
    x = min(max(x, lb), ub);  % 将 x 限制在上下界内
end