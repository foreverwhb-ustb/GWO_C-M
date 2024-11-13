function [Best_Pos,Best_Score,IterCurve]=GWOWithCM(pop,dim,ub,lb,fobj,MaxIter,crossProb,mutationProb,constraintsFunc)
% function [Best_Pos,Best_Score,IterCurve]=GWOWithCM(pop,dim,ub,lb,fobj,MaxIter,crossProb,mutationProb)
    Alpha_Pos = zeros(1, dim);  % 初始化Alpha狼群
    Alpha_Score = inf;
    Beta_Pos = zeros(1, dim);   % 初始化Beta狼群
    Beta_Score = inf;
    Delta_Pos = zeros(1, dim);  % 初始化Delta狼群
    Delta_Score = inf;
    
    x = initialization(pop, ub, lb, dim);  % 初始化种群
    Fitness = zeros(1, pop);  % 初始化适应度函数
    
    % 计算初始适应度（带约束处理）
    for i = 1:pop
        Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));  % 适应度值 + 约束处理
        % Fitness(i)=fobj(x(i,:));
    end
    
    % 排序适应度，选择Alpha, Beta, Delta
    [Fitness,IndexSort] = sort(Fitness);
    Alpha_Pos = x(IndexSort(1), :);
    Alpha_Score = Fitness(1);
    Beta_Pos = x(IndexSort(2), :);
    Beta_Score = Fitness(2);
    Delta_Pos = x(IndexSort(3), :);
    Delta_Score = Fitness(3);
    Group_Best_Pos = Alpha_Pos;
    Group_Best_Score = Alpha_Score;
    
    % 开始迭代
    for t = 1:MaxIter
        a = 2 - 2 * ((exp(t / MaxIter) - 1) / (exp(1) - 1));  % 计算动态参数a
        
        for i = 1:pop
            for j = 1:dim
                % 根据Alpha狼群更新位置X1
                r1 = rand;
                r2 = rand;
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_Alpha = abs(C1 * Alpha_Pos(j) - x(i, j));
                X1 = Alpha_Pos(j) - A1 * D_Alpha;
                
                % 根据Beta狼群更新位置X2
                r1 = rand;
                r2 = rand;
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_Beta = abs(C2 * Beta_Pos(j) - x(i, j));
                X2 = Beta_Pos(j) - A2 * D_Beta;
                
                % 根据Delta狼群更新位置X3
                r1 = rand;
                r2 = rand;
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_Delta = abs(C3 * Delta_Pos(j) - x(i, j));
                X3 = Delta_Pos(j) - A3 * D_Delta;
                
                % 更新位置
                x(i,j) = (X1 + X2 + X3) / 3;
            end

            % 交叉操作
            if rand < crossProb
                tmpx = x;
                [px, py] = size(x);
                
                % 选择交叉对象 (Alpha, Beta, Delta)
                cgreyIndex = round(rand * 3 + 1);
                if cgreyIndex == 1
                    cgrey = Alpha_Pos;
                elseif cgreyIndex == 2
                    cgrey = Beta_Pos;
                else
                    cgrey = Delta_Pos;
                end
                
                % 选择交叉点
                cpoint = randi([1, dim-1]);
                x(i, :) = [tmpx(i, 1:cpoint), cgrey(cpoint+1:end)];
            end
            
            % 变异操作
            if rand < mutationProb
                mpoint = randi([1, dim]);  % 随机选择一个变异点
                x(i, mpoint) = lb(mpoint) + rand * (ub(mpoint) - lb(mpoint));  % 在上下界范围内随机变异
            end
        end

        % 边界检查，确保所有解在上下界之间
        x = BoundrayCheck(x, ub, lb, dim);
        
        % 计算当前种群的适应度，并处理约束
        for i = 1:pop
            Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));  % 带约束的适应度处理
            % Fitness(i)=fobj(x(i,:));
        end
        
        % 更新Alpha, Beta, Delta狼群
        for i = 1:pop
            if Fitness(i) < Alpha_Score  % 替换 Alpha
                Delta_Pos = Beta_Pos; Delta_Score = Beta_Score;
                Beta_Pos = Alpha_Pos; Beta_Score = Alpha_Score;
                Alpha_Pos = x(i,:); Alpha_Score = Fitness(i);
            elseif Fitness(i) < Beta_Score  % 替换 Beta
                Delta_Pos = Beta_Pos; Delta_Score = Beta_Score;
                Beta_Pos = x(i,:); Beta_Score = Fitness(i);
            elseif Fitness(i) < Delta_Score  % 替换 Delta
                Delta_Pos = x(i,:); Delta_Score = Fitness(i);
            end
        end
        
        % 记录每次迭代中的最优解
        Group_Best_Pos = Alpha_Pos;
        Group_Best_Score = Alpha_Score;
        IterCurve(1,t) = Group_Best_Score;
    end
    
    % 最终返回最优解
    Best_Pos = Group_Best_Pos;
    Best_Score = Group_Best_Score;
end

% 边界检查函数：确保解在上下界内
function x = BoundrayCheck(x, ub, lb, dim)
    x = min(max(x, lb), ub);  % 使用min和max确保解在边界内
end