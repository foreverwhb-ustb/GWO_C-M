function [Best_Pos,Best_Score,IterCurve] = GWO(pop, dim, ub, lb, fobj, MaxIter, constraintsFunc)
% function [Best_Pos,Best_Score,IterCurve] = GWO(pop, dim, ub, lb, fobj, MaxIter)
    % 初始化 α, β, δ 的位置和得分
    Alpha_Pos = zeros(1, dim);
    Alpha_Score = inf;  % 初始化为正无穷
    
    Beta_Pos = zeros(1, dim);
    Beta_Score = inf;
    
    Delta_Pos = zeros(1, dim);
    Delta_Score = inf;
    
    % 初始化种群
    Positions = initialization(pop,ub,lb,dim);
    
    % 记录每次迭代的最优值
    IterCurve = zeros(1, MaxIter);
    
    % 主循环
    for iter = 1:MaxIter
        for i = 1:pop
            % 获取当前个体的位置
            current_Pos = Positions(i, :);
            
            % 计算当前个体的目标值和约束违反情况
            fitness = fobj(current_Pos);  % 目标函数值

            % =====
            constraints = constraintsFunc(current_Pos);  % 约束值            
            % 使用惩罚函数处理约束问题
            fitness = handle_constraints(fitness, constraints);
            
            % 更新 α, β, δ
            if fitness < Alpha_Score
                Delta_Pos = Beta_Pos;
                Delta_Score = Beta_Score;
                
                Beta_Pos = Alpha_Pos;
                Beta_Score = Alpha_Score;
                
                Alpha_Pos = current_Pos;
                Alpha_Score = fitness;
            elseif fitness < Beta_Score
                Delta_Pos = Beta_Pos;
                Delta_Score = Beta_Score;
                
                Beta_Pos = current_Pos;
                Beta_Score = fitness;
            elseif fitness < Delta_Score
                Delta_Pos = current_Pos;
                Delta_Score = fitness;
            end
        end
        
        % 更新种群的位置
        a = 2 - iter * (2 / MaxIter);  % 缩放因子 a，随迭代次数减少
        
        for i = 1:pop
            for j = 1:dim
                r1 = rand();
                r2 = rand();
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                
                D_alpha = abs(C1 * Alpha_Pos(j) - Positions(i, j));
                X1 = Alpha_Pos(j) - A1 * D_alpha;
                
                r1 = rand();
                r2 = rand();
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                
                D_beta = abs(C2 * Beta_Pos(j) - Positions(i, j));
                X2 = Beta_Pos(j) - A2 * D_beta;
                
                r1 = rand();
                r2 = rand();
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                
                D_delta = abs(C3 * Delta_Pos(j) - Positions(i, j));
                X3 = Delta_Pos(j) - A3 * D_delta;
                
                % 更新位置
                Positions(i, j) = (X1 + X2 + X3) / 3;
            end
            % 边界处理：确保新的位置在上下界内
            Positions(i, :) = min(max(Positions(i, :), lb), ub);
        end
        
        % 记录每次迭代的最优值
        IterCurve(iter) = Alpha_Score;
        
        % 打印当前迭代信息
        % fprintf('Iteration %d: Best Score = %.6f\n', iter, Alpha_Score);
    end
    
    % 返回最优位置和得分
    Best_Pos = Alpha_Pos;
    Best_Score = Alpha_Score;
end



