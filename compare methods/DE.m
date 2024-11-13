% function [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter, F, CR, constraintsFunc)
function [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter, F, CR)
    % 初始化种群
    x = initialization(pop, ub, lb, dim);
    Fitness = zeros(1, pop);  % 初始化适应度函数

    % 计算初始适应度并处理约束
    for i = 1:pop
        % Fitness(i) = handle_constraints(fobj(x(i, :)), constraintsFunc(x(i, :)));
        Fitness(i) = fobj(x(i, :)); % 计算每个个体的适应度值
        
    end

    [SortFitness, IndexSort] = sort(Fitness);
    Best_Pos = x(IndexSort(1), :);
    Best_Score = SortFitness(1);
    IterCurve = zeros(1, MaxIter);

    for t = 1:MaxIter
        for i = 1:pop
            % 选择不同于当前个体的三个随机个体
            r = randperm(pop, 3);
            while any(r == i)
                r = randperm(pop, 3);
            end
            a = x(r(1), :);
            b = x(r(2), :);
            c = x(r(3), :);

            % 差分变异操作
            v = a + F * (b - c);

            % 边界检查
            v = BoundrayCheck(v, ub, lb, dim);

            % 交叉操作
            u = x(i, :);
            for j = 1:dim
                if rand < CR
                    u(j) = v(j);
                end
            end

            % 计算适应度（带有约束处理）
            % uFitness = handle_constraints(fobj(u), constraintsFunc(u));
            uFitness = fobj(u);

            % 选择操作
            if uFitness < Fitness(i)
                x(i, :) = u;
                Fitness(i) = uFitness;
                if uFitness < Best_Score
                    Best_Pos = u;
                    Best_Score = uFitness;
                end
            end
        end

        % 记录当前最优解
        IterCurve(1, t) = Best_Score;
    end
end

% 边界检查函数：确保解在上下界内
function x = BoundrayCheck(x, ub, lb, dim)
    x = min(max(x, lb), ub);  % 使用min和max确保解在边界内
end

% % 约束处理函数：通过惩罚函数处理违反约束的解
% function penalized_fitness = handle_constraints(fitness, constraints)
%     penalty = 0;
%     for i = 1:length(constraints)
%         if constraints(i) > 0  % 如果违反约束
%             penalty = penalty + constraints(i);  % 增加违反约束的惩罚
%         end
%     end
%     penalized_fitness = fitness + 10000 * penalty;  % 添加惩罚项
% end





% function [Best_Pos, Best_Score, IterCurve] = DE(pop, dim, ub, lb, fobj, MaxIter, F, CR)
%     % 初始化种群
%     x = initialization(pop, ub, lb, dim);
%     Fitness = zeros(1, pop); % 初始化适应度函数
%     for i = 1:pop
%         Fitness(i) = fobj(x(i, :)); % 计算每个个体的适应度值
%     end
%     [SortFitness, IndexSort] = sort(Fitness);
%     Best_Pos = x(IndexSort(1), :);
%     Best_Score = SortFitness(1);
%     IterCurve = zeros(1, MaxIter);
% 
%     for t = 1:MaxIter
%         for i = 1:pop
%             % 选择不同于当前个体的三个随机个体
%             r = randperm(pop, 3);
%             while any(r == i)
%                 r = randperm(pop, 3);
%             end
%             a = x(r(1), :);
%             b = x(r(2), :);
%             c = x(r(3), :);
% 
%             % 差分变异操作
%             v = a + F * (b - c);
% 
%             % 边界检查
%             v = BoundrayCheck(v, ub, lb, dim);
% 
%             % 交叉操作
%             u = x(i, :);
%             for j = 1:dim
%                 if rand < CR
%                     u(j) = v(j);
%                 end
%             end
% 
%             % 选择操作
%             uFitness = fobj(u);
%             if uFitness < Fitness(i)
%                 x(i, :) = u;
%                 Fitness(i) = uFitness;
%                 if uFitness < Best_Score
%                     Best_Pos = u;
%                     Best_Score = uFitness;
%                 end
%             end
%         end
% 
%         IterCurve(1,t) = Best_Score;
%     end
% end
% 
% function x = BoundrayCheck(x, ub, lb, dim)
%     for i = 1:size(x, 1)
%         for j = 1:dim
%             if x(i, j) > ub(j)
%                 x(i, j) = ub(j);
%             end
%             if x(i, j) < lb(j)
%                 x(i, j) = lb(j);
%             end
%         end
%     end
% end