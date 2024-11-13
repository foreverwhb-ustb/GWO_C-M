function [newx, IterCurve_CMG] = GWOGroup(step, a, x, pop, dim, ub, lb, fobj, MaxIter, Alpha, Beta, Delta,crossProb,mutationProb, constraintsFunc)
% function [newx, IterCurve_CMG] = GWOGroup(step, a, x, pop, dim, ub, lb, fobj, MaxIter, Alpha, Beta, Delta,crossProb,mutationProb)
   % 初始化 α, β, δ 的位置和得分
    Alpha_Pos = Alpha;
    Beta_Pos = Beta;
    Delta_Pos = Delta;

    Alpha_Score=fobj(Alpha_Pos);
    Beta_Score=fobj(Beta_Pos);
    Delta_Score=fobj(Delta_Pos);
    
    % 初始化种群
    Positions = initialization(pop,ub,lb,dim);
    


    % 初始化
    initIter = (step - 1) * MaxIter;
    EndIter = step * MaxIter;
    IterCurve_CMG = zeros(1, EndIter);  % 初始化收敛曲线
    
    for t = initIter:EndIter-1  % 从当前迭代到结束迭代
        % 逐渐减少 a 值
        a = 2 - 2 * exp(t / (MaxIter * 4) - 1) / (exp(1) - 1);
        
        % 对每个个体更新位置
        for i = 1:pop
            for j = 1:dim
                % 根据 Alpha 更新 X1
                r1 = rand;
                r2 = rand;
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_Alpha = abs(C1 * Alpha_Pos(j) - x(i, j));
                X1 = Alpha_Pos(j) - A1 * D_Alpha;

                % 根据 Beta 更新 X2
                r1 = rand;
                r2 = rand;
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_Beta = abs(C2 * Beta_Pos(j) - x(i, j));
                X2 = Beta_Pos(j) - A2 * D_Beta;

                % 根据 Delta 更新 X3
                r1 = rand;
                r2 = rand;
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_Delta = abs(C3 * Delta_Pos(j) - x(i, j));
                X3 = Delta_Pos(j) - A3 * D_Delta;

                % 更新位置
                x(i,j) = (X1 + X2 + X3) / 3;
            end
            %交叉操作
            if rand<crossProb
                tmpx=x;
                [px,py] = size(x);
    
                %选择交叉对象
                cgreyIndex=round(rand*3+1);
                if cgreyIndex>=1 && cgreyIndex < 2
                    cgrey=Alpha_Pos;
                elseif cgreyIndex>=2 && cgreyIndex < 3
                    cgrey=Beta_Pos;
                else
                    cgrey=Delta_Pos;
                end 
    
                %选择交叉点
                cpoint=round(rand*dim);
                if cpoint==0
                    cpoint=1;
                elseif cpoint==pop
                    cpoint=pop-1;
                end
    
                x(i,:)=[tmpx(i,1:cpoint),cgrey(cpoint+1:py)];
                % x(IndexSort(i),:)=[tmpx(IndexSort(i),1:cpoint),cgrey(cpoint+1:py)];
    
            end
            %变异操作
            [px,py] = size(x);
            if(rand<mutationProb)
                mpoint=round(rand*py);
                if mpoint <= 0
                    mpoint = 1;
                end
                x(i,mpoint)=Alpha_Pos(mpoint);
                % x(IndexSort(i),mpoint)=Alpha_Pos(mpoint);
            end
        end
        
        % 边界检查，防止越界
        x = BoundrayCheck(x, ub, lb, dim);
        
        % 计算当前种群适应度，并处理约束
        Fitness = arrayfun(@(idx) handle_constraints(fobj(x(idx,:)), constraintsFunc(x(idx,:))), 1:pop);
        % 
        % 更新 Alpha、Beta、Delta
        for i = 1:pop
            Fitness(i) = handle_constraints(fobj(x(i,:)), constraintsFunc(x(i,:)));
            % Fitness(i)=fobj(x(i,:));
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
        
        % 记录每次迭代中的最优值
        IterCurve_CMG(1, t+1) = Alpha_Score;
    end

    % 计算每个个体的适应度值，并对解进行排序
    Fitness = arrayfun(@(idx) handle_constraints(fobj(x(idx,:)), constraintsFunc(x(idx,:))), 1:pop);
    % Fitness(i)=fobj(x(i,:));
    combinedX = [Fitness', x];  % 将适应度和个体连接
    sortedCombinedx = sortrows(combinedX, 1);  % 按适应度升序排列
    newx = sortedCombinedx(:, 2:end);  % 返回排序后的个体

    % for i=1:pop
    %     Fitness(i)=fobj(x(i,:)); %计算每个个体的适应度值
    % end
    % newFitness=transpose(Fitness);
    % 
    % combinedX=[newFitness, x];
    % sortedCombinedx=sortrows(combinedX,1);
    % newx=sortedCombinedx(:, 2:end);

end

% 边界检查函数
function x = BoundrayCheck(x, ub, lb, dim)
    x = min(max(x, lb), ub);  % 使用 min 和 max 限制在边界内
end