% function [gBestScore, cg_curve] = PSO(N, Max_iteration, lb, ub, dim, fobj, constraintsFunc)
function [gBestScore, cg_curve] = PSO(N, Max_iteration, lb, ub, dim, fobj)

    % PSO Infotmation
    Vmax = 6;
    noP = N;
    wMax = 0.9;
    wMin = 0.2;
    c1 = 2;
    c2 = 2;

    % Initializations
    iter = Max_iteration;
    vel = zeros(noP, dim);
    pBestScore = zeros(noP, 1);
    pBest = zeros(noP, dim);
    gBest = zeros(1, dim);
    cg_curve = zeros(1, iter);

    % Random initialization for agents.
    pos = initialization(noP, ub, lb, dim);

    % Initialize pBestScore and gBestScore for minimization problem
    pBestScore(:) = inf;
    gBestScore = inf;

    for loop = 1:iter
        % 遍历每个粒子，计算适应度并处理约束
        for i = 1:noP
            % Return back the particles that go beyond the boundaries of the search space
            Flag4ub = pos(i,:) > ub;
            Flag4lb = pos(i,:) < lb;
            pos(i,:) = (pos(i,:) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % Calculate objective function with constraint handling
            % fitness = handle_constraints(fobj(pos(i,:)), constraintsFunc(pos(i,:)));
            fitness = fobj(pos(i,:));


            % 更新pBest和gBest
            if pBestScore(i) > fitness
                pBestScore(i) = fitness;
                pBest(i,:) = pos(i,:);
            end
            if gBestScore > fitness
                gBestScore = fitness;
                gBest = pos(i,:);
            end
        end

        % Update the inertia weight
        w = wMax - loop * ((wMax - wMin) / iter);

        % Update the Velocity and Position of particles
        for i = 1:noP
            for j = 1:dim
                vel(i,j) = w * vel(i,j) + c1 * rand() * (pBest(i,j) - pos(i,j)) + c2 * rand() * (gBest(j) - pos(i,j));
                
                % 限制速度
                if vel(i,j) > Vmax
                    vel(i,j) = Vmax;
                elseif vel(i,j) < -Vmax
                    vel(i,j) = -Vmax;
                end
                
                % 更新位置
                pos(i,j) = pos(i,j) + vel(i,j);
            end
        end

        % 记录当前全局最优解
        cg_curve(loop) = gBestScore;
    end
end

% % 约束处理函数：通过惩罚函数处理违反约束的解
% function penalized_fitness = handle_constraints(fitness, constraints)
%     penalty = 0;
%     for i = 1:length(constraints)
%         if constraints(i) > 0  % 违反约束
%             penalty = penalty + constraints(i);  % 增加违反约束的惩罚
%         end
%     end
%     penalized_fitness = fitness + 10000 * penalty;  % 添加惩罚系数
% end




% % Particle Swarm Optimization
% function [gBestScore, cg_curve] = PSO(N,Max_iteration,lb,ub,dim,fobj)
% 
% %PSO Infotmation
% 
% Vmax=6;
% noP=N;
% wMax=0.9;
% wMin=0.2;
% c1=2;
% c2=2;
% 
% % Initializations
% iter=Max_iteration;
% vel=zeros(noP,dim);
% pBestScore=zeros(noP);
% pBest=zeros(noP,dim);
% gBest=zeros(1,dim);
% cg_curve=zeros(1,iter);
% 
% % Random initialization for agents.
% pos=initialization(noP,ub,lb,dim);
% 
% for i=1:noP
%     pBestScore(i)=inf;
% end
% 
% % Initialize gBestScore for a minimization problem
% gBestScore=inf;
% 
% 
% for loop=1:iter
% 
%     % Return back the particles that go beyond the boundaries of the search
%     % space
%     Flag4ub=pos(i,:)>ub;
%     Flag4lb=pos(i,:)<lb;
%     pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% 
%     for i=1:size(pos,1)
%         %Calculate objective function for each particle
%         fitness=fobj(pos(i,:));
% 
%         if(pBestScore(i)>fitness)
%             pBestScore(i)=fitness;
%             pBest(i,:)=pos(i,:);
%         end
%         if(gBestScore>fitness)
%             gBestScore=fitness;
%             gBest=pos(i,:);
%         end
%     end
% 
%     %Update the W of PSO
%     w=wMax-loop*((wMax-wMin)/iter);
%     %Update the Velocity and Position of particles
%     for i=1:size(pos,1)
%         for j=1:size(pos,2)
%             vel(i,j)=w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+c2*rand()*(gBest(j)-pos(i,j));
% 
%             if(vel(i,j)>Vmax)
%                 vel(i,j)=Vmax;
%             end
%             if(vel(i,j)<-Vmax)
%                 vel(i,j)=-Vmax;
%             end
%             pos(i,j)=pos(i,j)+vel(i,j);
%         end
%     end
%     cg_curve(1,loop)=gBestScore;
% end
% 
% end