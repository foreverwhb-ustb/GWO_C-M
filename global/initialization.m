%种群初始化函数
function x=initialization(pop,ub,lb,dim)
% rng(1);
for i=1:pop
    for j=1:dim
        x(i,j)=(ub(j)-lb(j))*rand()+lb(j);
    end
end
end
