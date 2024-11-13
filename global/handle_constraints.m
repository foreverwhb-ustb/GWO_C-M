% 处理约束函数的惩罚方法
function penalized_fitness = handle_constraints(fitness, constraints)
    penalty = 0;  % 初始化惩罚值
    for i = 1:length(constraints)
        if constraints(i) > 0
            penalty = penalty + constraints(i);  % 如果约束被违反，增加惩罚
        end
    end
    penalized_fitness = fitness + 10000 * penalty;  % 惩罚系数可以调整
end