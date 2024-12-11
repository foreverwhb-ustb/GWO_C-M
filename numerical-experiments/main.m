%pop——种群数量
%dim——问题维度
%ub——变量上界，[1,dim]矩阵
%lb——变量下界，[1,dim]矩阵
%fobj——适应度函数（指针）
%MaxIter——最大迭代次数
%Best_Pos——x的最佳值
%Best_Score——最优适应度值
clear;
close all;

nexp=50;
pop=200;
crossProb=0.3;
mutationProb=0.6;
MaxIter=1000;
%fobj=@(x)fitness1(x);%设置适应度函数


%VAGWO
k_max=0.9; % Upper bound of the inertia weight
k_min=0.4; % Lower bound of the inertia weight
a_max=sqrt(2); % Upper bound of the acceleration coefficient
a_min=0; % Lower bound of the acceleration coefficient
c_max=1; % Upper bound of the leading wolves multipliers
c_min=0; % Lower bound of the leading wolves multipliers
elitism=1;% elitism=1: Elitism is applied; elitism=any other number: Elitism is not applied 
limvel=0.1;

%分组策略
group_num=8;

time1(nexp)=inf;
time2(nexp)=inf;
time3(nexp)=inf;
time4(nexp)=inf;
time5(nexp)=inf;
time6(nexp)=inf;
time7(nexp)=inf;
time8(nexp)=inf;
time9(nexp)=inf;
time10(nexp)=inf;
time11(nexp)=inf;
time12(nexp)=inf;
time13(nexp)=inf;
time14(nexp)=inf;
time15(nexp)=inf;
time16(nexp)=inf;
time17(nexp)=inf;
time18(nexp)=inf;


%输出结果到文件中
filename='output.txt';
fid=fopen(filename,'a');

for fn=1:1
    %确定适应度函数
    function_name=['F', num2str(fn)];
    fobj=str2func(function_name);

    constraintsFun=str2func('constraints');
    consFun=str2func('constraintsFunc_G09');

    [ub,lb,dim]=setRangeAndDim(fn);

    %VAGWO
    varmax=ub.*ones(1,dim); % Upper bound defined for the positions which can generally be a desired vector
    varmin=lb.*ones(1,dim);
    velmax=limvel*(varmax(1,1:dim)-varmin(1,1:dim)); % Upper bound defined for the velocities
    velmin=-velmax; % Lower bound defined for the velocities

    %nexp次独立实验
    for i=1:nexp
    
         fprintf("nexp：%f\n",i);

        tic;
        [Best_Pos,Best_Score,IterCurve]=GWO(pop,dim,ub,lb,fobj,MaxIter);
        time1(i)=toc;
        tic;
        [Best_Pos_C,Best_Score_C,IterCurve_C]=GWOWithCross(pop,dim,ub,lb,fobj,MaxIter,crossProb);
        time2(i)=toc;
        tic;
        [Best_Pos_M,Best_Score_M,IterCurve_M]=GWOWithMutation(pop,dim,ub,lb,fobj,MaxIter,mutationProb);
        time3(i)=toc;
        tic;
        [Best_Pos_CM,Best_Score_CM,IterCurve_CM]=GWOWithCM(pop,dim,ub,lb,fobj,MaxIter,crossProb,mutationProb);
        time4(i)=toc;
        tic;
        [Best_Pos_G,Best_Score_G,IterCurve_G]=GWOWithCM_Group_noCon(pop,dim,ub,lb,fobj,MaxIter,crossProb,mutationProb,group_num);
        time5(i)=toc;
        tic;
        [Alpha_pos,Alpha_score,IterCurve_so]=SOGWO(pop,dim,ub,lb,fobj,MaxIter);
        time6(i)=toc;
        tic;
        [Best_score_pso, IterCurve_pso] = PSO(pop,MaxIter,lb,ub,dim,fobj);
        time7(i)=toc;
        tic;
        [Best_Pos_de, Best_score_de, IterCurve_de] = DE(pop, dim, ub, lb, fobj, MaxIter, mutationProb, crossProb);
        time8(i)=toc;
        tic;
        [IterCurve_va,Best_Score_va,Best_Pos_va] = VAGWO(pop,dim,MaxIter,varmax,varmin,velmax,velmin,k_max,k_min,a_max,a_min,c_max,c_min,fobj,elitism);
        % [Best_Pos_va,Best_Score_va,IterCurve_va] = SOGWO(pop,dim,ub,lb,fobj,MaxIter);
        time9(i)=toc;
        tic;
        [Best_Score_rsm,IterCurve_rsm,Best_Pos_rsm]=RSMGWO(fobj,lb,ub,dim,MaxIter,pop);
        time10(i)=toc;
        tic;
        [Best_Score_hp,IterCurve_hp, Best_Pos_hp]=HGWOP(fobj,lb,ub,dim,MaxIter,pop);
        time11(i)=toc;
        tic;
        [Best_score_swo,Best_SW,IterCurve_swo]=SWO(pop,MaxIter,ub,lb,dim,fobj);
        time12(i)=toc;
        tic;
        [bestSolution, Best_score_sh, IterCurve_sh] = SHADE(fobj, dim, ub, lb, pop, MaxIter);
        time13(i)=toc;
        tic;
        [best_solution_lsh, Best_score_lsh, IterCurve_lsh] = LSHADE(fobj, constraintsFun, lb, ub, pop, MaxIter, dim);
        time14(i)=toc;
        tic;
        [bestSolution_jso, Best_score_jso, IterCurve_jso] = JSO(fobj, dim, ub, lb, pop, MaxIter);
        time15(i)=toc;
        tic;
        [bestSolution_dish, Best_score_dish, IterCurve_dish] = DISH(fobj, dim, ub, lb, pop, MaxIter);
        time16(i)=toc;
        tic;
        [bestSolution_mo, Best_score_mo,IterCurve_mo] = MOEAISa(pop, MaxIter, dim, lb, ub, fobj);
        time17(i)=toc;
        tic;
        [Best_score_sa, IterCurve_sa] = SaPSO(pop, MaxIter, lb, ub, dim, fobj);
        time18(i)=toc;

%       第i次实验中，每轮迭代的最优适应度值
        Iter(i,:)=IterCurve;
        Iter_C(i,:)=IterCurve_C;
        Iter_M(i,:)=IterCurve_M;
        Iter_CM(i,:)=IterCurve_CM;
        Iter_G(i,:)=IterCurve_G;
        Iter_so(i,:)=IterCurve_so;
        Iter_pso(i,:)=IterCurve_pso;
        Iter_de(i,:)=IterCurve_de;
        Iter_va(i,:)=IterCurve_va;
        Iter_rsm(i,:)=IterCurve_rsm;
        Iter_hp(i,:)=IterCurve_hp;
        Iter_swo(i,:)=IterCurve_swo;
        Iter_sh(i,:)=IterCurve_sh;
        Iter_lsh(i,:)=IterCurve_lsh;
        Iter_jso(i,:)=IterCurve_jso;
        Iter_dish(i,:)=IterCurve_dish;
        Iter_mo(i,:)=IterCurve_mo;
        Iter_sa(i,:)=IterCurve_sa;


        Group_Best_Score(i)=Best_Score;
        Group_Best_Score_C(i)=Best_Score_C;
        Group_Best_Score_M(i)=Best_Score_M;
        Group_Best_Score_CM(i)=Best_Score_CM;
        Group_Best_Score_G(i)=Best_Score_G;
        Group_Best_Score_so(i)=Alpha_score;
        Group_Best_Score_pso(i)=Best_score_pso;
        Group_Best_Score_de(i)=Best_score_de;
        Group_Best_Score_va(i)=Best_Score_va;
        Group_Best_Score_rsm(i)=Best_Score_rsm;
        Group_Best_Score_hp(i)=Best_Score_hp;
        Group_Best_Score_swo(i)=Best_score_swo;
        Group_Best_Score_sh(i)=Best_score_sh;
        Group_Best_Score_lsh(i)=Best_score_lsh;
        Group_Best_Score_jso(i)=Best_score_jso;
        Group_Best_Score_dish(i)=Best_score_dish;
        Group_Best_Score_mo(i)=Best_score_mo;
        Group_Best_Score_sa(i)=Best_score_sa;


    end 

    fprintf("time1：%f\n",sum(time1));
    fprintf("time2：%f\n",sum(time2));
    fprintf("time3：%f\n",sum(time3));
    fprintf("time4：%f\n",sum(time4));
    fprintf("time5：%f\n",sum(time5));
    fprintf("time6：%f\n",sum(time6));
    fprintf("time7：%f\n",sum(time7));
    fprintf("time8：%f\n",sum(time8));
    fprintf("time9：%f\n",sum(time9));
    fprintf("time10：%f\n",sum(time10));
    fprintf("time11：%f\n",sum(time11));
    fprintf("time12：%f\n",sum(time12));
    fprintf("time13：%f\n",sum(time13));
    fprintf("time14：%f\n",sum(time14));
    fprintf("time15：%f\n",sum(time15));
    fprintf("time16：%f\n",sum(time16));
    fprintf("time17：%f\n",sum(time17));
    fprintf("time18：%f\n",sum(time18));
  
    
    %计算最优值
    best_gwo(1,:)=sort(Group_Best_Score);
    best_gwo(2,:)=sort(Group_Best_Score_C);
    best_gwo(3,:)=sort(Group_Best_Score_M);
    best_gwo(4,:)=sort(Group_Best_Score_CM);
    best_gwo(5,:)=sort(Group_Best_Score_G);
    best_gwo(6,:)=sort(Group_Best_Score_pso); 
    best_gwo(7,:)=sort(Group_Best_Score_de);
    best_gwo(8,:)=sort(Group_Best_Score_sh);
    best_gwo(9,:)=sort(Group_Best_Score_swo);  
    best_gwo(10,:)=sort(Group_Best_Score_so);  
    best_gwo(11,:)=sort(Group_Best_Score_va); 
    best_gwo(12,:)=sort(Group_Best_Score_rsm);  
    best_gwo(13,:)=sort(Group_Best_Score_hp);
    best_gwo(14,:)=sort(Group_Best_Score_lsh);  
    best_gwo(15,:)=sort(Group_Best_Score_jso);  
    best_gwo(16,:)=sort(Group_Best_Score_dish);  
    best_gwo(17,:)=sort(Group_Best_Score_mo);  
    best_gwo(18,:)=sort(Group_Best_Score_sa);  



    %计算标准差
    std_GWO(fn)=std(Group_Best_Score);
    std_GWO_C(fn)=std(Group_Best_Score_C);
    std_GWO_M(fn)=std(Group_Best_Score_M);
    std_GWO_CM(fn)=std(Group_Best_Score_CM);
    std_GWO_G(fn)=std(Group_Best_Score_G);
    std_GWO_so(fn)=std(Group_Best_Score_so);
    std_GWO_pso(fn)=std(Group_Best_Score_pso);
    std_GWO_de(fn)=std(Group_Best_Score_de);
    std_GWO_va(fn)=std(Group_Best_Score_va);
    std_GWO_rsm(fn)=std(Group_Best_Score_rsm);
    std_GWO_hp(fn)=std(Group_Best_Score_hp);
    std_GWO_swo(fn)=std(Group_Best_Score_swo);
    std_GWO_sh(fn)=std(Group_Best_Score_sh);
    std_GWO_lsh(fn)=std(Group_Best_Score_lsh);
    std_GWO_jso(fn)=std(Group_Best_Score_jso);
    std_GWO_dish(fn)=std(Group_Best_Score_dish);
    std_GWO_mo(fn)=std(Group_Best_Score_mo);
    std_GWO_sa(fn)=std(Group_Best_Score_sa);


    %求最优适应度值的均值
    mean_F_best_score_GWO(fn)=mean(Group_Best_Score);
    mean_F_best_score_GWO_C(fn)=mean(Group_Best_Score_C);
    mean_F_best_score_GWO_M(fn)=mean(Group_Best_Score_M);
    mean_F_best_score_GWO_CM(fn)=mean(Group_Best_Score_CM);
    mean_F_best_score_GWO_G(fn)=mean(Group_Best_Score_G);
    mean_F_best_score_GWO_so(fn)=mean(Group_Best_Score_so);
    mean_F_best_score_GWO_pso(fn)=mean(Group_Best_Score_pso);
    mean_F_best_score_GWO_de(fn)=mean(Group_Best_Score_de);
    mean_F_best_score_GWO_va(fn)=mean(Group_Best_Score_va);
    mean_F_best_score_GWO_rsm(fn)=mean(Group_Best_Score_rsm);
    mean_F_best_score_GWO_hp(fn)=mean(Group_Best_Score_hp);
    mean_F_best_score_GWO_swo(fn)=mean(Group_Best_Score_swo);
    mean_F_best_score_GWO_sh(fn)=mean(Group_Best_Score_sh);
    mean_F_best_score_GWO_lsh(fn)=mean(Group_Best_Score_lsh);
    mean_F_best_score_GWO_jso(fn)=mean(Group_Best_Score_jso);
    mean_F_best_score_GWO_dish(fn)=mean(Group_Best_Score_dish);
    mean_F_best_score_GWO_mo(fn)=mean(Group_Best_Score_mo);
    mean_F_best_score_GWO_sa(fn)=mean(Group_Best_Score_sa);

    %求30次最优适应度值的最优值
    best_F_best_score_GWO(fn)=min(Group_Best_Score);
    best_F_best_score_GWO_C(fn)=min(Group_Best_Score_C);
    best_F_best_score_GWO_M(fn)=min(Group_Best_Score_M);
    best_F_best_score_GWO_CM(fn)=min(Group_Best_Score_CM);
    best_F_best_score_GWO_G(fn)=min(Group_Best_Score_G);
    best_F_best_score_GWO_so(fn)=min(Group_Best_Score_so);
    best_F_best_score_GWO_pso(fn)=min(Group_Best_Score_pso);
    best_F_best_score_GWO_de(fn)=min(Group_Best_Score_de);
    best_F_best_score_GWO_va(fn)=min(Group_Best_Score_va);
    best_F_best_score_GWO_rsm(fn)=min(Group_Best_Score_rsm);
    best_F_best_score_GWO_hp(fn)=min(Group_Best_Score_hp);
    best_F_best_score_GWO_swo(fn)=min(Group_Best_Score_swo);
    best_F_best_score_GWO_sh(fn)=min(Group_Best_Score_sh);
    best_F_best_score_GWO_lsh(fn)=min(Group_Best_Score_lsh);
    best_F_best_score_GWO_jso(fn)=min(Group_Best_Score_jso);
    best_F_best_score_GWO_dish(fn)=min(Group_Best_Score_dish);
    best_F_best_score_GWO_mo(fn)=min(Group_Best_Score_mo);
    best_F_best_score_GWO_sa(fn)=min(Group_Best_Score_sa);

    %求30次最优适应度值的最差值
    worst_F_best_score_GWO(fn)=max(Group_Best_Score);
    worst_F_best_score_GWO_C(fn)=max(Group_Best_Score_C);
    worst_F_best_score_GWO_M(fn)=max(Group_Best_Score_M);
    worst_F_best_score_GWO_CM(fn)=max(Group_Best_Score_CM);
    worst_F_best_score_GWO_G(fn)=max(Group_Best_Score_G);
    worst_F_best_score_GWO_so(fn)=max(Group_Best_Score_so);
    worst_F_best_score_GWO_pso(fn)=max(Group_Best_Score_pso);
    worst_F_best_score_GWO_de(fn)=max(Group_Best_Score_de);
    worst_F_best_score_GWO_va(fn)=max(Group_Best_Score_va);
    worst_F_best_score_GWO_rsm(fn)=max(Group_Best_Score_rsm);
    worst_F_best_score_GWO_hp(fn)=max(Group_Best_Score_hp);
    worst_F_best_score_GWO_swo(fn)=max(Group_Best_Score_swo);
    worst_F_best_score_GWO_sh(fn)=max(Group_Best_Score_sh);
    worst_F_best_score_GWO_lsh(fn)=max(Group_Best_Score_lsh);
    worst_F_best_score_GWO_jso(fn)=max(Group_Best_Score_jso);
    worst_F_best_score_GWO_dish(fn)=max(Group_Best_Score_dish);
    worst_F_best_score_GWO_mo(fn)=max(Group_Best_Score_mo);
    worst_F_best_score_GWO_sa(fn)=max(Group_Best_Score_sa);

    % % 
    fprintf(fid,'%s==========GWO===========GWO_C============GWO_M===========GWO_CM=========GWO_G===========PSO============DE=========SHADE=========SWO============SOGWO============VAGWO=========RSMGWO==========HGWOP==========LSHADE==========JSO==========DISH==========MOEASIa==========SaPSO====%s,%s,%s,%s,%s\n'...
        ,function_name,num2str(MaxIter),num2str(nexp),num2str(pop),num2str(crossProb),num2str(mutationProb));
    % fprintf(fid,'%s \t GWO \t GWO_C \t GWO_M \t GWO_CM \t GWO_G \t SOGWO \t PSO \t DE \t VAGWO \t RSMGWO \t HGWOP \t SWO===%s,%s,%s,%s,%s\n'...
    %     ,function_name,num2str(MaxIter),num2str(nexp),num2str(pop),num2str(crossProb),num2str(mutationProb));

    fprintf(fid,'time \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', num2str(mean(time1)) ...
        ,num2str(mean(time2)),num2str(mean(time3)) ...
        ,num2str(mean(time4)),num2str(mean(time5))...
        ,num2str(mean(time6)),num2str(mean(time7))...
        ,num2str(mean(time8)),num2str(mean(time9))...
        ,num2str(mean(time10)),num2str(mean(time11))...
        ,num2str(mean(time12)),num2str(mean(time13))...
        ,num2str(mean(time14)),num2str(mean(time15))...
        ,num2str(mean(time16)),num2str(mean(time17))...
        ,num2str(mean(time18)));
    fprintf(fid,'\n');

    fprintf(fid,'best \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', num2str(best_F_best_score_GWO(fn)) ...
        ,num2str(best_F_best_score_GWO_C(fn)),num2str(best_F_best_score_GWO_M(fn)) ...
        ,num2str(best_F_best_score_GWO_CM(fn)),num2str(best_F_best_score_GWO_G(fn))...
        ,num2str(best_F_best_score_GWO_pso(fn)),num2str(best_F_best_score_GWO_de(fn))...
        ,num2str(best_F_best_score_GWO_sh(fn)),num2str(best_F_best_score_GWO_swo(fn))...
        ,num2str(best_F_best_score_GWO_so(fn)),num2str(best_F_best_score_GWO_va(fn))...
        ,num2str(best_F_best_score_GWO_rsm(fn)),num2str(best_F_best_score_GWO_hp(fn)) ...
        ,num2str(best_F_best_score_GWO_lsh(fn)),num2str(best_F_best_score_GWO_jso(fn)) ...
        ,num2str(best_F_best_score_GWO_dish(fn)),num2str(best_F_best_score_GWO_mo(fn)) ...
        ,num2str(best_F_best_score_GWO_sa(fn)));
    fprintf(fid,'\n');
    fprintf(fid,'worst \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', num2str(worst_F_best_score_GWO(fn)) ...
        ,num2str(worst_F_best_score_GWO_C(fn)),num2str(worst_F_best_score_GWO_M(fn)) ...
        ,num2str(worst_F_best_score_GWO_CM(fn)),num2str(worst_F_best_score_GWO_G(fn))...
        ,num2str(worst_F_best_score_GWO_pso(fn)),num2str(worst_F_best_score_GWO_de(fn))...
        ,num2str(worst_F_best_score_GWO_sh(fn)),num2str(worst_F_best_score_GWO_swo(fn))...
        ,num2str(worst_F_best_score_GWO_so(fn)),num2str(worst_F_best_score_GWO_va(fn))...
        ,num2str(worst_F_best_score_GWO_rsm(fn)),num2str(worst_F_best_score_GWO_hp(fn))...
        ,num2str(worst_F_best_score_GWO_lsh(fn)),num2str(worst_F_best_score_GWO_jso(fn))...
        ,num2str(worst_F_best_score_GWO_dish(fn)),num2str(worst_F_best_score_GWO_mo(fn))...
        ,num2str(worst_F_best_score_GWO_sa(fn)));
    fprintf(fid,'\n');
    fprintf(fid,'mean \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', num2str(mean_F_best_score_GWO(fn)) ...
        ,num2str(mean_F_best_score_GWO_C(fn)),num2str(mean_F_best_score_GWO_M(fn)) ...
        ,num2str(mean_F_best_score_GWO_CM(fn)),num2str(mean_F_best_score_GWO_G(fn)) ...
        ,num2str(mean_F_best_score_GWO_pso(fn)),num2str(mean_F_best_score_GWO_de(fn))...
        ,num2str(mean_F_best_score_GWO_sh(fn)),num2str(mean_F_best_score_GWO_swo(fn))...
        ,num2str(mean_F_best_score_GWO_so(fn)),num2str(mean_F_best_score_GWO_va(fn))...
        ,num2str(mean_F_best_score_GWO_rsm(fn)),num2str(mean_F_best_score_GWO_hp(fn))...
        ,num2str(mean_F_best_score_GWO_lsh(fn)),num2str(mean_F_best_score_GWO_jso(fn))...
        ,num2str(mean_F_best_score_GWO_dish(fn)),num2str(mean_F_best_score_GWO_mo(fn)) ...
        ,num2str(mean_F_best_score_GWO_sa(fn)));
    fprintf(fid,'\n');
    fprintf(fid,'std \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', num2str(std_GWO(fn)) ...
        ,num2str(std_GWO_C(fn)),num2str(std_GWO_M(fn)) ...
        ,num2str(std_GWO_CM(fn)),num2str(std_GWO_G(fn))...
        ,num2str(std_GWO_pso(fn)),num2str(std_GWO_de(fn))...
        ,num2str(std_GWO_sh(fn)),num2str(std_GWO_swo(fn))...
        ,num2str(std_GWO_so(fn)),num2str(std_GWO_va(fn))...
        ,num2str(std_GWO_rsm(fn)),num2str(std_GWO_hp(fn)) ...
        ,num2str(std_GWO_lsh(fn)),num2str(std_GWO_jso(fn)) ...
        ,num2str(std_GWO_dish(fn)),num2str(std_GWO_mo(fn)) ...
        ,num2str(std_GWO_sa(fn)));
    fprintf(fid,'\n');




    % ==========补充
    sIter(1:MaxIter)=0;
    sIter_C(1:MaxIter)=0;
    sIter_M(1:MaxIter)=0;
    sIter_CM(1:MaxIter)=0;
    sIter_G(1:MaxIter)=0;
    sIter_so(1:MaxIter)=0;
    sIter_pso(1:MaxIter)=0;
    sIter_de(1:MaxIter)=0;
    sIter_va(1:MaxIter)=0;
    sIter_rsm(1:MaxIter)=0;
    sIter_hp(1:MaxIter)=0;
    sIter_swo(1:MaxIter)=0;
    sIter_sh(1:MaxIter)=0;
    sIter_lsh(1:MaxIter)=0;
    sIter_jso(1:MaxIter)=0;
    sIter_dish(1:MaxIter)=0;
    sIter_mo(1:MaxIter)=0;
    sIter_sa(1:MaxIter)=0;

    for num=1:nexp
        sIter=sIter+Iter(num,:);
        sIter_C=sIter_C+Iter_C(num,:);
        sIter_M=sIter_M+Iter_M(num,:);
        sIter_CM=sIter_CM+Iter_CM(num,:);
        sIter_G=sIter_G+Iter_G(num,:);
        sIter_so=sIter_so+Iter_so(num,:);
        sIter_pso=sIter_pso+Iter_pso(num,:);
        sIter_de=sIter_de+Iter_de(num,:);
        sIter_va=sIter_va+Iter_va(num,:);
        sIter_rsm=sIter_rsm+Iter_rsm(num,:);
        sIter_hp=sIter_hp+Iter_hp(num,:);
        sIter_swo=sIter_swo+Iter_swo(num,:);
        sIter_sh=sIter_sh+Iter_sh(num,:);
        sIter_lsh=sIter_lsh+Iter_lsh(num,:);
        sIter_jso=sIter_jso+Iter_jso(num,:);
        sIter_dish=sIter_dish+Iter_dish(num,:);
        sIter_mo=sIter_mo+Iter_mo(num,:);
        sIter_sa=sIter_sa+Iter_sa(num,:);
    end
    sIter=sIter/nexp;
    sIter_C=sIter_C/nexp;
    sIter_M=sIter_M/nexp;
    sIter_CM=sIter_CM/nexp;   
    sIter_G=sIter_G/nexp; 
    sIter_so=sIter_so/nexp; 
    sIter_pso=sIter_pso/nexp; 
    sIter_de=sIter_de/nexp; 
    sIter_va=sIter_va/nexp; 
    sIter_rsm=sIter_rsm/nexp; 
    sIter_hp=sIter_hp/nexp; 
    sIter_swo=sIter_swo/nexp; 
    sIter_sh=sIter_sh/nexp;
    sIter_lsh=sIter_lsh/nexp; 
    sIter_jso=sIter_jso/nexp; 
    sIter_dish=sIter_dish/nexp; 
    sIter_mo=sIter_mo/nexp; 
    sIter_sa=sIter_sa/nexp;
end


figure(1);
temp=50;
xc=0:50:999;
for i=1:20
    IterCurve1(i)=sIter(i*temp-temp+1);
    IterCurve_C1(i)=sIter_C(i*temp-temp+1);
    IterCurve_M1(i)=sIter_M(i*temp-temp+1);
    IterCurve_CM1(i)=sIter_CM(i*temp-temp+1);
    IterCurve_G1(i)=sIter_G(i*temp-temp+1);
    IterCurve_so1(i)=sIter_so(i*temp-temp+1);
    IterCurve_pso1(i)=sIter_pso(i*temp-temp+1);
    IterCurve_de1(i)=sIter_de(i*temp-temp+1);
    IterCurve_va1(i)=sIter_va(i*temp-temp+1);
    IterCurve_rsm1(i)=sIter_rsm(i*temp-temp+1);
    IterCurve_hp1(i)=sIter_hp(i*temp-temp+1);
    IterCurve_swo1(i)=sIter_swo(i*temp-temp+1);
    IterCurve_sh1(i)=sIter_sh(i*temp-temp+1);
    IterCurve_lsh1(i)=sIter_lsh(i*temp-temp+1);
    IterCurve_jso1(i)=sIter_jso(i*temp-temp+1);
    IterCurve_dish1(i)=sIter_dish(i*temp-temp+1);
    IterCurve_mo1(i)=sIter_mo(i*temp-temp+1);
    IterCurve_sa1(i)=sIter_sa(i*temp-temp+1);
end

semilogy(xc,IterCurve1,'b*-');
hold on;
% plot(xc,IterCurve1,'b*-');
plot(xc,IterCurve_C1,'r^-');
plot(xc,IterCurve_M1,'ko:','LineWidth',1);
plot(xc,IterCurve_CM1,'cs-','LineWidth',1);
plot(xc,IterCurve_G1,'md-','LineWidth',1);
plot(xc,IterCurve_pso1,'rx-','LineWidth',1);
plot(xc,IterCurve_de1,'yo:','LineWidth',1);
plot(xc,IterCurve_sh1,'b+-','LineWidth',1);
plot(xc,IterCurve_swo1,'c*-','LineWidth',1);
plot(xc,IterCurve_so1,'g^-','LineWidth',1);
plot(xc,IterCurve_va1,'k-.','LineWidth',1);
plot(xc,IterCurve_rsm1,'mp-','LineWidth',1);
plot(xc,IterCurve_hp1,'g+-','LineWidth',1);
plot(xc,IterCurve_lsh1,'cv-','LineWidth',1);
plot(xc,IterCurve_jso1,'r*-','LineWidth',1);
plot(xc,IterCurve_dish1,'m^-','LineWidth',1);
plot(xc,IterCurve_mo1,'bp-','LineWidth',1);
plot(xc,IterCurve_sa1,'kd-','LineWidth',1);

xlabel('Iteration');
% ylabel('Function value');
ylabel('Function value(log)');
legend('GWO','GWO\_C','GWO\_M','GWO\_CM','GWO\_CMG','PSO','DE'...
    ,'SHADE','SWO','SOGWO','VAGWO','RSMGWO','HGWOP','LSHADE' ...
    ,'JSO','DISH','MOEAISa','SaPSO');
% legend('GWO','GWO With Cross','GWO With Mutation','GWO With Cross-Mutation','GWO With Group');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function o = constraints(x)
    g1 = 0;
    g2 = 0;
    D = size(x,2);
    for i = 1:D
        g1 = g1 + (-x(i) * sin(2 * x(i)));
        g2 = g2 + x(i) * sin(x(i));
    end
    o = max(g1, 0) + max(g2, 0);
    o = o / 2;
end

%% 子函数用于定义表达式
% F1

function o = F1(x)
    o=sum(x.^2);
end

% F2

function o = F2(x)
    o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
    dim=size(x,2);
    o=0;
    for i=1:dim
        o=o+sum(x(1:i))^2;
    end
end

% F4

function o = F4(x)
    o=max(abs(x));
end

% F5

function o = F5(x)
    dim=size(x,2);
    o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
    o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
    dim=size(x,2);
    o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
    o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
    dim=size(x,2);
    o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
    dim=size(x,2);
    o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
    dim=size(x,2);
    o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
    dim=size(x,2);
    o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
    (1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
    dim=size(x,2);
    o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
    ((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
    aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
    -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
    
    for j=1:25
        bS(j)=sum((x'-aS(:,j)).^6);
    end
    o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
    aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
    bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
    o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
    o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
    o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
    o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
    aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
    pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
    o=0;
    for i=1:4
        o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F20

function o = F20(x)
    aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
    cH=[1 1.2 3 3.2];
    pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
    .2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
    o=0;
    for i=1:4
        o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F21

function o = F21(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:5
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F22

function o = F22(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:7
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F23

function o = F23(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    o=0;
    for i=1:10
        o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end
    
function o=Ufun(x,a,k,m)
    o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end

%CEC2019
function o = F24(x)
    % Chebyshev Function - Generalized Version (based on Storn's implementation)
    % Inputs:
    %   x - input vector of size D
    %   D - dimension of the problem
    % Output:
    %   f - function value (fitness)

    % Initialize variables
    D=9;
    a = 1.0; % Initial value for a
    b = 1.2; % Initial value for b
    sample = 32 * D; % Number of samples
    sumValue = 0.0; % Sum accumulator

    % Compute dx for boundary condition
    for j = 1:(D-2)
        dx = 2.4 * b - a;
        a = b;
        b = dx;
    end

    % Compute dy (step size in y direction)
    dy = 2.0 / sample;
    y = -1.0; % Start from -1

    % Main summation loop over samples
    for i = 0:sample
        px = x(1); % Initialize px with the first element of x
        for j = 2:D
            px = y * px + x(j); % Compute px recursively
        end
        if px < -1 || px > 1
            sumValue = sumValue + (1.0 - abs(px))^2;
        end
        y = y + dy; % Increment y
    end

    % Boundary condition for px with ±1.2 scaling
    for i = [-1, 1]
        px = x(1); % Initialize px with the first element of x
        for j = 2:D
            px = 1.2 * px + x(j); % Compute px recursively
        end
        if px < dx
            sumValue = sumValue + px^2;
        end
    end

    % Final fitness value
    o = sumValue+1;
end

function o=F25(x)
    D=size(x,2);
    n=sqrt(D);
    sum=0;
    H=zeros(n,n);
    for i=1:n
        for k=1:n
            H(i,k)=1/(i+k-1);
        end
    end

    Z=zeros(n,n);
    for i=1:n
        for k=1:n
            Z(i,k)=x(i+n*(k-1));
        end
    end

    I=eye(n);
    W=H.*Z-I;
    for i=1:n
        for k=1:n
            sum=sum+abs(W(i,k));
        end
    end
    o=sum+1;
end

function o=F26(x)
    D=size(x,2);
    n=D/3;
    sum=0;
    if n < 2  
        n = 2;
        D=6;
    end
    for i=1:n-1
        for j=i+1:n
            a=3*i-2;
            b=3*j-2;
            xd=x(a)-x(b);
            yd=x(a+1)-x(b+1);
            zd=x(a+2)-x(b+2);
            ed=xd^2+yd^2+zd^2;
            ud=ed^3;
            if ud>1.0e-10
                sum=sum+(1/ud-2)/ud;
            else
                sum=sum+1.0e20;
            end
        end
    end
    o=sum+12.7120622568;
end

function o=F27(x)
    D=size(x,2);
    sum=0;
    for i=1:D
        sum=sum+(x(i).^2-10*cos(2*pi*x(i))+10);
    end
    o=sum+1;
end

function o=F28(x)
    D=size(x,2);
    sum1=0;
    sum2=1;
    for i=1:D
        sum1=sum1+x(i)^2/4000;
        sum2=sum2*cos(x(i)/sqrt(i));
    end
    o=(sum1-sum2+1)+1;
end


function o = F29(x)
    % Weierstrass function
    a = 0.5;    % parameter a
    b = 3;      % parameter b
    k = 20;     % iteration depth
    D = size(x,2);  % dimension of the problem
    term_1 = 0;
    for i = 1:D
        term_2 = 0;
        for k = 0:k
            term_2 = term_2 + a^k * cos(2*pi*b^k*(x(i) + 0.5));
        end
        term_1 = term_1 + term_2;
    end
    term_3 = 0;
    for k = 0:k
        term_3 = term_3 + a^k * cos(pi*b^k);
    end
    o = (term_1 - D*term_3)+1;
end

function o=F30(x)
    D=size(x,2);
    tmp=0;
    for i=1:D
        z(i) = x(i)+420.9687462275036;
        if z(i)>500
            tmp=tmp+(500.0-mod(z(i),500))*sin(sqrt(abs(500.0-mod(z(i),500))));
            tmp2=-(z(i)-500)^2/(10000*D);
            tmp= tmp+tmp2;
        elseif z(i)<-500
            tmp=tmp+(-500.0+mod(abs(z(i)),500))*sin(sqrt(abs(mod(z(i),500)-500)));
            tmp2=-(z(i)+500)^2/10000*D;
            tmp= tmp+tmp2;
        else
            tmp=tmp+z(i)*sin(sqrt(abs(z(i))));
        end
    end
    o=(418.9829*D-tmp)+1;
end


function f = F31(x)
    % 调用sr_func进行平移和旋转操作
    nx=size(x,2);
    Os=ones(1,nx);
    Mr=ones(nx,nx);
    z = sr_func1(x, nx, Os, Mr, 1.0, 1, 1);

    f = 0;
    for i = 1 : nx - 1
        temp1 = sin(sqrt(z(i)^2 + z(i + 1)^2));
        temp1 = temp1^2;
        temp2 = 1 + 0.001*(z(i)^2 + z(i + 1)^2);
        f = f + 0.5+(temp1 - 0.5)/(temp2^2);
    end
    temp1 = sin(sqrt(z(nx)^2 + z(1)^2));
    temp1 = temp1^2;
    temp2 = 1 + 0.001*(z(nx)^2 + z(1)^2);
    f = f + 0.5+(temp1 - 0.5)/(temp2^2)+1;
end

function f = F32(x)
    % 输入参数:
    %   - x: 输入的自变量数组（维度与nx相关）
    %   - nx: 自变量的维度数量相关参数
    %   - Os: 用于平移操作相关的参数数组
    %   - Mr: 用于旋转操作相关的参数数组
    %   - s_flag: 标志位，用于控制是否进行平移操作的逻辑
    %   - r_flag: 标志位，用于控制是否进行旋转操作的逻辑
    % 输出参数:
    %   - f: 计算得到的HappyCat函数值
    nx=size(x,2);
    Os=ones(1,nx);
    Mr=ones(nx,nx);
    
    alpha = 1.0 / 8.0;
    
    % 调用sr_func函数进行平移和旋转等相关操作，得到处理后的结果
    sr_x = sr_func1(x, nx, Os, Mr, 5.0 / 100.0, 1, 1);
    
    r2 = 0;
    sum_z = 0;
    for i = 1 : nx % 在Matlab中循环索引通常从1开始
        sr_x(i) = sr_x(i) - 1.0; % 将坐标值向原点平移
        r2 = r2 + sr_x(i) ^ 2; % 计算坐标的平方和，Matlab中用^表示幂运算
        sum_z = sum_z + sr_x(i);
    end
    
    f = ((abs(r2 - nx) ^ (2 * alpha)) + (0.5 * r2 + sum_z) / nx + 0.5)/10 +1; % 根据公式计算HappyCat函数值
end

function f = F33(x)
    % ackley_func 函数用于计算Ackley函数的值
    % 输入参数:
    %   - x: 输入的自变量数组，维度与nx相关
    %   - nx: 自变量的维度数量
    %   - Os: 用于传递给sr_func函数的参数，和相关平移操作有关
    %   - Mr: 用于传递给sr_func函数的参数，和相关旋转操作有关
    %   - s_flag: 标志位，控制sr_func中平移相关操作逻辑
    %   - r_flag: 标志位，控制sr_func中旋转相关操作逻辑
    % 输出参数:
    %   - f: 计算得到的Ackley函数值
    
    % 定义常量
    nx=size(x,2);
    Os=ones(1,nx);
    Mr=ones(nx,nx);
    
    sum1 = 0;
    sum2 = 0;
    
    % 调用sr_func函数进行平移和旋转操作
    z = sr_func1(x, nx, Os, Mr, 1.0, 1, 1);
    
    for i = 1 : nx % 在Matlab中循环索引通常从1开始
        sum1 = sum1 + z(i)^2; % 计算平方和
        sum2 = sum2 + cos(2 * pi * z(i)); % 计算余弦项的和
    end
    
    sum1 = -0.2 * sqrt(sum1 / nx);
    sum2 = sum2 / nx;
    
    f = exp(1) - 20 * exp(sum1) - exp(sum2) + 20+1; % 根据Ackley函数公式计算最终值
end

function xshift = shiftfunc1(x, nx, Os)
    % shiftfunc 实现平移功能，将输入的x根据Os进行平移操作
    % 输入参数:
    %   - x: 原始输入数据数组
    %   - nx: 数据维度数量相关参数
    %   - Os: 平移的偏移量数组
    % 输出参数:
    %   - xshift: 平移后的数组
    
    xshift = zeros(1, nx); % 初始化平移后的结果数组
    for i = 1 : nx
        xshift(i) = x(i) - Os(i);
    end
end

function xrot = rotatefunc1(x, nx, Mr)
    % rotatefunc 实现旋转功能，根据旋转矩阵Mr对输入的x进行旋转操作
    % 输入参数:
    %   - x: 原始输入数据数组
    %   - nx: 数据维度数量相关参数
    %   - Mr: 旋转矩阵（以特定的展开形式传入，和原C语言代码对应）
    % 输出参数:
    %   - xrot: 旋转后的数组
    
    xrot = zeros(1, nx);
    for i = 1 : nx
        for j = 1 : nx
            xrot(i) = xrot(i) + x(j) * Mr((i - 1) * nx + j); % 根据旋转逻辑计算，注意Matlab索引从1开始调整计算方式
        end
    end
end

function sr_x = sr_func1(x, nx, Os, Mr, sh_rate, s_flag, r_flag)
    % sr_func 执行平移和旋转以及可能的缩放相关操作
    % 输入参数:
    %   - x: 原始输入数据
    %   - nx: 维度相关数量
    %   - Os: 平移相关参数数组
    %   - Mr: 旋转相关参数数组（以展开形式传入）
    %   - sh_rate: 缩放比例参数
    %   - s_flag: 标志位，控制是否进行平移操作
    %   - r_flag: 标志位，控制是否进行旋转操作
    % 输出参数:
    %   - sr_x: 经过一系列操作后的结果数组
    
    if s_flag == 1
        if r_flag == 1
            y = shiftfunc1(x, nx, Os); % 先进行平移操作
            for i = 1 : nx
                y(i) = y(i) * sh_rate; % 进行缩放操作
            end
            sr_x = rotatefunc1(y, nx, Mr); % 最后进行旋转操作
        else
            sr_x = shiftfunc1(x, nx, Os); % 只进行平移操作
            for i = 1 : nx
                sr_x(i) = sr_x(i) * sh_rate; % 进行缩放操作
            end
        end
    else
        if r_flag == 1
            y = x;
            for i = 1 : nx
                y(i) = y(i) * sh_rate; % 先进行缩放操作
            end
            sr_x = rotatefunc1(y, nx, Mr); % 再进行旋转操作
        else
            sr_x = x;
            for i = 1 : nx
                sr_x(i) = sr_x(i) * sh_rate; % 直接进行缩放操作
            end
        end
    end
end

% 任务调度
function o=F50(x)
    task1=[714,500,886,500,886,826, 886, 714, 500, 886, 671, 714, 886, 500, 500, 714, 714, 826, 714, 886, 886, 886, 500, 826, 886, 886, 886, 826, 886, 886, 826, 800, 714, 500, 714, 714, 886, 886, 886, 886, 886, 886, 500, 826, 886, 714, 886, 886, 886, 714, 886, 886, 886, 671, 714, 500, 500, 886, 886, 886];
    task2=[];
    vmTasks = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    dim = numel(task1);
    for i = 1:dim
    % 得到每个虚拟机分得的任务序号，存在 vmTasks Map中
        if ~isKey(vmTasks, x(i))
            taskList = [];
            taskList(i)=i;
            vmTasks(x(i)) = taskList;
        else
            vmTasks(x(i)) = [vmTasks(x(i)), i];
        end
    end
    
    makespan = 0;
    sumTime = 0;
    vmTaskKeys = keys(vmTasks);
    for k = 1:numel(vmTaskKeys)
        vmTask = vmTasks(vmTaskKeys{k});
        length = 0;
        for j = 1:numel(vmTask)
            taskId = vmTask(j);
            if taskId==0
                continue;
            end
            disp(taskId);
            length = length + task1(taskId);
        end
    
        % 计算每个虚拟机完成分配的任务所用的时间, 任务总长度/虚拟机的Mips
        runtime = length / 100;
        sumTime = sumTime + runtime;
        if makespan < runtime
            makespan = runtime;
        end
    end

    vmNum=numel(vmTaskKeys);
    
    RU=sumTime/makespan*vmNum;

    w1=0.6;
    w2=0.1;
    w3=0.3;
    o=w1*makespan+w2*sumTime+w3*RU;
end

% CEC2022

function sr_x = sr_func(x, Os, Mr,sh_rate, s_flag, r_flag)
    nx = size(x,2);
    sr_x = zeros(1, nx);

    if s_flag == 1
        if r_flag == 1
            y = shiftfunc(x, Os);
            y = y * sh_rate; % shrink to the original search range
            sr_x = rotatefunc(y, Mr);
        else
            sr_x = shiftfunc(x, Os);
            sr_x = sr_x * sh_rate; % shrink to the original search range
        end
    else
        if r_flag == 1
            y = x * sh_rate; % shrink to the original search range
            sr_x = rotatefunc(y, Mr);
        else
            sr_x = x * sh_rate; % shrink to the original search range
        end
    end
end

function xshift = shiftfunc(x, Os)
    nx = size(x,2);
    xshift = zeros(1, nx);
    for i = 1:nx
        xshift(i) = x(i) - Os(i);
    end
end

function xrot = rotatefunc(x, Mr)
    nx = size(x,2);
    xrot = zeros(1, nx);
    for i = 1:nx
        xrot(i) = 0;
        for j = 1:nx
            xrot(i) = xrot(i) + x(j) * Mr((i - 1) * nx + j);
        end
    end
end

% Shifted and full Rotated Zakharov Function
function f = F34(x)
    nx=size(x,2);
    z = zeros(1, nx);
    Os=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate

    f = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    for i = 1:nx
        xi = x(i);
        sum1 = sum1 + xi^2;
        sum2 = sum2 + 0.5 * (i + 1) * xi;
    end

    f = sum1 + sum2^2 + sum2^4 + 300;
end

% Shifted and full Rotated Rosenbrock's Function
function f = F35(x)
    nx=size(x,2);
    z = zeros(1, nx);
    Os=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;
    sr_func(x, Os, Mr, 2.048/100.0, s_flag, r_flag); % shift and rotate
    x(1) = x(1) + 1.0; % shift to origin

    f = 0.0;
    for i = 1:nx-1
        x(i+1) = x(i+1) + 1.0; % shift to origin
        tmp1 = x(i)^2 - x(i+1);
        tmp2 = x(i) - 1.0;
        f = f + 100.0 * tmp1^2 + tmp2^2 + 400;
    end
end

% Shifted and full Rotated Expanded Schaffer's f6 Function
function f = F36(x)
    nx=size(x,2);
    z = zeros(1, nx);
    Os=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;
    f = 0.0;
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate

    for i = 1:nx-1
        x(i) = sqrt(x(i)^2 + x(i+1)^2);
        tmp = sin(50.0 * x(i)^0.2);
        f = f + sqrt(x(i)) + sqrt(x(i)) * tmp^2;
    end

    f = (f^2) / (nx-1)^2 + 600;
end

% Shifted and full Rotated Non-Continuous Rastrigin's Function
function f = F37(x)
    nx=size(x,2);
    z = zeros(1, nx);
    Os=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;
    f = 0.0;
    sr_func(x, Os, Mr, 5.12/100.0, s_flag, r_flag); % shift and rotate

    for i = 1:nx
        f = f + (x(i)^2 - 10.0 * cos(2.0 * pi * x(i)) + 10.0) + 800;
    end
end

% Shifted and full Rotated Levy Function
function f = F38(x)
    nx = size(x,2);
    f = 0;
    z = zeros(1, nx);
    Os=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;
    sr_func(x, Os, Mr,1, s_flag, r_flag); % shift and rotate

    w = zeros(1, nx);
    sum1 = 0.0;
    for i = 1:nx
        w(i) = 1.0 + (x(i) - 1.0) / 4.0;
    end

    term1 = sin(pi * w(1))^2;
    term3 = (w(nx) - 1)^2 * (1 + sin(2 * pi * w(nx))^2);

    sum = 0.0;
    for i = 1:nx-1
        wi = w(i);
        newv = (wi - 1)^2 * (1 + 10 * sin(pi * wi + 1)^2);
        sum = sum + newv;
    end

    f = term1 + sum + term3 + 900;
end

% Hybrid Function 1(N=3)  No.6 in CEC2022
% Hybrid Function 1(N=3)  No.6 in CEC2022
function f = F39(x)

    nx = size(x,2);
    f = 0;
    z = zeros(1, nx);
    Os=zeros(1, nx);
    S=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=1;
    r_flag=1;

    cf_num = 3;
    fit = zeros(1, cf_num);
    Gp = [0.4,0.4,0.2];

    tmp = 0;
    for i = 1:cf_num-1
        G_nx(i) = ceil(Gp(i) * nx);
        tmp = tmp + G_nx(i);
    end
    G_nx(cf_num) = nx - tmp;

    G(1) = 0;
    for i = 2:cf_num
        G(i) = G(i-1) + G_nx(i-1);
    end

    z = zeros(1, nx);
    z=sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate

    for i = 1:nx
        y(i) = x(i);
    end

    i = 1;
    bent_cigar_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 2;
    hgbat_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 3;
    rastrigin_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);

    f(1) = sum(fit)+1800;
end

function f=F40(x)
    cf_num = 6;
    fit = zeros(1, cf_num);
    G = zeros(1, cf_num);
    Gp = [0.1, 0.2, 0.2, 0.2, 0.1, 0.2];
    nx = size(x,2);
    f = 0;
    z = zeros(1, nx*nx);
    Os=zeros(1, nx);
    S=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=0;
    r_flag=0;
    
    tmp = 0;
    for i = 1:cf_num-1
        G_nx(i) = ceil(Gp(i) * nx);
        tmp = tmp + G_nx(i);
    end
    G_nx(cf_num) = nx - tmp;
    
    G(1) = 1;
    for i = 2:cf_num
        G(i) = G(i-1) + G_nx(i-1);
    end
    
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    
    for i = 1:nx*nx
        y(i) = z(i);
    end
    
    i = 1;
    hgbat_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 2;
    katsuura_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 3;
    ackley_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 4;
    rastrigin_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 5;
    schwefel_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 6;
    schaffer_F7_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    
    f(1) = 0.0;
    for i = 1:cf_num
        f(1) = f(1) + fit(i);
    end
    f=f(1)+2000;
end

function f=F41(x)
    cf_num = 5;
    fit = zeros(1, cf_num);
    G = zeros(1, cf_num);
    G_nx = zeros(1, cf_num);
    Gp = [0.3, 0.2, 0.2, 0.1, 0.2];

    nx = size(x,2);
    f = 0;
    z = zeros(1, nx*nx);
    Os=zeros(1, nx);
    S=zeros(1, nx);
    Mr=zeros(1, nx*nx);
    s_flag=0;
    r_flag=0;
    
    tmp = 0;
    for i = 1:cf_num-1
        G_nx(i) = ceil(Gp(i) * nx);
        tmp = tmp + G_nx(i);
    end
    G_nx(cf_num) = nx - tmp;
    
    G(1) = 1;
    for i = 2:cf_num
        G(i) = G(i-1) + G_nx(i-1);
    end
    
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    
    for i = 1:nx*nx
        y(i) = z(i);
    end
    
    i = 1;
    katsuura_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 2;
    happycat_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 3;
    grie_rosen_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 4;
    schwefel_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    i = 5;
    ackley_func(y(G(i)+1:G(i)+G_nx(i)), fit(i), G_nx(i), Os, Mr, 0, 0);
    
    f(1) = 0.0;
    for i = 1:cf_num
        f(1) = f(1) + fit(i);
    end
    f=f(1)+2200;
end

function f=F42(x)
    cf_num = 5;
    fit = zeros(1, cf_num);
    delta = [10, 20, 30, 40, 50];
    bias = [0, 200, 300, 100, 400];

    nx = size(x,2);
    f = 0;
    Os=zeros(1, 5*nx*nx);
    Mr=zeros(1, 5*nx*nx);
    r_flag=0;

    i = 1;
    rosenbrock_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+4;
    
    i = 2;
    ellips_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+10;
    
    i = 3;
    bent_cigar_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+30;
    
    i = 4;
    discus_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+10;
    
    i = 5;
    ellips_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, 0);
    fit(i) = 10000 * fit(i) / 1e+10;

    f=cf_cal(x, nx, Os, delta, bias, fit, cf_num)+2100;

end

function f=F43(x)
    cf_num = 3;
    fit = zeros(1, cf_num);
    delta = [20, 10, 10];
    bias = [0, 200, 100];
    nx = size(x,2);
    f = 0;
    Os=zeros(1, 3*nx*nx);
    Mr=zeros(1, 3*nx*nx);
    r_flag=0;

    i = 1;
    schwefel_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, 0);
    
    i = 2;
    rastrigin_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    
    i = 3;
    hgbat_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    
    f=cf_cal(x, nx, Os, delta, bias, fit, cf_num)+2400;
end

function f=F44(x)
    cf_num = 5;
    fit = zeros(1, cf_num);
    delta = [20, 20, 30, 30, 20];
    bias = [0, 200, 300, 400, 200];

    nx = size(x,2);
    f = 0;
    Os=zeros(1, 5*nx*nx);
    Mr=zeros(1, 5*nx*nx);
    r_flag=0;
    
    i = 1;
    escaffer6_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 2e+7;
    
    i = 2;
    schwefel_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    
    i = 3;
    griewank_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 1000 * fit(i) / 100;
    
    i = 4;
    rosenbrock_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    
    i = 5;
    rastrigin_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+3;
    
    f=cf_cal(x, nx, Os, delta, bias, fit, cf_num)+2600;
end

function f=F45(x)
    cf_num = 6;
    fit = zeros(1, cf_num);
    delta = [10, 20, 30, 40, 50, 60];
    bias = [0, 300, 500, 100, 400, 200];

    nx = size(x,2);
    f = 0;
    Os=zeros(1, 6*nx*nx);
    Mr=zeros(1, 6*nx*nx);
    r_flag=0;
    
    i = 1;
    hgbat_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1000;
    
    i = 2;
    rastrigin_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+3;
    
    i = 3;
    schwefel_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 4e+3;
    
    i = 4;
    bent_cigar_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+30;
    
    i = 5;
    ellips_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 1e+10;
    
    i = 6;
    escaffer6_func(x, fit(i), nx, Os((i-1)*nx+1:i*nx), Mr((i-1)*nx*nx+1:i*nx*nx), 1, r_flag);
    fit(i) = 10000 * fit(i) / 2e+7;
    
    f=cf_cal(x, nx, Os, delta, bias, fit, cf_num)+2700;
end

function f=cf_cal(x, nx, Os, delta, bias, fit, cf_num)
    w = zeros(1, cf_num);
    w_max = 0;
    w_sum = 0;
    
    for i = 1:cf_num
        fit(i) = fit(i) + bias(i);
        w(i) = 0;
        for j = 1:nx
            w(i) = w(i) + (x(j) - Os((i-1)*nx + j))^2;
        end
        
        if w(i) ~= 0
            w(i) = 1/sqrt(w(i)) * exp(-w(i)/(2*nx*delta(i)^2));
        else
            w(i) = Inf;
        end
        
        if w(i) > w_max
            w_max = w(i);
        end
    end
    
    for i = 1:cf_num
        w_sum = w_sum + w(i);
    end
    
    if w_max == 0
        for i = 1:cf_num
            w(i) = 1;
        end
        w_sum = cf_num;
    end
    
    f(1) = 0.0;
    for i = 1:cf_num
        f(1) = f(1) + w(i)/w_sum * fit(i);
    end
end

function escaffer6_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z = zeros(1, nx);
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    
    f(1) = 0.0;
    for i = 1:nx-1
        temp1 = sin(sqrt(z(i)^2 + z(i+1)^2));
        temp1 = temp1^2;
        temp2 = 1.0 + 0.001 * (z(i)^2 + z(i+1)^2);
        f(1) = f(1) + 0.5 + (temp1 - 0.5) / (temp2^2);
    end
    
    temp1 = sin(sqrt(z(nx)^2 + z(1)^2));
    temp1 = temp1^2;
    temp2 = 1.0 + 0.001 * (z(nx)^2 + z(1)^2);
    f(1) = f(1) + 0.5 + (temp1 - 0.5) / (temp2^2);
end

function griewank_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z = zeros(1, nx);
    s = 0.0;
    p = 1.0;
    
    sr_func(x, Os, Mr, 600.0/100.0, s_flag, r_flag); % shift and rotate
    
    for i = 1:nx
        s = s + z(i)^2;
        p = p * cos(z(i) / sqrt(1.0 + i));
    end
    
    f(1) = 1.0 + s/4000.0 - p;
end


function rosenbrock_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z = zeros(1, nx);
    f(1) = 0.0;
    % sr_func(x, z, nx, Os, Mr, 2.048/100.0, s_flag, r_flag); % shift and rotate
    z(1) = z(1) + 1.0; % shift to origin
    for i = 1:nx-1
        z(i + 1) = z(i + 1) + 1.0; % shift to origin
        tmp1 = z(i) * z(i) - z(i + 1);
        tmp2 = z(i) - 1.0;
        f(1) = f(1) + 100.0 * tmp1^2 + tmp2^2;
    end
end

function ellips_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z = zeros(1, nx);
    f(1) = 0.0;
    % sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    for i = 1:nx
        f(1) = f(1) + 10^(6*i/(nx-1)) * z(i)^2;
    end
end

function discus_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z = zeros(1, nx);
    % sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    f(1) = 10^6 * z(1)^2;
    for i = 2:nx
        f(1) = f(1) + z(i)^2;
    end
end

function happycat_func(x, f, nx, Os, Mr, s_flag, r_flag)
    alpha = 1.0 / 8.0;
    
    sr_func(x, Os, Mr, 5.0 / 100.0, s_flag, r_flag); % shift and rotate
    z=x;

    r2 = 0.0;
    sum_z = 0.0;
    for i = 1:nx
        z(i) = z(i) - 1.0; % shift to origin
        r2 = r2 + z(i) * z(i);
        sum_z = sum_z + z(i);
    end
    f(1) = abs(r2 - nx)^(2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
end

function grie_rosen_func(x, f, nx, Os, Mr, s_flag, r_flag)
    f(1) = 0.0;

    sr_func(x, Os, Mr, 5.0 / 100.0, s_flag, r_flag); % shift and rotate
    z=x;

    z(1) = z(1) + 1.0; % shift to origin
    for i = 1:nx - 1
        z(i + 1) = z(i + 1) + 1.0; % shift to origin
        tmp1 = z(i) * z(i) - z(i + 1);
        tmp2 = z(i) - 1.0;
        temp = 100.0 * tmp1^2 + tmp2^2;
        f(1) = f(1) + (temp^2) / 4000.0 - cos(temp) + 1.0;
    end
    tmp1 = z(nx) * z(nx) - z(1);
    tmp2 = z(nx) - 1.0;
    temp = 100.0 * tmp1^2 + tmp2^2;
    f(1) = f(1) + (temp^2) / 4000.0 - cos(temp) + 1.0;
end

function katsuura_func(x, f, nx, Os, Mr, s_flag, r_flag)
    sr_func(x, Os, Mr, 5.0/100.0, s_flag, r_flag); % shift and rotate
    z=x;
    
    t1 = 0.0;
    t2 = 1.0;
    for i = 1:nx
        t1 = t1 + abs(z(i))^(i+1);
    end
    for i = 1:nx-1
        t2 = t2 * abs(z(i) - z(i+1)^(i+1));
    end
    f(1) = t1 + t2;
end

function ackley_func(x, f, nx, Os, Mr, s_flag, r_flag)
    sr_func(x, Os, Mr, 1.0, s_flag, r_flag); % shift and rotate
    z=x;
    sum1 = 0.0;
    sum2 = 0.0;
    for i = 1:nx
        sum1 = sum1 + z(i)^2;
        sum2 = sum2 + cos(2*pi*z(i));
    end
    f(1) = -20*exp(-0.2*sqrt(sum1/nx)) - exp(sum2/nx) + 20 + exp(1);
end



function schwefel_func(x, f, nx, Os, Mr, s_flag, r_flag)
    sr_func(x, Os, Mr, 1000.0/100.0, s_flag, r_flag); % shift and rotate
    z=x;
    sum = 0.0;
    for i = 1:nx
        sum = sum - z(i) * sin(sqrt(abs(z(i))));
    end
    f(1) = sum;
end

function schaffer_F7_func(x, f, nx, Os, Mr, s_flag, r_flag)
    sr_func(x, Os, Mr, 5.0/100.0, s_flag, r_flag); % shift and rotate
    z=x;
    sum1 = 0.0;
    sum2 = 0.0;
    for i = 1:nx-1
        sum1 = sum1 + (z(i)^2 + z(i+1)^2)^0.5;
        sum2 = sum2 + sin(sum1) * (sin((z(i)^2 + z(i+1)^2)^0.5)^2 - 0.5) / (1 + 0.001*(z(i)^2 + z(i+1)^2))^2;
    end
    f(1) = sum2;
end

function bent_cigar_func(x, f, nx, Os, Mr, s_flag, r_flag)
    f(1) = x(1)^2 + sum(10^6 * x(2:end).^2);
end

function hgbat_func(x, f, nx, Os, Mr, s_flag, r_flag)
    alpha = 1.0 / 4.0;
    z=x; % shift and rotate

    r2 = 0.0;
    sum_z = 0.0;
    for i = 1:nx
        z(i) = z(i) - 1.0;
        r2 = r2 + z(i)^2;
        sum_z = sum_z + z(i);
    end

    f(1) = (abs(r2^2 - sum_z^2))^(2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
end

function rastrigin_func(x, f, nx, Os, Mr, s_flag, r_flag)
    z=x; % shift and rotate

    f(1) = 0.0;
    for i = 1:nx
        f(1) = f(1) + (z(i)^2 - 10*cos(2*pi*z(i)) + 10);
    end
end

% CEC2006 G9 G10
% G9
% G09目标函数
function f = F99(x)
    f = (x(1) - 10)^2 + 5 * (x(2) - 12)^2 + x(3)^4 + 3 * (x(4) - 11)^2 + ...
        10 * x(5)^6 + 7 * x(6)^2 + x(7)^4 - 4 * x(6) * x(7) - 10 * x(6) - 8 * x(7);
end

% G09约束条件
function g = constraintsFunc_G09(x)
    g(1) = -127 + 2 * x(1)^2 + 3 * x(2)^4 + x(3) + 4 * x(4)^2 + 5 * x(5);
    g(2) = -282 + 7 * x(1) + 3 * x(2) + 10 * x(3)^2 + x(4) - x(5);
    g(3) = -196 + 23 * x(1) + x(2)^2 + 6 * x(6)^2 - 8 * x(7);
    g(4) = 4 * x(1)^2 + x(2)^2 - 3 * x(1) * x(2) + 2 * x(3)^2 + 5 * x(6) - 11 * x(7);
end


% G10
% G10目标函数
function f = F100(x)
    f = x(1) + x(2) + x(3);  % 目标函数
end

% G10约束条件
function g = constraintsFunc_G10(x)
    g(1) = -1 + 0.0025 * (x(4) + x(6));
    g(2) = -1 + 0.0025 * (x(5) + x(7) - x(4));
    g(3) = -1 + 0.01 * (x(8) - x(5));
    g(4) = -x(1) * x(6) + 833.33252 * x(4) + 100 * x(1) - 83333.333;
    g(5) = -x(2) * x(7) + 1250 * x(5) + x(2) * x(4) - 1250 * x(4);
    g(6) = -x(3) * x(8) + 1250000 + x(3) * x(5) - 2500 * x(5);
end
