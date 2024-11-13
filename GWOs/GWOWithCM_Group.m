%…………………………………………灰狼算法主体………………………………………
function [Best_Pos,Best_Score,IterCurve_CMG]=GWOWithCM_Group(pop,dim,ub,lb,fobj,MaxIter,crossProb,mutationProb,max_group_num,consFun)

Alpha_Pos=zeros(1,dim);%初始化Alpha狼群
Alpha_Score=inf;
Beta_Pos=zeros(1,dim);%初始化Beta狼群
Beta_Score=inf;
Delta_Pos=zeros(1,dim);%初始化化Delta狼群
Delta_Score=inf;
% crossProb=0.1;
% mutationProb=0.01;
% rng(1);
x=initialization(pop,ub,lb,dim);%初始化种群

Fitness=zeros(1,pop);%初始化适应度函数
for i=1:pop
    Fitness(i)=fobj(x(i,:)); %计算每个个体的适应度值
end
[SortFitness,IndexSort]=sort(Fitness);
Alpha_Pos=x(IndexSort(1),:);
Alpha_Score=SortFitness(1);
Beta_Pos=x(IndexSort(2),:);
Beta_Score=SortFitness(2);
Delta_Pos=x(IndexSort(3),:);
Delta_Score=SortFitness(3);
Group_Best_Pos=Alpha_Pos;
Group_Best_Score=Alpha_Score;



group_iter=MaxIter/4;
pop=64;


%第一阶段
a=2;
[x1,IterCurve1]=GWOGroup(1,a,x,pop,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

x2=x1(1:2:end,:);
x3=x1(2:2:end,:); 

%更新三只头狼
Alpha_Pos=x1(1,:);
Beta_Pos=x1(2,:);
Delta_Pos=x1(3,:);
% disp(fobj(Alpha_Pos));
% disp(fobj(Beta_Pos));
% disp(fobj(Delta_Pos));
% disp(Alpha_Pos);
% disp(Beta_Pos);
% disp(Delta_Pos);

%第二阶段
% a=2-250*(2/1000);
a=2-2*(exp(group_iter*1/MaxIter)-1)/(exp(1)-1);
[newx2,IterCurve21]=GWOGroup(2,a,x2,pop/2,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x4=newx2(1:2:end,:);
x5=newx2(2:2:end,:);

[newx3,IterCurve22]=GWOGroup(2,a,x3,pop/2,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x6=newx3(1:2:end,:);
x7=newx3(2:2:end,:);

IterCurve2=[IterCurve21,IterCurve22];
%排序，取前250
IterCurve2=sort(IterCurve2,'descend');
IterCurve2=IterCurve2(1:500);

%更新三只头狼
step2=[newx2(1:3,:);newx3(1:3,:)];
Fitness_step2=zeros(1,6);
for i=1:6
    Fitness_step2(i)=fobj(step2(i,:));
end
[SortFitness,IndexSort]=sort(Fitness_step2);
Alpha_Pos=step2(IndexSort(1),:);
Beta_Pos=step2(IndexSort(2),:);
Delta_Pos=step2(IndexSort(3),:);

%第三阶段
% a=2-500*(2/1000);
a=2-2*(exp(group_iter*2/MaxIter)-1)/(exp(1)-1);
[newx4,IterCurve31]=GWOGroup(3,a,x4,pop/4,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x8=newx4(1:2:end,:);
x9=newx4(2:2:end,:);

[newx5,IterCurve32]=GWOGroup(3,a,x5,pop/4,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x10=newx5(1:2:end,:);
x11=newx5(2:2:end,:);

[newx6,IterCurve33]=GWOGroup(3,a,x6,pop/4,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x12=newx6(1:2:end,:);
x13=newx6(2:2:end,:);

[newx7,IterCurve34]=GWOGroup(3,a,x7,pop/4,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);
x14=newx7(1:2:end,:);
x15=newx7(2:2:end,:);
IterCurve3=[IterCurve31,IterCurve32,IterCurve33,IterCurve34];
IterCurve3=sort(IterCurve3,'descend');
IterCurve3=IterCurve3(1:250);
%更新三只头狼
step3=[newx4(1:3,:);newx5(1:3,:);newx6(1:3,:);newx7(1:3,:)];
Fitness_step3=zeros(1,12);
for i=1:12
    Fitness_step3(i)=fobj(step3(i,:));
end
[SortFitness,IndexSort]=sort(Fitness_step3);
Alpha_Pos=step3(IndexSort(1),:);
Beta_Pos=step3(IndexSort(2),:);
Delta_Pos=step3(IndexSort(3),:);



%第四阶段
% a=2-750*(2/1000);
a=2-2*(exp(group_iter*3/MaxIter)-1)/(exp(1)-1);
[newx8,IterCurve41]=GWOGroup(4,a,x8,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx9,IterCurve42]=GWOGroup(4,a,x9,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx10,IterCurve43]=GWOGroup(4,a,x10,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx11,IterCurve44]=GWOGroup(4,a,x11,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx12,IterCurve45]=GWOGroup(4,a,x12,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx13,IterCurve46]=GWOGroup(4,a,x13,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx14,IterCurve47]=GWOGroup(4,a,x14,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

[newx15,IterCurve48]=GWOGroup(4,a,x15,pop/8,dim,ub,lb,fobj,group_iter,Alpha_Pos,Beta_Pos,Delta_Pos,crossProb,mutationProb,consFun);

IterCurve4=[IterCurve41,IterCurve42,IterCurve43,IterCurve44,IterCurve45,IterCurve46,IterCurve47,IterCurve48];
IterCurve4=sort(IterCurve4,'descend');
IterCurve4=IterCurve4(1:250);

res=[newx8(1:3,:);newx9(1:3,:);newx10(1:3,:);newx11(1:3,:);newx12(1:3,:);newx13(1:3,:);newx14(1:3,:);newx15(1:3,:)];
%disp(res);
Fitness_step4=zeros(1,24);
for i=1:24
    Fitness_step4(i)=fobj(res(i,:));
end
[SortFitness,IndexSort]=sort(Fitness_step4);
Best_Pos=res(IndexSort(1),:);
Best_Score=SortFitness(1);

IterCurve_CMG=[IterCurve1,IterCurve2,IterCurve3,IterCurve4];
IterCurve_CMG(1,:)=sort(IterCurve_CMG,'descend');
IterCurve_CMG=IterCurve_CMG(1,1:1000);
% IterCurve_CMG=IterCurve4(1,1:1000);

end