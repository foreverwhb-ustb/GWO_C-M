function [ub,lb,dim] = setRangeAndDim(fn)
dim=30;

%单峰测试函数 1-7
if fn==1
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
% if fn==1
%     dim=100;
%     ub=ones(1,dim)*100;
%     lb=ones(1,dim)*(-100);
elseif fn==2
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==3
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==4
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==5
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==6
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==7
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);

%多峰测试函数 8-13
elseif fn==8
    ub=ones(1,dim)*500;
    lb=ones(1,dim)*(-500);
elseif fn==9
    ub=ones(1,dim)*5.12;
    lb=ones(1,dim)*(-5.12);
elseif fn==10
    ub=ones(1,dim)*32;
    lb=ones(1,dim)*(-32);
elseif fn==11
    ub=ones(1,dim)*600;
    lb=ones(1,dim)*(-600);
elseif fn==12
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);
elseif fn==13
    ub=ones(1,dim)*50;
    lb=ones(1,dim)*(-50);

%复合测试函数 14-23
elseif fn==14
    dim=2;
    ub=ones(1,dim)*65;
    lb=ones(1,dim)*(-65);
elseif fn==15
    dim=4;
    ub=ones(1,dim)*5;
    lb=ones(1,dim)*(-5);
elseif fn==16
    dim=2;
    ub=ones(1,dim)*5;
    lb=ones(1,dim)*(-5);
elseif fn==17
    dim=2;
    ub=ones(1,dim)*5;
    lb=ones(1,dim)*(-5);
elseif fn==18
    dim=2;
    ub=ones(1,dim)*2;
    lb=ones(1,dim)*(-2);
elseif fn==19
    dim=3;
    ub=ones(1,dim)*3;
    lb=ones(1,dim)*1;
elseif fn==20
    dim=6;
    ub=ones(1,dim)*1;
    lb=ones(1,dim)*0;
elseif fn==21
    dim=4;
    ub=ones(1,dim)*10;
    lb=ones(1,dim)*0;
elseif fn==22
    dim=4;
    ub=ones(1,dim)*10;
    lb=ones(1,dim)*0;
elseif fn==23
    dim=4;
    ub=ones(1,dim)*10;
    lb=ones(1,dim)*0;

%CEC2019 24-33
elseif fn==24
    dim=9;
    ub=ones(1,dim)*8192;
    lb=ones(1,dim)*(-8192);
elseif fn==25
    dim=16;
    ub=ones(1,dim)*16384;
    lb=ones(1,dim)*(-16384);
elseif fn==26
    dim=18;
    ub=ones(1,dim)*4;
    lb=ones(1,dim)*(-4);
elseif fn==27
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==28
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==29
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==30
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==31
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==32
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==33
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
% CEC2022  F34-F45
elseif fn==34
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==35
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==36
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==37
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==38
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==39
    dim=2;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==40
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==41
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==42
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==43
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==44
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);
elseif fn==45
    dim=10;
    ub=ones(1,dim)*100;
    lb=ones(1,dim)*(-100);

% CEC2006 G9 G10
elseif fn==99
    dim=7;
    ub=ones(1,dim)*10;
    lb=ones(1,dim)*(-10);
elseif fn==100
    dim=8;
    lb = [100, 1000, 1000, 10, 10, 10, 10, 10];  % 下界
    ub = [10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000];  % 上界
end
end

