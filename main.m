clear all

% test_suite = 'cec13';
%CEC2013 28
%fhd=str2func('cec13_func'); Fnum = 28;
%CEC2014 30
% fhd=str2func('cec14_func'); Fnum = 30;
%CEC2017 29 去掉第二个
% fhd=str2func('cec17_func'); Fnum = 30;

Fnum=1

SearchAgents_no = 80;
D=[10,20,30,40,50,60,70,80,90];
dim = D(3);
lb = -100;
ub = 100;
Max_Iteration=1000;

ul_res = zeros(Fnum,Max_Iteration);
for f = 1:Fnum
    % new1
    res = zeros(30,Max_Iteration);
    parfor i = 1:30
        [~,~,His_Fit] = GWO(fhd,SearchAgents_no,Max_Iteration,dim,lb,ub,f);
%         display(['The best fitness is: ', num2str(His_Fit(Max_Iteration))]);
        res(i,:) = His_Fit;
    end
    ul_res(f,:) = sum(res,1)./30;
end
writematrix(ul_res,'./data_GWO.xlsx');

ul_res = zeros(Fnum,Max_Iteration);
for f = 1:Fnum
    % new1
    res = zeros(30,Max_Iteration);
    parfor i = 1:30
        [~,~,His_Fit] = ship_aid(fhd,SearchAgents_no,Max_Iteration,dim,lb,ub,f);
%         display(['The best fitness is: ', num2str(His_Fit(Max_Iteration))]);
        res(i,:) = His_Fit;
    end
    ul_res(f,:) = sum(res,1)./30;
end
writematrix(ul_res,'./data_new.xlsx');

ul_res = zeros(Fnum,Max_Iteration);
for f = 1:Fnum
    % new1
    res = zeros(30,Max_Iteration);
    parfor i = 1:30
        [~,~,His_Fit] = PSO(fhd,Atom_Num,Max_Iteration,dim,lb,ub,f);
%         display(['The best fitness is: ', num2str(His_Fit(Max_Iteration))]);
        res(i,:) = His_Fit;
    end
    ul_res(f,:) = sum(res,1)./30;
end

writematrix(ul_res,filename,'./data_pso.xlsx');

ul_res = zeros(Fnum,Max_Iteration);
for f = 1:Fnum
    % new1
    res = zeros(30,Max_Iteration);
    parfor i = 1:30
        [~,~,His_Fit] = FMO(fhd,Atom_Num,Max_Iteration,dim,lb,ub,f);
%         display(['The best fitness is: ', num2str(His_Fit(Max_Iteration))]);
        res(i,:) = His_Fit;
    end
    ul_res(f,:) = sum(res,1)./30;
end

writematrix(ul_res,filename,'./data_fmo.xlsx');


