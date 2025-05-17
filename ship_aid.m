function [xbest,fbest,conv] = ship_aid(fhd,SearchAgents_no,Max_Iteration,dim,lb,ub,varargin)

if size(ub,2)==1
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end

ship_pos = initialization(SearchAgents_no,dim,ub,lb);
ship_v = initialization(SearchAgents_no,dim,ub,lb);
C = 2;
ship_fitness = inf*ones(SearchAgents_no,1);
ship_pos_new = ship_pos;
for iter = 1:Max_Iteration
    for i = 1:SearchAgents_no
        Flag4ub=ship_pos_new(i,:)>ub;
        Flag4lb=ship_pos_new(i,:)<lb;
        ship_pos_new(i,:)=(ship_pos_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    ship_fitness_new = feval(fhd,ship_pos_new',varargin{:});
    ship_fitness_new = ship_fitness_new';
    for i = 1:SearchAgents_no
        if ship_fitness_new(i) < ship_fitness(i)
            ship_pos(i,:) = ship_pos_new(i,:);
            ship_fitness(i) = ship_fitness_new(i);
        end
    end
    [fbest,I] = min(ship_fitness);
    xbest = ship_pos(I,:);
    F_wind = 3*rand(1,dim);
    m = ones(1,dim);
    a = F_wind./m;
    beta=3/2;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    for i = 1:SearchAgents_no
        ship_v(i,:) = ship_v(i,:) + a;
        p = rand;
        if p < 0.3
           % 静水情况 
           ship_pos_new(i,:) = ship_pos(i,:) + ship_v(i,:);
        elseif p >= 0.3 && p < 0.6
            u=randn(size(1,SearchAgents_no))*sigma;
            v=randn(size(1,SearchAgents_no));
            step=u./abs(v).^(1/beta);
            stepsize=0.01*step.*(ship_pos(i,:)-xbest);
            ship_pos_new(i,:) = ship_pos(i,:) + ship_v(i,:) + stepsize.*randn(size(1,SearchAgents_no));
        else
            u=randn(size(1,SearchAgents_no))*sigma;
            v=randn(size(1,SearchAgents_no));
            step=u./abs(v).^(1/beta);
            stepsize=0.01*step.*(ship_pos(i,:)-xbest);
            ship_pos_new(i,:) = ship_pos(i,:) + ship_v(i,:) - stepsize.*randn(size(1,SearchAgents_no));
        end
            
    end
    
    groups=4;
    groupsize=round(SearchAgents_no/4);
    
    % 创建分组
    for g = 1:groups
        group(g).fitness = inf * ones(groupsize,1);
        group(g).pos = ship_pos_new((g-1)*groupsize+1:(g-1)*groupsize+groupsize,:);
        group(g).fitness = ship_fitness((g-1)*groupsize+1:(g-1)*groupsize+groupsize); 
    end
    
    % 评选组内最有
    for g = 1:groups
        for i = 1:groupsize
            [~,I] = min(group(g).fitness);
            group(g).xbest = group(g).pos(I,:);
            group(g).fbest = group(g).fitness(I);
            [~,I] = max(group(g).fitness);
            group(g).xworse = group(g).pos(I,:);
            group(g).fworse = group(g).fitness(I);
        end
    end
    
    % 组内单独进化
    for g = 1:groups
        for i = 1:groupsize
            for tg = 1:groups
                dist_ib(tg) = dist(group(g).pos(i,:),group(tg).fbest);
                rssi(tg) = 1/dist_ib(tg);
            end
            [sorted_rssi,index] = sort(rssi);
            if 1/dist(group(g).pos(i,:),group(g).fbest) < (sorted_rssi(1)+sorted_rssi(2))/2
                VT = rand(1,dim).*(1/iter)^((group(g).fitness(i)-group(g).fbest)/(group(g).fbest-group(g).fworse));
                group(g).pos(i,:) = group(g).pos(i,:) + C.*rand.*(group(g).xbest-group(g).pos(i,:)) + VT;
            else
                VT = rand(1,dim).*(1/iter)^((group(g).fitness(i)-group(g).fbest)/(group(g).fbest-group(g).fworse));
                group(g).pos(i,:) = group(g).pos(i,:) + C.*rand.*(group(index(1)).xbest-group(g).pos(i,:)) + VT;
            end
            
        end
    end
    if rem(iter,20) == 0
        
        for g = 1:groups
            [~,index] = sort(group(g).fitness);
            for a = round(2*groupsize/3):groupsize
                group(g).pos(index(a),:) = group(mod(g+1,4)+1).xbest;
            end
        end
        
    end
    for g = 1:groups
        ship_pos_new((g-1)*groupsize+1:(g-1)*groupsize+groupsize,:) = group(g).pos;
    end
    conv(1,iter) = fbest;
end

end