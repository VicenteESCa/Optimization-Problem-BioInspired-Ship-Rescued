function [xbest, fbest, conv] = ship_aid(fhd, SearchAgents_no, Max_Iteration, dim, lb, ub, varargin)
    % --- 1. Ajuste de límites ---
    if size(ub, 2) == 1 || isscalar(ub)
        ub = ub .* ones(1, dim);
        lb = lb .* ones(1, dim);
    end
    
    % --- 2. Parámetros del algoritmo ---
    k = 2e-5;               % Constante de ajuste
    T = Max_Iteration;      % Coeficiente de seguimiento
    K_com = 2;              % Constante de comunicación
    groups = 4;             % Número de grupos
    C = 1.0;               % Factor de aprendizaje
    
    % --- Inicialización de grupos y sus índices ---
    groupsize = floor(SearchAgents_no / groups);
    remainder = mod(SearchAgents_no, groups);
    groupIndices = cell(1, groups);
    
    % Asignar índices a cada grupo
    start_idx = 1;
    for g = 1:groups
        current_size = groupsize;
        if g == groups
            current_size = current_size + remainder;
        end
        end_idx = start_idx + current_size - 1;
        groupIndices{g} = start_idx:end_idx;
        start_idx = end_idx + 1;
    end
    
    % --- 3. Inicialización de variables ---
    F = 2 * randi([0, 1], SearchAgents_no, 1) - 1;  % Dirección aleatoria
    omega = zeros(SearchAgents_no, 1);  
    delta = zeros(SearchAgents_no, dim);  
    % C_coef = zeros(SearchAgents_no, dim);  
    angle = zeros(SearchAgents_no, dim);  
    ship_vel = zeros(SearchAgents_no, dim);  
    % ship_v = zeros(SearchAgents_no, dim);
    ship_fitness_new = zeros(SearchAgents_no, 1);
    
    % --- 4. Inicialización de posiciones y fitness ---
    ship_pos = initialization(SearchAgents_no, dim, ub, lb);  
    ship_fitness = zeros(SearchAgents_no, 1);
    for i = 1:SearchAgents_no
        ship_fitness(i) = fhd(ship_pos(i,:));
    end

    ship_pos_new = ship_pos;
    
    % --- 5. Inicializar MEJOR solución global ---
    [fbest, idx] = min(ship_fitness);
    xbest = ship_pos(idx, :);
    
    % --- 6. Control de convergencia ---
    conv = zeros(1, Max_Iteration);
    conv(1) = fbest;
    no_improvement = 0;
    max_no_improvement = 50;
    
    % --- Bucle principal ---
    for iter = 1:Max_Iteration
        try
            % ============================================
            % 1. MOVIMIENTO DE BARCOS
            % ============================================
            % Cálculo de ángulos
            for i = 1:SearchAgents_no
                norm_ship = norm(ship_pos(i,:));
                norm_xbest = norm(xbest);
                if norm_ship > 0 && norm_xbest > 0
                    cos_angle = dot(ship_pos(i,:), xbest) / (norm_ship * norm_xbest);
                    cos_angle = min(max(cos_angle, -1), 1);
                    angle(i,:) = acos(cos_angle);
                end
            end
            
            % Actualización de parámetros
            c = -2 + 4 * rand(SearchAgents_no, 1);
            delta = c .* F .* angle;
            %C_coef = (omega - k * delta(:,1)) ./ exp(-iter / T);
            omega = omega + k * delta(:,1);
            
            % Actualización de velocidad y posición
            ship_vel = ship_vel + omega .* (ub - lb) .* randn(SearchAgents_no, dim);
            ship_pos_new = ship_pos + C * (ship_vel + delta .* randn(SearchAgents_no, dim) .* (xbest - ship_pos));
            ship_pos_new = max(min(ship_pos_new, ub), lb);
            
            % ============================================
            % 2. EVALUACIÓN DE FITNESS
            % ============================================
            for i = 1:SearchAgents_no
                ship_fitness_new(i) = fhd(ship_pos_new(i,:));
            end
            
            improved = ship_fitness_new < ship_fitness;
            ship_pos(improved, :) = ship_pos_new(improved, :);
            ship_fitness(improved) = ship_fitness_new(improved);
            
            [current_best_fitness, idx] = min(ship_fitness);
            if current_best_fitness < fbest
                fbest = current_best_fitness;
                xbest = ship_pos(idx, :);
                no_improvement = 0;
            else
                no_improvement = no_improvement + 1;
            end
            
            % ============================================
            % 3. COMUNICACIÓN ENTRE GRUPOS
            % ============================================
            if rem(iter, 20) == 0
                for g = 1:groups
                    current_indices = groupIndices{g};
                    group(g).pos = ship_pos(current_indices, :);
                    group(g).fitness = ship_fitness(current_indices);
                    
                    [group(g).fbest, bestIndice] = min(group(g).fitness);
                    group(g).xbest = group(g).pos(bestIndice, :);
                    
                    [~, worstIdx] = sort(group(g).fitness, 'descend');
                    worst = worstIdx(1:ceil(length(current_indices)/3));
                    
                    Com1 = K_com * rand;
                    c2 = K_com * rand;
                    for w = 1:length(worst)
                        group(g).pos(worst(w), :) = group(g).pos(worst(w), :) + ...
                            Com1 * (group(g).xbest - group(g).pos(worst(w), :)) + ...
                            c2 * (xbest - group(g).pos(worst(w), :));
                        group(g).fitness(worst(w)) = fhd(group(g).pos(worst(w), :));
                    end
                    
                    ship_pos(current_indices, :) = group(g).pos;
                    ship_fitness(current_indices) = group(g).fitness;
                end
            end
            
            % Actualizar convergencia
            conv(iter) = fbest;
            
            % ============================================
            % 4. MANEJO DE CONVERGENCIA
            % ============================================
            if no_improvement > max_no_improvement
                fprintf('Convergencia alcanzada en iteración %d\n', iter);
                conv(iter+1:end) = fbest;
                break;
            end
            
        catch ME
            warning('Error en iteración %d: %s', iter, ME.message);
            disp(ME.stack(1));
            break;
        end
    end
end