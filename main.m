clear all
clc

% 1. Configuración del problema (Speed Reducer)
try
    [nDim, LB, UB, ~, GloMin, Obj] = ProbInfo(1);
    if ~isa(Obj, 'function_handle')
        error('Obj no es un function_handle válido');
    end
    
    % Conversión explícita de límites
    lb = LB(:)';  % Asegura vector fila 1x7
    ub = UB(:)';  % Asegura vector fila 1x7
    dim = nDim(1);   % 7 para Speed Reducer
    
    % Factor de violación para las 11 restricciones
    Vio = ones(1, 11);  % Vector de penalizaciones
    
    % Función de costo adaptada
    fhd = @(x) CostFunction(x, Vio, Obj);
    
catch ME
    error('Error en configuración inicial: %s', ME.message);
end

% 2. Parámetros del algoritmo
SearchAgents_no = 80;
Max_Iteration = 1000;
% 3. Preparación de resultados
ul_res = zeros(1, Max_Iteration);  % Solo 1 fila para 1 problema
res = inf(30, Max_Iteration);    % 30 ejecuciones independientes
best_global = inf;               % Mejor fitness global

% 4. Ejecución principal
fprintf('Iniciando optimización para Speed Reducer (dim=%d)...\n', dim);
for i = 1:30  % Corregido a 30 iteraciones
    try
        [~, ~, His_Fit] = ship_aid(fhd, SearchAgents_no, Max_Iteration, dim, lb, ub);
        res(i, :) = His_Fit;
        current_best = His_Fit(end);
        fprintf('Ejecución %d/30 completada - Mejor fitness: %f\n', i, current_best);
        
        % Actualizar mejor fitness global
        if current_best < best_global
            best_global = current_best;
        end
    catch ME
        warning('Error en iteración %d: %s', i, ME.message);
        res(i, :) = NaN(1, Max_Iteration);
        disp(ME.stack(1));
    end
end

% 5. Procesamiento de resultados
ul_res(1, :) = mean(res, 1, 'omitnan');
disp('Resultados finales:');
disp(['Mejor fitness promedio: ' num2str(ul_res(end))]);
disp(['Mejor fitness global: ' num2str(best_global)]);