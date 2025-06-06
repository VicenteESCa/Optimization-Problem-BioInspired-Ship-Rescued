function z = CostFunction(x, VioFactor, Obj)
    if size(x,2) ~= 7
        if size(x,1) == 7
            x = x';
        else
            error('Dimensiones de x incorrectas.');
        end
    end

    [f, g, h] = Obj(x);
    
    v = zeros(size(x,1), 1);  % Asegurar que v tenga la dimensión correcta
    
    if ~isempty(g)
        ng = size(g,2);
        if length(VioFactor) < ng
            error('VioFactor debe tener al menos %d elementos', ng);
        end

        for i = 1:size(x,1)
            v(i) = sum(VioFactor(1:ng) .* max(0, g(i, :)).^2);  % Penalización cuadrática
        end
    end
    z = f + v;
end