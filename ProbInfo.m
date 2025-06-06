function [nDim, LB, UB, Vio, GloMin, Obj] = ProbInfo(n)
    if n == 1  % Speed Reducer
        nDim = 7;  % Dimensiones del problema
        LB = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5.0];
        UB = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
        % Vector de factores de violación para las 11 restricciones
        Vio = ones(1, 11);  % Crear vector de unos de tamaño 11
        GloMin = 2994.0;  % Valor mínimo global conocido
        Obj = @f1;  % Función objetivo
    end
end

function [z, g, h] = f1(x)
    z = 0.7854*x(:,1).*x(:,2).^2.*(3.3333.*x(:,3).^2 + 14.9334.*x(:,3) - 43.0934) ...
        -1.508.*x(:,1).*(x(:,6).^2 + x(:,7).^2) + 7.477.*(x(:,6).^3 + x(:,7).^3) ...
        + 0.7854.*(x(:,4).*x(:,6).^2 + x(:,5).*x(:,7).^2);
    
    g = zeros(size(x,1), 11);
    g(:,1) = -x(:,1).*x(:,2).^2.*x(:,3) + 27;
    g(:,2) = -x(:,1).*x(:,2).^2.*x(:,3).^2 + 397.5;
    g(:,3) = -x(:,2).*x(:,6).^4.*x(:,3).*x(:,4).^(-3) + 1.93;
    g(:,4) = -x(:,2).*x(:,7).^4.*x(:,3)./x(:,5).^3 + 1.93;
    g(:,5) = 10.*x(:,6).^(-3).*sqrt(16.91.*10^6+(745.*x(:,4)./(x(:,2).*x(:,3))).^2)-1100;
    g(:,6) = 10.*x(:,7).^(-3).*sqrt(157.5.*10^6+(745.*x(:,5)./(x(:,2).*x(:,3))).^2)-850;
    g(:,7) = x(:,2).*x(:,3)-40;
    g(:,8) = -x(:,1)./x(:,2)+5;
    g(:,9) = x(:,1)./x(:,2)-12;
    g(:,10) = 1.5.*x(:,6)-x(:,4)+1.9;
    g(:,11) = 1.1.*x(:,7)-x(:,5)+1.9;
     % --- Penalización de restricciones (Método de Barrera) ---
    penalty = 1e6;  % Peso alto para penalizar violaciones
    z = z + penalty * sum(max(0, g).^2, 2);  % Suma cuadrática de violaciones

    h = [];
end