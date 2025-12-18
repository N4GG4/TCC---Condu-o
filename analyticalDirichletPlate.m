function Tana = analyticalDirichletPlate(mesh, T_south, T_north, T_east, T_west, maxIter)
% Solução analítica estacionária de placa retangular 2D
% com temperaturas constantes nas 4 bordas (Dirichlet).
%
%   y = Ly
%     T_north
%  ^  +------------------+
%  |  |                  |
%  |  |                  |
%  |  |                  |
%  |  +------------------+
%  |  T_south
%  +--------------------------> x = Lx
%
%  x=0   -> T_west
%  x=Lx  -> T_east

Nx = mesh.Nx;
Ny = mesh.Ny;
x  = mesh.x;
y  = mesh.y;
Lx = mesh.Lx;
Ly = mesh.Ly;

if nargin < 7
    maxIter = 500;    % número de termos da série
end

Tana = zeros(Ny, Nx);

for i = 1:Ny
    for j = 1:Nx
        xx = x(j);
        yy = y(i);

        % --- contribuição da borda SUL (y = 0)
        Ts = 0;
        for k = 1:maxIter
            beta = k*pi/Lx;
            coef = (2*T_south/pi) * (((-1)^(k+1) + 1)/k);
            SUM  = coef * ( ...
                   -sinh(beta*yy)/tanh(beta*Ly) + cosh(beta*yy) ) ...
                   * sin(beta*xx);
            if ~isnan(SUM) && ~isinf(SUM)
                Ts = Ts + SUM;
            end
        end

        % --- contribuição da borda LESTE (x = Lx)
        Te = 0;
        for k = 1:maxIter
            beta = k*pi/Ly;
            coef = (2*T_east/pi) * (((-1)^(k+1) + 1)/k);
            SUM  = coef * sinh(beta*xx) .* sin(beta*yy) / sinh(beta*Lx);
            if ~isnan(SUM) && ~isinf(SUM)
                Te = Te + SUM;
            end
        end

        % --- contribuição da borda NORTE (y = Ly)
        Tn = 0;
        for k = 1:maxIter
            beta = k*pi/Lx;
            coef = (2*T_north/pi) * (((-1)^(k+1) + 1)/k);
            SUM  = coef * sinh(beta*yy) .* sin(beta*xx) / sinh(beta*Ly);
            if ~isnan(SUM) && ~isinf(SUM)
                Tn = Tn + SUM;
            end
        end

        % --- contribuição da borda OESTE (x = 0)
        Tw = 0;
        for k = 1:maxIter
            beta = k*pi/Ly;
            coef = (2*T_west/pi) * (((-1)^(k+1) + 1)/k);
            SUM  = coef * ( ...
                   -sinh(beta*xx)/tanh(beta*Lx) + cosh(beta*xx) ) ...
                   .* sin(beta*yy);
            if ~isnan(SUM) && ~isinf(SUM)
                Tw = Tw + SUM;
            end
        end

        Tana(i,j) = Ts + Te + Tn + Tw;
    end
end

% Impõe exatamente as bordas (só pra ficar bonitinho nos nós de fronteira)
Tana(1,:)   = T_south;
Tana(end,:) = T_north;
Tana(:,1)   = T_west;
Tana(:,end) = T_east;
end
