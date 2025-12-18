function T = applyBC(T, mesh, params, bc)
% Aplica BCs Dirichlet/Neumann nas quatro bordas

Nx = mesh.Nx;
Ny = mesh.Ny;
dx = mesh.dx;
dy = mesh.dy;
k  = params.k;

% LEFT (x = 0, coluna 1)
if bc.left.type == 'D'
    T(:,1) = bc.left.value;
elseif bc.left.type == 'N'
    qn = bc.left.value;        % W/m^2, positivo entrando
    T(:,1) = T(:,2) - qn*dx/k; % -k (T1 - T2)/dx = qn
end

% RIGHT (x = Lx, coluna Nx)
if bc.right.type == 'D'
    T(:,Nx) = bc.right.value;
elseif bc.right.type == 'N'
    qn = bc.right.value;
    T(:,Nx) = T(:,Nx-1) - qn*dx/k;
end

% BOTTOM (y = 0, linha 1)
if bc.bottom.type == 'D'
    T(1,:) = bc.bottom.value;
elseif bc.bottom.type == 'N'
    qn = bc.bottom.value;
    T(1,:) = T(2,:) - qn*dy/k;
end

% TOP (y = Ly, linha Ny)
if bc.top.type == 'D'
    T(Ny,:) = bc.top.value;
elseif bc.top.type == 'N'
    qn = bc.top.value;
    T(Ny,:) = T(Ny-1,:) - qn*dy/k;
end
end
