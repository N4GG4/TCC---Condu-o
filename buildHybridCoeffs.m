function coeffs = buildHybridCoeffs(mesh, params, bc)
% Coeficientes para esquema híbrido 2D (bordas implícitas, interior explícito)
% Consistente com Eq. (3.77) do TCC.

Nx = mesh.Nx;
Ny = mesh.Ny;
dx = mesh.dx;
dy = mesh.dy;

k     = params.k;
rhoCp = params.rhoCp;

if isfield(params, 'qdot')
    qdot = params.qdot;
else
    qdot = 0.0;
end

% Arrays de coeficientes
aE  = zeros(Ny, Nx);
aW  = zeros(Ny, Nx);
aN  = zeros(Ny, Nx);
aS  = zeros(Ny, Nx);
SP  = zeros(Ny, Nx);
SU  = zeros(Ny, Nx);
theta = zeros(Ny, Nx);

% Coeficientes difusivos básicos (malha cartesiana uniforme)
aE(:) = k * dy / dx;
aW(:) = k * dy / dx;
aN(:) = k * dx / dy;
aS(:) = k * dx / dy;

% Fonte volumétrica linearizada Su = qdot*V, Sp = 0
SU(:) = qdot * dx * dy;
SP(:) = 0.0;

% --- Tratamento de Dirichlet nas bordas via Sp, Su (Patankar-like) ---

% West (i=1): T = T_west
Tw = bc.left.value;
for j = 1:Ny
    i = 1;
    aW(j,i) = 0.0;
    SP(j,i) = SP(j,i) - 2*k * dy / dx;
    SU(j,i) = SU(j,i) + 2*k * dy / dx * Tw;
end

% East (i=Nx): T = T_east
Te = bc.right.value;
for j = 1:Ny
    i = Nx;
    aE(j,i) = 0.0;
    SP(j,i) = SP(j,i) - 2*k * dy / dx;
    SU(j,i) = SU(j,i) + 2*k * dy / dx * Te;
end

% South (j=1): T = T_south
Ts = bc.bottom.value;
for i = 1:Nx
    j = 1;
    aS(j,i) = 0.0;
    SP(j,i) = SP(j,i) - 2*k * dx / dy;
    SU(j,i) = SU(j,i) + 2*k * dx / dy * Ts;
end

% North (j=Ny): T = T_north
Tn = bc.top.value;
for i = 1:Nx
    j = Ny;
    aN(j,i) = 0.0;
    SP(j,i) = SP(j,i) - 2*k * dx / dy;
    SU(j,i) = SU(j,i) + 2*k * dx / dy * Tn;
end

% Parâmetro temporal
dt  = mesh.dt;
aP0 = rhoCp * dx * dy / dt;

% Máscara do esquema híbrido: bordas implícitas, interior explícito
for j = 1:Ny
    for i = 1:Nx
        if (i == 1) || (i == Nx) || (j == 1) || (j == Ny)
            theta(j,i) = 1.0;  % célula de contorno -> implícita
        else
            theta(j,i) = 0.0;  % célula interna -> explícita
        end
    end
end

% aP de acordo com Eq. (3.78)
aP = theta .* (aE + aW + aN + aS - SP) + aP0;

% Guarda tudo em uma struct
coeffs.aE    = aE;
coeffs.aW    = aW;
coeffs.aN    = aN;
coeffs.aS    = aS;
coeffs.SP    = SP;
coeffs.SU    = SU;
coeffs.aP0   = aP0;
coeffs.aP    = aP;
coeffs.theta = theta;
coeffs.Nx    = Nx;
coeffs.Ny    = Ny;
coeffs.dx    = dx;
coeffs.dy    = dy;

end
