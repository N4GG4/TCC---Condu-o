function [mesh, dt, nSteps] = buildMesh(params)

Nx    = params.Nx;
Ny    = params.Ny;
Lx    = params.Lx;
Ly    = params.Ly;
alpha = params.alpha;

dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);

mesh.Nx = Nx;
mesh.Ny = Ny;
mesh.dx = dx;
mesh.dy = dy;
mesh.x  = linspace(0, Lx, Nx);
mesh.y  = linspace(0, Ly, Ny);
mesh.Lx = Lx;
mesh.Ly = Ly;

% "lambda" difusivo que aparece no CFL
lambda = 2*alpha*(1/dx^2 + 1/dy^2);

% -------------------------------------------------------------------------
% ESCOLHA DO PASSO DE TEMPO
% -------------------------------------------------------------------------
if isfield(params,'useImplicit') && params.useImplicit
    % ---------------------- CASO IMPLÍCITO ----------------------
    %
    % PRIORIDADE:
    % 1) Se params.dt existe  -> usa dt fixo
    % 2) Senão, se params.nSteps existe -> divide tEnd em nSteps
    % 3) Senão, se params.CFLimp existe -> usa CFL "implícito" (pode ser > 1)
    % 4) Senão, cai num padrão (ex: CFLimp = 5)
    
    if isfield(params, 'dt')
        dt     = params.dt;
        nSteps = ceil(params.tEnd / dt);
        
    elseif isfield(params, 'nSteps')
        nSteps = params.nSteps;
        dt     = params.tEnd / nSteps;
        
    else
        % Usa um CFL para o implícito (pode ser bem maior que 1)
        if isfield(params, 'CFLimp')
            CFLimp = params.CFLimp;
        else
            CFLimp = 5;   % padrão (você pode mudar)
        end
        
        dt     = CFLimp / lambda;
        nSteps = ceil(params.tEnd / dt);
    end
    
else
    % ---------------------- CASO EXPLÍCITO / HÍBRIDO ----------------------
    %
    % Aqui segue o CFL clássico. Para o explícito puro, o limite é CFL <= 1.
    % Para o híbrido, você pode decidir se respeita isso ou não.
    
    CFL = params.CFL;   % você define lá no script principal
    dt  = CFL / lambda;
    nSteps = ceil(params.tEnd / dt);
end

end

