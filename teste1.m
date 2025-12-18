clear; clc;

%% --- Parâmetros físicos
params.k     = 14.9;          % W/m.K
params.rhoCp = 1.0e7;         % J/m^3.K
params.alpha = params.k/params.rhoCp;
params.qdot  = 0.0;           % fonte volumétrica

%% --- Domínio
params.Lx = 0.02;             % m
params.Ly = 0.01;             % m

%% --- Malha
params.Nx  = 40;              % nós em x
params.Ny  = 40;              % nós em y

%% --- Tempo global
params.tEnd = 30;             % tempo final (s)

%% ----------------------------------------------------------
%% --- Escolha do esquema (muito importante)
%% ----------------------------------------------------------
scheme = "hybrid";   % "explicit", "implicit" ou "hybrid"

useImplicit = (scheme == "implicit");
useHybrid   = (scheme == "hybrid");
params.useImplicit = useImplicit;

switch scheme
    case "implicit"
        methodName = 'Implícito';
    case "explicit"
        methodName = 'Explícito FVM';
    case "hybrid"
        methodName = 'Híbrido';
end

%% ----------------------------------------------------------
%% --- Configuração do passo de tempo dependendo do esquema
%% ----------------------------------------------------------

if scheme == "explicit"
    % segue o CFL normal
    params.CFL = 1.0013;  % estável (<=1)

elseif scheme == "hybrid"
    % híbrido segue CFL também (explícito no interior)
    params.CFL = 1.0016;  

elseif scheme == "implicit"
    % Aqui você escolhe UMA das opções:
    
    % 1) controlar por nSteps (recomendado p/ ver diferenças)
    % params.nSteps = 200;

    % 2) controlar por dt fixo
    % params.dt = 0.1;

    % 3) controlar por CFL implícito (para aumentar dt)
    params.CFLimp = 1.0016;   % funciona muito bem para comparar
end

%% ----------------------------------------------------------
%% --- Monta malha e obtém dt e nSteps
%% ----------------------------------------------------------

[mesh, dt, nSteps] = buildMesh(params);
mesh.dt = dt;

fprintf('dx = %.4e, dy = %.4e, dt = %.4e, nSteps = %d\n', ...
    mesh.dx, mesh.dy, dt, nSteps);

%% --- Condição inicial
T0 = 25.0;                              
T  = T0*ones(mesh.Ny, mesh.Nx);

%% --- Condições de contorno
T_south = 50;    
T_north = 75;
T_west  = 100;
T_east  = 25;

bc.left.type   = 'D';  bc.left.value   = T_west;
bc.right.type  = 'D';  bc.right.value  = T_east;
bc.bottom.type = 'D';  bc.bottom.value = T_south;
bc.top.type    = 'D';  bc.top.value    = T_north;

T = applyBC(T, mesh, params, bc);

%% --- Analítico estacionário
maxIter = 300;
T_ana = analyticalDirichletPlate(mesh, T_south, T_north, T_east, T_west, maxIter);
normAnaL2 = norm(T_ana(:), 2);

%% --- Monta sistemas
switch scheme
    case "implicit"
        sysImplicitFVM = buildImplicitSystemFVM(mesh, params, bc);

    case "hybrid"
        sysHybrid = buildHybridSystemMatrix(mesh, params, bc);

    case "explicit"
        coeffsExp = buildExplicitFvmCoeffs(mesh, params, bc);
end

%% --- Probes
ixMid = round(mesh.Nx/2);
iyMid = round(mesh.Ny/2);

iyTop    = mesh.Ny - 1;
iyBottom = 2;

time     = (0:nSteps)' * dt;
probeCenter  = zeros(nSteps+1, 1);
probeTop     = zeros(nSteps+1, 1);
probeBottom  = zeros(nSteps+1, 1);
Tmean        = zeros(nSteps+1, 1);
errRelL2     = zeros(nSteps+1, 1);

probeCenter(1)  = T(iyMid, ixMid);
probeTop(1)     = T(iyTop, ixMid);
probeBottom(1)  = T(iyBottom, ixMid);
Tmean(1)        = mean(T(:));
errRelL2(1)     = norm(T(:) - T_ana(:), 2) / normAnaL2;

%% --- Loop no tempo
tic;
for n = 1:nSteps

    switch scheme
        case "implicit"
            T = solveHeatStepImplicitFVM(T, sysImplicitFVM);

        case "hybrid"
            T = solveHeatStepHybridMatrix(T, sysHybrid);

        case "explicit"
            T = solveHeatStepExplicitFVM(T, coeffsExp);
    end

    probeCenter(n+1)  = T(iyMid, ixMid);
    probeTop(n+1)     = T(iyTop, ixMid);
    probeBottom(n+1)  = T(iyBottom, ixMid);
    Tmean(n+1)        = mean(T(:));
    errRelL2(n+1)     = norm(T(:) - T_ana(:), 2) / normAnaL2;
end
tempoExecucao = toc;

fprintf('\nEsquema = %s | Tempo execução = %.4f s | Erro final = %.3e\n',...
    scheme, tempoExecucao, errRelL2(end));

%% --- Erro final da simulação ---
errFinal = errRelL2(end);
fprintf('Esquema: %s | Erro relativo L2 final = %.3e\n', methodName, errFinal);

%% --- Campo de temperatura numérico (pixel a pixel)
figure;
plotField(mesh, T);
title(['Campo de temperatura numérico (FVM - ' methodName ')']);

%% --- Mapa de cores (Temperatura) + Isocontornos Numérico vs Analítico
[X, Y] = meshgrid(mesh.x, mesh.y);

% níveis comuns (iguais para os dois)
nLevels = 12;
levels  = linspace(min(T_ana(:)), max(T_ana(:)), nLevels);

% escala de cor comum (pegando ambos pra ser justo)
Tmin = min([T(:); T_ana(:)]);
Tmax = max([T(:); T_ana(:)]);

figure;

% --- FUNDO: campo de temperatura (volume a volume)
pcolor(X, Y, T);
shading flat;              % ESSENCIAL → sem interpolação
axis equal tight;
xlabel('x (m)');
ylabel('y (m)');

% >>> TÍTULO ADICIONADO AQUI <<<
title(['Comparação bidimensional: Numérico (FVM) × Analítico — ' methodName]);

cb = colorbar;
cb.Label.String = 'T (°C)';
clim([Tmin Tmax]);          % mesma escala p/ numérico e analítico

hold on;

% --- contornos Numérico (preto, contínuo)
[Cnum, hNum] = contour(X, Y, T, levels, 'k', 'LineWidth', 1.6);

% --- contornos Analítico (vermelho, tracejado)
[Cana, hAna] = contour(X, Y, T_ana, levels, 'r--', 'LineWidth', 1.6);

% --- rótulos nas linhas (valores)
clabel(Cnum, hNum, 'FontSize', 8, 'Color', 'k', 'LabelSpacing', 300);
clabel(Cana, hAna, 'FontSize', 8, 'Color', 'r', 'LabelSpacing', 300);

legend([hNum hAna], {'Numérico (FVM)', 'Analítico'}, 'Location', 'best');
grid on;

% garante contornos em cima
uistack(hNum,'top'); 
uistack(hAna,'top');

%% --- Evolução da temperatura no centro e em outros pontos
figure;
plot(time, probeCenter, 'k-', 'LineWidth', 1.5); hold on;
plot(time, probeTop,    'r--', 'LineWidth', 1.2);
plot(time, probeBottom, 'b-.', 'LineWidth', 1.2);
xlabel('t (s)');
ylabel('T (^{\circ}C)');
legend('Centro','Próx. topo','Próx. base','Location','best');
grid on;
title(['Evolução da temperatura em diferentes pontos - ' methodName]);

%% --- Temperatura média do domínio vs tempo
figure;
plot(time, Tmean, 'LineWidth', 1.5);
xlabel('t (s)');
ylabel('T_{média} (^{\circ}C)');
grid on;
title(['Evolução da temperatura média do domínio - ' methodName]);

%% --- Erro relativo L2 (FVM x Analítico) vs tempo
figure;
semilogy(time, errRelL2, 'LineWidth', 1.5);
xlabel('t (s)');
ylabel('Erro relativo L_2');
grid on;
title(['Convergência da solução numérica para a solução analítica - ' methodName]);

%% --- Comparação em um corte horizontal (y = meio)
iyMid = round(mesh.Ny/2);

figure;
plot(mesh.x, T_ana(iyMid,:), 'r-', 'LineWidth', 2); hold on;
plot(mesh.x, T(iyMid,:),     'ko', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('x (m)');
ylabel('T (^{\circ}C)');
legend('Analítico','FVM','Location','best');
grid on;
title(['Corte horizontal (y = L_y/2) - ' methodName]);

%% --- Comparação em um corte vertical (x = meio)
ixMid = round(mesh.Nx/2);

figure;
plot(T_ana(:,ixMid), mesh.y, 'r-', 'LineWidth', 2); hold on;
plot(T(:,ixMid),     mesh.y, 'ko', 'MarkerSize', 5, 'LineWidth', 1);
set(gca,'YDir','normal');
xlabel('T (^{\circ}C)');
ylabel('y (m)');
legend('Analítico','FVM','Location','best');
grid on;
title(['Corte vertical (x = L_x/2) - ' methodName]);

%% ==========================================================
%% --- Gradiente da temperatura (FVM)
%% ==========================================================

% malha
dx = mesh.dx;
dy = mesh.dy;
[X, Y] = meshgrid(mesh.x, mesh.y);

% gradientes (diferenças centrais no interior)
[dTdy, dTdx] = gradient(T, dy, dx);

% módulo do gradiente
gradT_mag = sqrt(dTdx.^2 + dTdy.^2);

figure;

% fundo: |∇T| (volume a volume)
pcolor(X, Y, gradT_mag);
shading flat;              % ESSENCIAL (sem interpolação)
axis equal tight;
xlabel('x (m)');
ylabel('y (m)');

title(['Módulo do gradiente de temperatura |∇T| (' methodName ')']);

cb = colorbar;
cb.Label.String = '|∇T| (°C/m)';

grid on;


