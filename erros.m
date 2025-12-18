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

%% --- Condição inicial
T0 = 25.0;

%% --- Condições de contorno
T_south = 50;    
T_north = 75;
T_west  = 100;
T_east  = 25;

bc.left.type   = 'D';  bc.left.value   = T_west;
bc.right.type  = 'D';  bc.right.value  = T_east;
bc.bottom.type = 'D';  bc.bottom.value = T_south;
bc.top.type    = 'D';  bc.top.value    = T_north;

%% ==========================================================
%% --- Roda os 3 métodos e guarda erros
%% ==========================================================
schemes = ["explicit","implicit","hybrid"];

results = struct();

for s = 1:numel(schemes)

    scheme = schemes(s);

    % --- flags
    params.useImplicit = (scheme == "implicit");

    switch scheme
        case "implicit"
            methodName = 'Implícito';
        case "explicit"
            methodName = 'Explícito FVM';
        case "hybrid"
            methodName = 'Híbrido';
    end

    %% --- Configuração do passo de tempo dependendo do esquema
    % Ajuste aqui como você quiser:
    if scheme == "explicit"
        params.CFL = 1.0016;

    elseif scheme == "hybrid"
        params.CFL = 1.0016;

    elseif scheme == "implicit"
        % para o implícito, se quiser menos passos, aumente CFLimp:
        params.CFLimp = 1.0016;   % ex: 10, 50, 100 ...
        % ou use params.dt / params.nSteps se seu buildMesh suportar
    end

    %% --- Monta malha e obtém dt e nSteps
    [mesh, dt, nSteps] = buildMesh(params);
    mesh.dt = dt;

    fprintf('--- %s --- dx=%.3e dy=%.3e dt=%.3e nSteps=%d\n',...
        scheme, mesh.dx, mesh.dy, dt, nSteps);

    %% --- Inicializa campo
    T = T0*ones(mesh.Ny, mesh.Nx);
    T = applyBC(T, mesh, params, bc);

    %% --- Analítico estacionário (para erro)
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

    %% --- Vetores de erro
    time     = (0:nSteps)' * dt;
    errRelL2 = zeros(nSteps+1, 1);
    errRelL2(1) = norm(T(:) - T_ana(:), 2) / normAnaL2;

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

        errRelL2(n+1) = norm(T(:) - T_ana(:), 2) / normAnaL2;
    end
    tempoExecucao = toc;

    fprintf('   Tempo = %.3f s | Erro final = %.3e\n\n', tempoExecucao, errRelL2(end));

    %% --- Salva resultados
    results.(scheme).time = time;
    results.(scheme).err  = errRelL2;
    results.(scheme).name = methodName;
end

%% ==========================================================
%% --- Gráfico: erro dos 3 métodos juntos
%% ==========================================================
figure; hold on;

semilogy(results.explicit.time, results.explicit.err, 'LineWidth', 1.8);
semilogy(results.implicit.time, results.implicit.err, 'LineWidth', 1.8);
semilogy(results.hybrid.time,   results.hybrid.err,   'LineWidth', 1.8);

xlabel('t (s)');
ylabel('Erro relativo L_2');
grid on;

legend(results.explicit.name, results.implicit.name, results.hybrid.name, ...
    'Location','best');

title('Convergência do erro relativo L_2 — comparação entre esquemas FVM');
