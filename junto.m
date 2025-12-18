%% Condução 1D – 3 métodos (Explícito, Implícito, Híbrido)
% 3 FIGURAS TOTAIS:
% (1) Perfis: todos os métodos + analítico, para t=40/80/120 no MESMO gráfico
% (2) Erro relativo vs x: todos os métodos e tempos no MESMO gráfico
% (3) Erro relativo com barras (Richardson dt/2): tudo no MESMO gráfico

clear; clc;

%% ---- Dados do problema
L     = 0.02;
k     = 10;
rhoc  = 1e7;
alpha = k/rhoc;

T0    = 200;
TB    = 50;

%% ---- Malha
N  = 20;
dx = L/N;
xc = (dx/2:dx:(L - dx/2)).';

%% ---- Tempo
dt    = 0.5;
t_end = 120;
times_target = [40 80 120];

%% ---- Métodos (só 3)
methods = { ...
    struct('name','Explícito', 'type','theta',  'theta',0.0), ...
    struct('name','Implícito', 'type','theta',  'theta',1.0), ...
    struct('name','Híbrido',   'type','hybrid', 'theta',NaN) ...
};

%% ---- Analítico (série truncada)
analytical = @(x,t) TB + (T0 - TB)*(4/pi)*sum(arrayfun(@(n) ...
    ((-1)^(n+1))/(2*n-1)*exp(-alpha*((2*n-1)*pi/(2*L))^2*t).* ...
    cos(((2*n-1)*pi/(2*L))*x), 1:200));

TA = struct();
for tval = times_target
    TA.(sprintf('t_%d',tval)) = arrayfun(@(x) analytical(x,tval), xc);
end

scale = max(abs(T0 - TB), eps);

%% ---- Funções de execução (dt) e (dt/2)
run_theta_dt   = @(theta) local_run_theta(N, dx, dt,  t_end, times_target, k, rhoc, T0, TB, theta);
run_hybrid_dt  = @()      local_run_hybrid(N, dx, dt,  t_end, times_target, k, rhoc, T0, TB);

run_theta_dt2  = @(theta) local_run_theta(N, dx, dt/2,t_end, times_target, k, rhoc, T0, TB, theta);
run_hybrid_dt2 = @()      local_run_hybrid(N, dx, dt/2,t_end, times_target, k, rhoc, T0, TB);

%% ---- Estabilidade explícito / híbrido (interior explícito)
dt_lim_exp = rhoc*dx^2/(2*k);

%% ---- Armazena resultados: results.(methodTag).(tTag)
results_dt  = struct();
results_dt2 = struct();
ran = false(size(methods));

for m = 1:numel(methods)
    mt  = methods{m};
    tag = matlab.lang.makeValidName(mt.name);

    % ---- roda em dt
    if mt.type == "theta"
        if abs(mt.theta) < 1e-12 && dt > dt_lim_exp
            fprintf('[AVISO] Pulando %s: dt=%.3g > dt_lim=%.3g (explícito instável)\n', mt.name, dt, dt_lim_exp);
            continue;
        end
        if abs(mt.theta - 1) < 1e-12
            % implícito: ok
        end
        keep_dt = run_theta_dt(mt.theta);
        keep_dt2 = run_theta_dt2(mt.theta);
        p = 1; % explícito/implícito BTCS -> 1ª ordem no tempo
    else
        if dt > dt_lim_exp
            fprintf('[AVISO] %s: dt=%.3g > dt_lim=%.3g (interior explícito instável)\n', mt.name, dt, dt_lim_exp);
            % Se quiser impedir rodar, descomente:
            % continue;
        end
        keep_dt  = run_hybrid_dt();
        keep_dt2 = run_hybrid_dt2();
        p = 1;
    end

    ran(m) = true;

    for tval = times_target
        results_dt.(tag).(sprintf('t_%d',tval))  = keep_dt(tval);
        results_dt2.(tag).(sprintf('t_%d',tval)) = keep_dt2(tval);
    end

    % guarda p/ barras (Richardson)
    results_dt.(tag).p = p;

    fprintf('OK: %s\n', mt.name);
end

%% ============================================================
%  (1) PERFIS: tudo em um único gráfico
% ============================================================
figure; hold on; grid on; box on;

% analítico (3 tempos)
for tval = times_target
    plot(xc, TA.(sprintf('t_%d',tval)), '-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Analítico (t=%ds)', tval));
end

% numéricos: 3 métodos x 3 tempos
for m = 1:numel(methods)
    if ~ran(m), continue; end
    mt  = methods{m};
    tag = matlab.lang.makeValidName(mt.name);

    for tval = times_target
        Tm = results_dt.(tag).(sprintf('t_%d',tval));
        plot(xc, Tm, 'o-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s (t=%ds)', mt.name, tval));
    end
end

xlabel('x (m)');
ylabel('T (°C)');
title(sprintf('Perfis de temperatura – N=%d, dt=%.3g (t=40,80,120 s)', N, dt));
legend('Location','bestoutside');

%% ============================================================
%  (2) ERRO RELATIVO vs x: tudo em um único gráfico
%      100*|Tnum - Tana|/|T0-TB|
% ============================================================
figure; hold on; grid on; box on;

for m = 1:numel(methods)
    if ~ran(m), continue; end
    mt  = methods{m};
    tag = matlab.lang.makeValidName(mt.name);

    for tval = times_target
        num = results_dt.(tag).(sprintf('t_%d',tval));
        ana = TA.(sprintf('t_%d',tval));

        err_rel = 100*abs(num - ana)./scale;

        plot(xc, err_rel, 's-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s (t=%ds)', mt.name, tval));
    end
end

xlabel('x (m)');
ylabel('Erro relativo (%)');
title(sprintf('Erro relativo vs posição – N=%d, dt=%.3g (t=40,80,120 s)', N, dt));
legend('Location','bestoutside');

%% ============================================================
%  (3) ERRO COM BARRAS (Richardson dt/2): tudo em um único gráfico
%      barra ≈ |T(dt/2) - T(dt)|/(2^p - 1)
% ============================================================
figure; hold on; grid on; box on;

for m = 1:numel(methods)
    if ~ran(m), continue; end
    mt  = methods{m};
    tag = matlab.lang.makeValidName(mt.name);

    p = results_dt.(tag).p;
    fac = 2^p - 1;

    for tval = times_target
        num  = results_dt.(tag).(sprintf('t_%d',tval));
        num2 = results_dt2.(tag).(sprintf('t_%d',tval));
        ana  = TA.(sprintf('t_%d',tval));

        err_rel = 100*abs(num - ana)./scale;

        err_time_est = abs(num2 - num)./fac;     % °C (estimativa do erro temporal)
        bar_rel      = 100*err_time_est./scale;  % %

        errorbar(xc, err_rel, bar_rel, 'o-', ...
            'LineWidth', 1.1, 'CapSize', 6, ...
            'DisplayName', sprintf('%s (t=%ds)', mt.name, tval));
    end
end

xlabel('x (m)');
ylabel('Erro relativo (%)');
title(sprintf('Erro relativo com barras (Richardson em dt) – N=%d, dt=%.3g', N, dt));
legend('Location','bestoutside');

%% ============================================================
%  =============== FUNÇÕES AUXILIARES =========================
% ============================================================
function keep = local_run_theta(N, dx, dt, t_end, times_target, k, rhoc, T0, TB, theta)
    aP0 = rhoc*dx/dt;
    aW  = k/dx;
    aE  = k/dx;
    G   = 2*k/dx;

    A_imp = zeros(N,N);
    b_imp = zeros(N,1);

    % Nó 1: Neumann homogêneo (isolado)
    A_imp(1,1) = aP0 + aE;
    A_imp(1,2) = -aE;

    % Internos
    for P = 2:N-1
        A_imp(P,P-1) = -aW;
        A_imp(P,P)   = aP0 + aW + aE;
        A_imp(P,P+1) = -aE;
    end

    % Nó N: Dirichlet via fonte equivalente
    A_imp(N,N-1) = -aW;
    A_imp(N,N)   = aP0 + aW + G;
    b_imp(N)     = G*TB;

    I      = eye(N);
    A_main = (1-theta)*aP0*I + theta*A_imp;
    M_main = (2-theta)*aP0*I - (1-theta)*A_imp;

    T = T0*ones(N,1);
    keep = containers.Map('KeyType','double','ValueType','any');

    nt = ceil(t_end/dt);
    for istep = 1:nt
        To  = T;
        RHS = M_main*To + b_imp;
        T   = A_main \ RHS;

        tnow = istep*dt;
        if any(abs(tnow - times_target) < 1e-12)
            keep(tnow) = T;
        end
    end
end

function keep = local_run_hybrid(N, dx, dt, t_end, times_target, k, rhoc, T0, TB)
    % bordas implícitas (θ=1) e interior explícito (θ=0)
    thetaP    = zeros(N,1);
    thetaP(1) = 1;
    thetaP(N) = 1;

    aP0 = rhoc*dx/dt;
    aW  = (k/dx)*ones(N,1);
    aE  = (k/dx)*ones(N,1);
    SP  = zeros(N,1);
    SU  = zeros(N,1);

    % Fronteiras
    aW(1) = 0;                 % oeste isolado
    aE(N) = 0;                 % leste Dirichlet por fonte eq.
    SP(N) = -2*k/dx;
    SU(N) =  2*k/dx*TB;

    % Matriz A constante
    A = zeros(N,N);
    for P = 1:N
        th  = thetaP(P);
        aw  = aW(P);
        ae  = aE(P);
        Sp  = SP(P);

        aP = aP0 + th*(aw + ae - Sp);
        A(P,P) = aP;

        if P > 1, A(P,P-1) = -th*aw; end
        if P < N, A(P,P+1) = -th*ae; end
    end

    T = T0*ones(N,1);
    keep = containers.Map('KeyType','double','ValueType','any');

    nt = ceil(t_end/dt);
    for istep = 1:nt
        To  = T;
        RHS = zeros(N,1);

        for P = 1:N
            th  = thetaP(P);
            aw  = aW(P);
            ae  = aE(P);
            Sp  = SP(P);
            Su  = SU(P);

            Tw = 0; Te = 0;
            if P > 1, Tw = To(P-1); end
            if P < N, Te = To(P+1); end

            RHS(P) = (1-th)*aw*Tw + (1-th)*ae*Te ...
                   + (aP0 - (1-th)*(aw + ae - Sp))*To(P) + Su;
        end

        T = A \ RHS;

        tnow = istep*dt;
        if any(abs(tnow - times_target) < 1e-12)
            keep(tnow) = T;
        end
    end
end
