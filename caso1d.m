%% Condução 1D em FVM com θ-esquema geral
% theta = 0   -> Explícito (FTCS)
% theta = 0.5 -> Crank–Nicolson
% theta = 1   -> Implícito (BTCS)

clear; clc;

% ----- Dados do problema
L     = 0.02;           % m  (espessura)
k     = 10;             % W/m.K
rhoc  = 1e7;            % J/m^3.K
alpha = k/rhoc;         % m^2/s
T0    = 200;            % °C (inicial)
TB    = 50;              % °C (fronteira leste após t=0)

% ----- Escolha do esquema temporal
theta = 1;            % 0 = explícito, 0.5 = CN, 1 = implícito

% ----- Malha 1D
N  = 20;                                % nº de volumes
dx = L/N;
xc = (dx/2:dx:(L - dx/2)).';           % centros

% ----- Passo de tempo
dt = 0.5;                               % s (mude à vontade)
% ----- Critérios de estabilidade / consistência em função de theta
Fo = k*dt/(rhoc*dx^2);    % número de Fourier

if abs(theta) < 1e-12
    % ----- theta = 0  -> explícito
    dt_exp = rhoc*dx^2/(2*k);    % condição clássica de estabilidade
    Fo_exp = k*dt_exp/(rhoc*dx^2);  % = 0.5

    fprintf('Esquema EXPLÍCITO (theta = 0)\n');
    fprintf('Fo   = %.3g  |  Fo_lim = %.3g\n', Fo, Fo_exp);
    fprintf('Δt   = %.3g s | Δt_lim(explicito) = %.3g s\n', dt, dt_exp);

    if dt > dt_exp
        warning('Δt > Δt_lim do esquema explícito: solução tende a divergir.');
    end

elseif abs(theta - 0.5) < 1e-12
    % ----- theta = 0.5 -> Crank–Nicolson
    dt_CN  = rhoc*dx^2/k;       % sua Eq. (3.71)
    Fo_CN  = k*dt_CN/(rhoc*dx^2);  % = 1

    fprintf('Esquema CRANK–NICOLSON (theta = 0.5)\n');
    fprintf('Fo   = %.3g\n', Fo);
    fprintf('Δt   = %.3g s | Δt_max p/ coeficientes positivos (Eq. 3.71) = %.3g s\n', ...
            dt, dt_CN);

    if dt > dt_CN
        warning(['Δt maior que o recomendado pela Eq. (3.71): ' ...
                 'coeficientes podem ficar não físicos (negativos).']);
    end

else
    % ----- theta > 0.5 -> implícito genérico (inclui theta = 1)
    fprintf('Esquema IMPLÍCITO (theta = %.2f)\n', theta);
    fprintf('Fo   = %.3g (estável para qualquer Δt neste problema difusivo).\n', Fo);
end

t_end = 120;
nt    = ceil(t_end/dt);

% ----- Condição inicial
T = T0*ones(N,1);

% Tempos para armazenar a solução
times_target = [40 80 120];
keep = containers.Map('KeyType','double','ValueType','any');

% ----- Coeficientes básicos
aP0 = rhoc*dx/dt;   % termo transiente
aW  = k/dx;
aE  = k/dx;
G   = 2*k/dx;       % coeficiente do Dirichlet na face leste

% =========================================================
% 1) Montar o sistema IMPLÍCITO (theta = 1):
%     A_imp * T^{n+1} = aP0 * T^n + b_imp
% =========================================================
A_imp = zeros(N,N);
b_imp = zeros(N,1);

% Nó 1 (oeste isolado, q''=0 => Neumann homogênea)
A_imp(1,1) = aP0 + aE;
A_imp(1,2) = -aE;
% b_imp(1) = 0;

% Nós internos 2..N-1
for P = 2:N-1
    A_imp(P,P-1) = -aW;
    A_imp(P,P)   = aP0 + aW + aE;
    A_imp(P,P+1) = -aE;
end

% Nó N (leste com T = TB; Dirichlet por fonte equivalente)
A_imp(N,N-1) = -aW;
A_imp(N,N)   = aP0 + aW + G;    % G = 2k/dx
b_imp(N)     = G*TB;            % Su = 2k/dx * TB

% =========================================================
% 2) Construir A_theta e M_theta a partir de A_imp (θ-esquema)
% =========================================================
I      = eye(N);
A_theta = (1-theta)*aP0*I + theta*A_imp;
M_theta = (2-theta)*aP0*I - (1-theta)*A_imp;
b_theta = b_imp;   % independe de theta

% ----- Loop no tempo
tic;
for istep = 1:nt
    To  = T;                         % T^n
    RHS = M_theta * To + b_theta;   % lado direito

    % Resolver sistema A_theta * T^{n+1} = RHS
    T = A_theta \ RHS;

    tnow = istep*dt;
    if any(abs(tnow - times_target) < 1e-12)
        keep(tnow) = T;
    end
end
tempo_theta = toc;  % termina a contagem
fprintf('Tempo de execução (esquema theta) = %.6f s\n', tempo_theta);
% =========================================================
% 3) Solução analítica (para comparação) - COM TB
% =========================================================
analytical = @(x,t) TB + (T0 - TB)*(4/pi)*sum(arrayfun(@(n) ...
    ((-1)^(n+1))/(2*n-1)*exp(-alpha*((2*n-1)*pi/(2*L))^2*t).* ...
    cos(((2*n-1)*pi/(2*L))*x), 1:200));


TA = struct();
for tval = times_target
    TA.(sprintf('t_%d',tval)) = arrayfun(@(x) analytical(x,tval), xc);
end

% ----- Tabelas (numérico x analítico)
fprintf('\nResultados (°C) nos centros x (m)\n');
for tval = times_target
    num = keep(tval);
    ana = TA.(sprintf('t_%d',tval));
    err = 100 * abs(num - ana) ./ max(abs(T0 - TB), eps);
    fprintf('\nTempo t = %d s\n', tval);
    Ttab = table((1:N).', xc, num, ana, err, ...
        'VariableNames', {'No','x_m','Numerico','Analitico','Erro_percent'});
    disp(Ttab);
end

% =========================================================
% Erro global no console (norma L2 relativa)
% =========================================================
fprintf('\nErro relativo global (norma L2):\n');

for tval = times_target
    num = keep(tval);
    ana = TA.(sprintf('t_%d',tval));

    errL2 = norm(num - ana, 2) / norm(ana, 2);

    fprintf('t = %3d s | erro L2 relativo = %.4e\n', tval, errL2);
end


% ----- Gráfico conjunto
figure;
hold on; grid on; box on;
for tval = times_target
    plot(xc, keep(tval), 's-', 'DisplayName', ...
        sprintf('Num (\\theta=%.2f) t = %ds', theta, tval));
    plot(xc, TA.(sprintf('t_%d',tval)), '-', 'DisplayName', ...
        sprintf('Analítico t = %ds', tval));
end
xlabel('x (m)'); ylabel('T (°C)');
title(sprintf('FVM 1D – esquema \\theta (N = %d, \\theta = %.2f, \\Deltat = %.3g s)', ...
      N, theta, dt));
legend('Location','southwest');

%% --------- Gráfico de erro vs posição (NORMALIZADO por |T0-TB|) - juntos
figure; hold on; grid on; box on;

scale = max(abs(T0 - TB), eps);

for tval = times_target
    num = keep(tval);
    ana = TA.(sprintf('t_%d',tval));

    erro_rel = 100 * abs(num - ana) ./ scale;

    plot(xc, erro_rel, 's-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('t = %ds', tval));
end

xlabel('x (m)');
ylabel('100·|T_{num} - T_{ana}| / |T_0 - T_B|  (%)');
title('Erro normalizado vs posição — t = 40, 80 e 120 s');
legend('Location','best');

%% --------- Gráfico de erro vs posição (COM BARRAS por refino em dt) - juntos
% Ideia: estimar erro temporal comparando solução com dt e dt/2
% Barra = |T(dt/2) - T(dt)| / (2^p - 1), com p = ordem temporal do esquema

% --- ordem temporal do theta-esquema
if abs(theta - 0.5) < 1e-12
    p = 2;   % Crank–Nicolson
else
    p = 1;   % explícito ou implícito de 1a ordem
end
fac = 2^p - 1;

% --- rodar novamente com dt/2 (mesma física e malha)
dt2 = dt/2;

t_end2 = 120;
nt2    = ceil(t_end2/dt2);

% (recalcular coeficientes dependentes de dt2)
aP0_2 = rhoc*dx/dt2;

% montar A_imp e b_imp (iguais, só muda aP0)
A_imp2 = zeros(N,N);
b_imp2 = zeros(N,1);

A_imp2(1,1) = aP0_2 + aE;
A_imp2(1,2) = -aE;

for P = 2:N-1
    A_imp2(P,P-1) = -aW;
    A_imp2(P,P)   = aP0_2 + aW + aE;
    A_imp2(P,P+1) = -aE;
end

A_imp2(N,N-1) = -aW;
A_imp2(N,N)   = aP0_2 + aW + G;
b_imp2(N)     = G*TB;

I = eye(N);
A_theta2 = (1-theta)*aP0_2*I + theta*A_imp2;
M_theta2 = (2-theta)*aP0_2*I - (1-theta)*A_imp2;
b_theta2 = b_imp2;

% condição inicial
T2 = T0*ones(N,1);

% guardar tempos-alvo para dt/2
keep2 = containers.Map('KeyType','double','ValueType','any');

for istep = 1:nt2
    To2  = T2;
    RHS2 = M_theta2 * To2 + b_theta2;
    T2   = A_theta2 \ RHS2;

    tnow2 = istep*dt2;

    % salva quando passar do alvo (robusto)
    for i = 1:numel(times_target)
        tt = times_target(i);
        if ~isKey(keep2, tt) && (tnow2 >= tt)
            keep2(tt) = T2;
        end
    end
end

% --- agora plota erro relativo com barras
figure; hold on; grid on; box on;

scale = max(abs(T0 - TB), eps);

for tval = times_target
    num  = keep(tval);     % solução com dt (você já tem)
    num2 = keep2(tval);    % solução com dt/2
    ana  = TA.(sprintf('t_%d',tval));

    % erro relativo principal (%)
    err_rel = 100 * abs(num - ana) ./ scale;

    % barra: estimativa do erro temporal (Richardson) -> em %
    err_time_est = abs(num2 - num) ./ fac;   % em °C
    bar_rel      = 100 * err_time_est ./ scale;

    errorbar(xc, err_rel, bar_rel, 'o-', ...
        'LineWidth', 1.2, 'CapSize', 6, ...
        'DisplayName', sprintf('t = %ds', tval));
end

xlabel('x (m)');
ylabel('Erro relativo (%)');
title(sprintf('Erro relativo vs posição (com barras por refino em \\Deltat) — \\theta=%.2f', theta));
legend('Location','best');




