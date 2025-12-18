%% Condução 1D – FVM híbrido (bordas implícitas, interior explícito)
clear; clc;

% ----- Dados do problema
L     = 0.02;
k     = 10;
rhoc  = 1e7;
alpha = k/rhoc;
T0    = 200;
TB    = 50;

% ----- Malha
N  = 20;
dx = L/N;
xc = (dx/2:dx:(L - dx/2)).';

% ----- Passo de tempo
dt = 0.5;   % tem que respeitar o limite do EXPLÍCITO, pois interior é θ=0
dt_exp = rhoc*dx^2/(2*k);
fprintf('dt = %.3g  |  dt_limite_exp = %.3g\n', dt, dt_exp);

% ----- Esquema híbrido: θ por nó
thetaP         = zeros(N,1);   % default: explícito
thetaP(1)      = 1;            % nó 1 implícito
thetaP(N)      = 1;            % nó N implícito

% ----- Coeficientes básicos
aP0 = rhoc*dx/dt;
aW  = k/dx * ones(N,1);
aE  = k/dx * ones(N,1);
SP  = zeros(N,1);
SU  = zeros(N,1);

% Ajustes de fronteira:
% Nó 1 – oeste isolado -> aW = 0
aW(1) = 0;
% Nó N – leste em T = TB  -> Sp = -2k/dx, Su = 2k/dx*TB, aE = 0
aE(N) = 0;
SP(N) = -2*k/dx;
SU(N) =  2*k/dx*TB;

% ----- Montar matriz A (N x N) usando θ local (constante no tempo)
A = zeros(N,N);
for P = 1:N
    th  = thetaP(P);
    aw  = aW(P);
    ae  = aE(P);
    Sp  = SP(P);

    % coeficiente central
    aP = aP0 + th*(aw + ae - Sp);
    A(P,P) = aP;

    % vizinho oeste
    if P > 1
        A(P,P-1) = -th*aw;
    end
    % vizinho leste
    if P < N
        A(P,P+1) = -th*ae;
    end
end

% ----- Condição inicial
T = T0*ones(N,1);

% Tempos para armazenar
times_target = [40 80 120];
keep = containers.Map('KeyType','double','ValueType','any');

% ----- Loop no tempo (híbrido)
t_end = 120;
nt    = ceil(t_end/dt);
tic;

for istep = 1:nt
    To  = T;
    RHS = zeros(N,1);

    for P = 1:N
        th  = thetaP(P);
        aw  = aW(P);
        ae  = aE(P);
        Sp  = SP(P);
        Su  = SU(P);

        % valores vizinhos em n (se não existir, tanto faz: coeficiente = 0)
        if P > 1
            Tw = To(P-1);
        else
            Tw = 0;
        end
        if P < N
            Te = To(P+1);
        else
            Te = 0;
        end

        % Equação geral θ por nó (com Sp, Su):
        RHS(P) = (1-th)*aw*Tw + (1-th)*ae*Te ...
               + (aP0 - (1-th)*(aw + ae - Sp))*To(P) + Su;
    end

    % Resolver sistema: A*T^{n+1} = RHS
    T = A \ RHS;

    tnow = istep*dt;
    if any(abs(tnow - times_target) < 1e-12)
        keep(tnow) = T;
    end
end

tempo_hibrido = toc;
fprintf('Tempo de execução (esquema híbrido) = %.6f s\n', tempo_hibrido);

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

% ====== TABELAS DE RESULTADOS (HÍBRIDO) ======
fprintf('\nResultados (°C) nos centros x (m) – esquema híbrido\n');
scale = max(abs(T0 - TB), eps);  % escala p/ erro relativo estável

for tval = times_target
    num = keep(tval);                         % numérico híbrido
    ana = TA.(sprintf('t_%d',tval));          % analítico

    % erro relativo normalizado por |T0 - TB| (%)
    err = 100 * (num - ana) ./ scale;

    fprintf('\nTempo t = %d s\n', tval);
    Ttab = table((1:N).', xc, num, ana, err, ...
        'VariableNames', {'No','x_m','Numerico','Analitico','Erro_percent'});
    disp(Ttab);
end
% =============================================

% ----- Gráfico conjunto
figure;
hold on; grid on; box on;
for tval = times_target
    plot(xc, keep(tval), 's-', 'DisplayName', ...
        sprintf('Num híbrido t = %ds', tval));
    plot(xc, TA.(sprintf('t_%d',tval)), '-', 'DisplayName', ...
        sprintf('Analítico t = %ds', tval));
end
xlabel('x (m)'); ylabel('T (°C)');
title(sprintf('FVM híbrido (bordas implícitas, interior explícito) – N=%d, Δt=%.3g', ...
      N, dt));
legend('Location','southwest');

%% --------- Gráfico de erro vs posição (RELATIVO) - 40, 80 e 120 JUNTOS
figure; hold on; grid on; box on;

for tval = times_target
    num = keep(tval);
    ana = TA.(sprintf('t_%d',tval));

    erro_rel = 100 * abs(num - ana) ./ scale;   % (%)

    plot(xc, erro_rel, 's-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('t = %ds', tval));
end

xlabel('x (m)');
ylabel('100·|T_{num} - T_{ana}| / |T_0 - T_B|  (%)');
title('Erro relativo (normalizado) vs posição — t = 40, 80 e 120 s');
legend('Location','best');

% =========================================================
% ERRO RELATIVO GLOBAL (NORMA L2) - híbrido
% normalizado por ||T_ana||_2, igual ao caso 2D
% =========================================================
fprintf('\nErro relativo global (norma L2) — híbrido:\n');
for tval = times_target
    num = keep(tval);
    ana = TA.(sprintf('t_%d',tval));

    errL2 = norm(num - ana, 2) / max(norm(ana, 2), eps);
    fprintf('t = %3d s | erro L2 relativo = %.4e\n', tval, errL2);
end

% =========================================================
% RODADA COM dt/2 PARA ESTIMAR ERRO TEMPORAL (RICHARDSON)
% =========================================================
dt2   = dt/2;
aP0_2 = rhoc*dx/dt2;

% --- Recria A2 (depende de aP0)
A2 = zeros(N,N);
for P = 1:N
    th  = thetaP(P);
    aw  = aW(P);
    ae  = aE(P);
    Sp  = SP(P);

    aP = aP0_2 + th*(aw + ae - Sp);
    A2(P,P) = aP;

    if P > 1, A2(P,P-1) = -th*aw; end
    if P < N, A2(P,P+1) = -th*ae; end
end

% --- Integra com dt/2
T2    = T0*ones(N,1);
nt2   = ceil(t_end/dt2);
keep2 = containers.Map('KeyType','double','ValueType','any');

for istep = 1:nt2
    To2  = T2;
    RHS2 = zeros(N,1);

    for P = 1:N
        th  = thetaP(P);
        aw  = aW(P);
        ae  = aE(P);
        Sp  = SP(P);
        Su  = SU(P);

        if P > 1, Tw = To2(P-1); else, Tw = 0; end
        if P < N, Te = To2(P+1); else, Te = 0; end

        RHS2(P) = (1-th)*aw*Tw + (1-th)*ae*Te ...
                + (aP0_2 - (1-th)*(aw + ae - Sp))*To2(P) + Su;
    end

    T2 = A2 \ RHS2;

    tnow2 = istep*dt2;

    % salva de forma robusta quando ultrapassar o tempo-alvo
    for i = 1:numel(times_target)
        tt = times_target(i);
        if ~isKey(keep2, tt) && (tnow2 >= tt)
            keep2(tt) = T2;
        end
    end
end

% =========================================================
% GRÁFICO: erro relativo vs posição COM BARRAS (Richardson)
% =========================================================
p   = 1;          % híbrido (como implementado) é 1ª ordem no tempo
fac = 2^p - 1;

figure; hold on; grid on; box on;

for tval = times_target
    num  = keep(tval);      % dt
    num2 = keep2(tval);     % dt/2
    ana  = TA.(sprintf('t_%d',tval));

    % erro relativo principal (normalizado por |T0-TB|)
    err_rel = 100 * abs(num - ana) ./ scale;

    % estimativa do erro temporal (Richardson) em %
    err_time_est = abs(num2 - num) ./ fac;      % °C
    bar_rel      = 100 * err_time_est ./ scale; % %

    errorbar(xc, err_rel, bar_rel, 'o-', ...
        'LineWidth', 1.2, 'CapSize', 6, ...
        'DisplayName', sprintf('t = %ds', tval));
end

xlabel('x (m)');
ylabel('Erro relativo (%)');
title('Erro relativo vs posição (com barras por refino em \Deltat) — híbrido');
legend('Location','best');



