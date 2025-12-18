function sysH = buildHybridSystemMatrix(mesh, params, bc)
% Monta sistema linear do esquema híbrido 2D (θ-local) em forma matricial:
% A * T^{n+1} = RHS(T^n)

    % Coeficientes FVM + theta (do seu TCC)
    coeffs = buildHybridCoeffs(mesh, params, bc);

    Nx = coeffs.Nx;
    Ny = coeffs.Ny;
    N  = Nx * Ny;

    aE    = coeffs.aE;
    aW    = coeffs.aW;
    aN    = coeffs.aN;
    aS    = coeffs.aS;
    SP    = coeffs.SP;
    SU    = coeffs.SU;
    aP0   = coeffs.aP0;
    theta = coeffs.theta;

    % aP segundo a fórmula do TCC
    aP = theta .* (aE + aW + aN + aS - SP) + aP0;
    coeffs.aP = aP;   % guarda também

    % Vamos montar A como esparsa (cada linha tem até 5 entradas)
    ii = zeros(5*N,1);
    jj = zeros(5*N,1);
    vv = zeros(5*N,1);
    k  = 0;

    for i = 1:Nx
        for j = 1:Ny

            % índice global em ordem COLUNA-major (p = j + (i-1)*Ny)
            p = j + (i-1)*Ny;

            th  = theta(j,i);
            aE_ = aE(j,i);
            aW_ = aW(j,i);
            aN_ = aN(j,i);
            aS_ = aS(j,i);
            aP_ = aP(j,i);

            % Diagonal sempre entra
            k = k + 1;
            ii(k) = p;
            jj(k) = p;
            vv(k) = aP_;

            % Se θ = 0, célula é totalmente explícita:
            % não há vizinhos em T^{n+1} -> só a diagonal
            if th == 0
                continue;
            end

            % Vizinho Leste (mesma linha, coluna i+1)
            if i < Nx
                pE = j + (i    )*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pE;
                vv(k) = -th * aE_;
            end

            % Vizinho Oeste (mesma linha, coluna i-1)
            if i > 1
                pW = j + (i-2)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pW;
                vv(k) = -th * aW_;
            end

            % Vizinho Norte (linha j+1)
            if j < Ny
                pN = (j+1) + (i-1)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pN;
                vv(k) = -th * aN_;
            end

            % Vizinho Sul (linha j-1)
            if j > 1
                pS = (j-1) + (i-1)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pS;
                vv(k) = -th * aS_;
            end
        end
    end

    % Monta matriz esparsa usando só as entradas preenchidas
    A = sparse(ii(1:k), jj(1:k), vv(1:k), N, N);

    % Guarda tudo que vamos precisar depois
    sysH.A      = A;
    sysH.coeffs = coeffs;
end

