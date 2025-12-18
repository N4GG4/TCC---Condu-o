function sysI = buildImplicitSystemFVM(mesh, params, bc)
% Monta sistema linear do esquema FVM totalmente implícito (θ = 1)
% A * T^{n+1} = RHS(T^n)

    % Coeficientes FVM gerais (do TCC)
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

    % Implícito total: θ = 1 em todos os nós
    theta        = ones(Ny, Nx);
    coeffs.theta = theta;

    % aP segundo a fórmula do TCC com θ = 1
    aP = (aE + aW + aN + aS - SP) + aP0;
    coeffs.aP = aP;

    % Monta A esparsa (mesma indexação coluna-major que usamos no híbrido)
    ii = zeros(5*N,1);
    jj = zeros(5*N,1);
    vv = zeros(5*N,1);
    k  = 0;

    for i = 1:Nx
        for j = 1:Ny

            % índice global em ordem COLUNA-major
            % p = j + (i-1)*Ny
            p = j + (i-1)*Ny;

            aE_ = aE(j,i);
            aW_ = aW(j,i);
            aN_ = aN(j,i);
            aS_ = aS(j,i);
            aP_ = aP(j,i);

            % Diagonal
            k = k + 1;
            ii(k) = p;
            jj(k) = p;
            vv(k) = aP_;

            % Vizinho Leste
            if i < Nx
                pE = j + (i    )*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pE;
                vv(k) = -aE_;
            end

            % Vizinho Oeste
            if i > 1
                pW = j + (i-2)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pW;
                vv(k) = -aW_;
            end

            % Vizinho Norte
            if j < Ny
                pN = (j+1) + (i-1)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pN;
                vv(k) = -aN_;
            end

            % Vizinho Sul
            if j > 1
                pS = (j-1) + (i-1)*Ny;
                k = k + 1;
                ii(k) = p;
                jj(k) = pS;
                vv(k) = -aS_;
            end

        end
    end

    A = sparse(ii(1:k), jj(1:k), vv(1:k), N, N);

    % guarda o que vamos precisar
    sysI.A      = A;
    sysI.coeffs = coeffs;
end
