function Tnew = solveHeatStepHybridMatrix(Told, sysH)
% A * T^{n+1} = RHS(T^n) para o esquema híbrido θ-local

    A      = sysH.A;
    c      = sysH.coeffs;

    Nx    = c.Nx;
    Ny    = c.Ny;
    aE    = c.aE;
    aW    = c.aW;
    aN    = c.aN;
    aS    = c.aS;
    SP    = c.SP;
    SU    = c.SU;
    aP0   = c.aP0;
    theta = c.theta;

    rhs = zeros(Nx*Ny,1);

    for i = 1:Nx
        for j = 1:Ny

            p = j + (i-1)*Ny;   % MESMA indexação da montagem da matriz

            th  = theta(j,i);
            aE_ = aE(j,i);
            aW_ = aW(j,i);
            aN_ = aN(j,i);
            aS_ = aS(j,i);
            SP_ = SP(j,i);
            SU_ = SU(j,i);

            TP_old = Told(j,i);

            % vizinhos em n
            if i < Nx
                TE_old = Told(j,  i+1);
            else
                TE_old = 0;
            end
            if i > 1
                TW_old = Told(j,  i-1);
            else
                TW_old = 0;
            end
            if j < Ny
                TN_old = Told(j+1,i);
            else
                TN_old = 0;
            end
            if j > 1
                TS_old = Told(j-1,i);
            else
                TS_old = 0;
            end

            sum_a_nb   = aE_ + aW_ + aN_ + aS_;
            % parte com vizinhos em n (explícita)
            sum_nb_old = (1-th)*( aE_*TE_old + aW_*TW_old + aN_*TN_old + aS_*TS_old );
            % termo de armazenamento geral
            storage    = (aP0 - (1-th)*(sum_a_nb - SP_)) * TP_old;

            rhs(p) = sum_nb_old + storage + SU_;
        end
    end

    % Resolve sistema híbrido
    Tnew_vec = A \ rhs;

    % Reconstrói Ny x Nx (compatível com p = j + (i-1)*Ny)
    Tnew = reshape(Tnew_vec, Ny, Nx);
end

