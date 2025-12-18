function Tnew = solveHeatStepImplicitFVM(Told, sysI)
% Resolve um passo de tempo FVM totalmente implícito (θ = 1):
% A * T^{n+1} = RHS(T^n)
%
% Para θ = 1, a equação geral do TCC dá:
% RHS_p = aP0 * T_P^n + SU

    A = sysI.A;
    c = sysI.coeffs;

    Nx  = c.Nx;
    Ny  = c.Ny;
    aP0 = c.aP0;
    SU  = c.SU;

    N   = Nx * Ny;
    rhs = zeros(N,1);

    for i = 1:Nx
        for j = 1:Ny

            p = j + (i-1)*Ny;  % MESMA indexação da montagem de A

            TP_old = Told(j,i);
            rhs(p) = aP0 * TP_old + SU(j,i);
        end
    end

    % Resolve sistema
    Tnew_vec = A \ rhs;

    % Reconstrói Ny x Nx
    Tnew = reshape(Tnew_vec, Ny, Nx);
end
