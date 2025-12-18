function Tnew = solveHeatStepExplicitFVM(Told, coeffs)
% Esquema FVM 2D totalmente explícito (θ = 0 em todos os nós)
% Mesma formulação do TCC, só que sem resolver sistema (um passo direto).

Nx    = coeffs.Nx;
Ny    = coeffs.Ny;
aE    = coeffs.aE;
aW    = coeffs.aW;
aN    = coeffs.aN;
aS    = coeffs.aS;
SP    = coeffs.SP;
SU    = coeffs.SU;
aP0   = coeffs.aP0;   % escalar
aP    = coeffs.aP;    % aqui vale aP0

Tnew = Told;  % só para ter o tamanho certo

for j = 1:Ny
    for i = 1:Nx

        aE_ = aE(j,i);
        aW_ = aW(j,i);
        aN_ = aN(j,i);
        aS_ = aS(j,i);
        SP_ = SP(j,i);
        SU_ = SU(j,i);
        aP_ = aP(j,i);   % = aP0

        % vizinhos em n (explícito)
        if i < Nx
            TE_old = Told(j,i+1);
        else
            TE_old = 0;
        end

        if i > 1
            TW_old = Told(j,i-1);
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

        sum_nb   = aE_*TE_old + aW_*TW_old + aN_*TN_old + aS_*TS_old;
        sum_a_nb = (aE_ + aW_ + aN_ + aS_);

        % termo de armazenamento para θ = 0
        storage = (aP0 - (sum_a_nb - SP_)) * Told(j,i);

        RHS = sum_nb + storage + SU_;

        Tnew(j,i) = RHS / aP_;
    end
end

end