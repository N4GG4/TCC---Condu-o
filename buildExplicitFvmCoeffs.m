function coeffsExp = buildExplicitFvmCoeffs(mesh, params, bc)
% Usa a mesma função do híbrido para montar os coeficientes
% e depois força θ = 0 (explícito) em todo o domínio.

coeffsExp = buildHybridCoeffs(mesh, params, bc);  % sua função híbrida

% explícito FVM puro: θ = 0 em todos os nós
coeffsExp.theta(:) = 0.0;

% para θ = 0, aP = aP0 (como na equação geral do TCC)
coeffsExp.aP = coeffsExp.aP0 * ones(size(coeffsExp.aP));

end
