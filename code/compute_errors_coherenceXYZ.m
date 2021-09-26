function [eposX, eposY, eposZ, eamps] = compute_errors_coherenceXYZ(Xest, Xtrue, Aest, Atrue, Pmic, k)

D = dictionary(Pmic, Xest, k);
norms = sqrt(sum(abs(D).^2, 1))';

Dnorm = D ./ norms';

coh = abs(Dnorm'*Dnorm) - eye(size(Dnorm, 2));




[m, u] = max(coh);
[m, v] = max(m);

u = u(v);
thres = 0.98;

if min(min(coh + eye(size(Dnorm, 2)))) > thres
    Amoy = sum(Aest(:))/3;
    Aest = Amoy * ones(size(Aest));    
elseif m > thres
    Amoy = (Aest(u) + Aest(v))/2;
    Aest(u) = Amoy;
    Aest(v) = Amoy;
end
% compute position and amplitude errors by matching sources

dists = (Xest(:,1) - Xtrue(:,1)').^2 + (Xest(:,2) - Xtrue(:,2)').^2 + (Xest(:,3) - Xtrue(:,3)').^2;

M = matchpairs(dists, 1e10);

eposX = 0;
eposY = 0;
eposZ = 0;

eamps = 0;

for u = 1:size(M, 1)
    eposX = eposX + norm(Xest(M(u, 1), 1) - Xtrue(M(u, 2), 1))^2;
    eposY = eposY + norm(Xest(M(u, 1), 2) - Xtrue(M(u, 2), 2))^2;
    eposZ = eposZ + norm(Xest(M(u, 1), 3) - Xtrue(M(u, 2), 3))^2;

    eamps = eamps + norm(abs(Aest(M(u, 1), :)) - abs(Atrue(M(u, 2), :)), 'fro')^2 / size(Aest, 2);
end



end