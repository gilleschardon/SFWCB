function [epos, eamps] = compute_errors(Xest, Xtrue, Aest, Atrue)

dists = (Xest(:,1) - Xtrue(:,1)').^2 + (Xest(:,2) - Xtrue(:,2)').^2 + (Xest(:,3) - Xtrue(:,3)').^2;

M = matchpairs(dists, 1e10);

epos = 0;
eamps = 0;

for u = 1:size(M, 1)
    epos = epos + norm(Xest(M(u, 1), :) - Xtrue(M(u, 2), :))^2;
    eamps = eamps + norm(abs(Aest(M(u, 1), :)) - abs(Atrue(M(u, 2), :)), 'fro')^2 / size(Aest, 2);
end



end