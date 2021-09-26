function [epos, eamps] = compute_errors_music(Xest, Xtrue, Y, Atrue, Dic)

% compute position and amplitude errors by matching sources

dists = (Xest(:,1) - Xtrue(:,1)').^2 + (Xest(:,2) - Xtrue(:,2)').^2 + (Xest(:,3) - Xtrue(:,3)').^2;

M = matchpairs(dists, 1e10);

epos = 0;
eamps = 0;

Dloc = [];

for u = 1:size(M, 1)
    epos = epos + norm(Xest(M(u, 1), :) - Xtrue(M(u, 2), :))^2;
    Dloc = [Dloc Dic(:, M(u, 1))]; 
end

if size(Dloc, 2) > 0

Amps = sqrt(abs(Dloc\Y).^2);

for u = 1:size(M, 1)
    eamps = eamps + norm(abs(Amps(u, :)) - abs(Atrue(M(u, 2), :)), 'fro')^2 / size(Amps, 2);
end

end
end