function Xest = MUSIC_local(Y, nbSources,XX,Pmic, k)

C = Y * Y';

[Dom] = dictionary(Pmic, XX, k);

Domnorm = Dom ./ sqrt( sum(abs(Dom).^2, 1));

[V, D] = eigs(C, nbSources, 'largestabs');

ps = sum( abs(V'*Domnorm).^2, 1);

size(ps)
ps = reshape(ps, 51, 151, 51);



mmax = movmax(movmax(movmax(ps, 3, 1), 3, 2), 3, 3);

ps(ps ~= mmax) = 0;

ps(ps~= 0)
XX(ps~=0, :)

idx = find(ps > 0.3);

XX(idx, :);
Xest = [];


if length(idx) > 0
for u = 1:length(idx)
    
    
    xopt = fminunc(@(x) obj(x, Pmic, k, V), XX(idx(u), :));

    Xest = [Xest ; xopt];
    

end
end

uu = 1;
while uu <= size(Xest, 1)
    cur = Xest(uu, :);
    vv = uu+1;
    while vv <= size(Xest, 1)

        if norm(cur - Xest(vv, :)) < 1e-4
            Xest(vv, :) = [];
            vv = vv - 1;
        end
        vv = vv + 1;
    end
    uu = uu + 1;
end



end

function f = obj(x, Pmic, k, V)

    d = dictionary(Pmic, x, k);
    d = d / norm(d);

    f = - sum( abs(V'*d).^2, 1);
    
end