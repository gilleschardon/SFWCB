function [Xs, RE, IM, nu, err] = sfw_multi_greedy_complex(Xm, k, Data, Xgrid, tolpos, tolamp, Niter, LX, UX, Xs, amps, phases)


options_nu = optimoptions(@fmincon,'Display', 'off', 'Algorithm','sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false, 'OptimalityTolerance', 1e-12, 'MaxFunctionEvaluations', 1e7);
options_amps = optimoptions(@fmincon,'Display', 'off', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false,  'OptimalityTolerance', 1e-12, 'MaxFunctionEvaluations', 1e7, 'StepTolerance', 1e-14);
options_all= optimoptions(@fmincon,'Display', 'off', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient',false, 'CheckGradient', false,  'OptimalityTolerance', 1e-12, 'MaxFunctionEvaluations', 1e7, 'StepTolerance', 1e-14);


Dgrid = dictionary(Xm, Xgrid, k);

Dgrid = Dgrid ./ sqrt(sum(abs(Dgrid).^2, 1));

Nsnap = size(Data, 2);

if nargin < 11
    Xs = zeros(0, 3);
    RE = zeros(0, Nsnap);
    IM = zeros(0, Nsnap);
end

C = zeros(size(Data));
tolnu = 1.01;

lambda = 0;

for u = 1:Niter
    
    % residual
    Dloc = dictionary(Xm, Xs, k);
    A = RE + 1i * IM;
    C = Dloc * A;
    R = Data - C;
    
    fprintf("Iteration %u\t %u sources\n", u, size(Xs, 1))
    
    %% Adding a source
    
    [xnew, nu] = maximize_nu(Xm, k, lambda, Xgrid, Dgrid, R, LX, UX, options_nu);
    

    % adding the source to the set
    Xs = [Xs ; xnew];

    %% Optimization of the amplitudes
    % local dictionary
    Dloc =  dictionary(Xm, Xs, k);
    
    dd = dictionary(Xm, xnew, k);
    aa = dd'*R / norm(dd);
    
    RE = [RE ; real(aa)];
    IM = [IM ; imag(aa)];


    [RE, IM] = optimize_amplitudes(Dloc, Data, lambda, RE, IM, options_amps);

    %% Joint optimization of the amplitudes and positions
    [Xs, RE, IM] = optimize_all(Xm, k, lambda, Data, Xs, RE, IM, LX, UX, options_all);
    
    %% Cleaning (optional)
    %[Xs, RE, IM] = clean(Xs, RE, IM, tolamp, tolpos);
    % wrapping the phases
    %phases = mod(phases, 2*pi);


end
%% Max iter. reached
fprintf("Max iter, stopping\n")
    Dloc = dictionary(Xm, Xs, k);
    A = RE + 1i * IM;
    C = Dloc * A;
    R = Data - C;
    err = norm(R, 'fro');


end


%% maximize nu(x)
% Xm sensors
% Xgrid grid
% Dgrid dictionary on the grid
% R residual
% options options
%
% Xnu argmax
% nu max

function [Xnu, nu] = maximize_nu(Xm, k, lambda, Xgrid, Dgrid, R, LX, UX, options)

nugrid = sum(abs(Dgrid'*R).^2, 2);
[~, idx] = max(nugrid);
xnewgrid = Xgrid(idx, :);


nuf = @(X) nux_multi(Xm, k, R, X);
[Xnu, numin] = fmincon(nuf, xnewgrid, [], [], [], [], LX, UX, [], options);

nu = sqrt(-(numin))/lambda;
end

function [RE, IM] = optimize_amplitudes(Dloc, Data, lambda, RE, IM, options)

Nsnap = size(Data, 2);

[amps] = fmincon(@(x) obj_amplitudes_multi(Dloc, Data, lambda, x), [RE; IM], [], [], [], [], [] ,[], [], options);
Ns = size(amps, 1)/2;

RE = amps(1:Ns, :);
IM = amps(Ns+1:end, :);



end

function [Xs, RE, IM] = optimize_all(Xm, k, lambda, Data, Xs, RE, IM, LX, UX, options)



xopt = [Xs RE IM];

% bounds
[Ns, Nsnap] = size(RE);


lbounds = [ones(Ns,1)*LX(1)  ones(Ns, 1)*LX(2) ones(Ns, 1)*LX(3) -Inf(Ns, Nsnap) -Inf(Ns, Nsnap)];
ubounds = [ones(Ns,1)*UX(1)  ones(Ns, 1)*UX(2) ones(Ns, 1)*UX(3) Inf(Ns, Nsnap) Inf(Ns, Nsnap)]; % no upper bounds on amplitudes

%obj_amplitudes_positions_multi(Xm, k, Data, lambda, xopt)

ZZ = fmincon(@(x) obj_amplitudes_positions_multi(Xm, k, Data, lambda, x), xopt, [], [], [], [], lbounds, ubounds, [], options);
%obj_amplitudes_positions_multi(Xm, k, Data, lambda, ZZ)

% extaction of the amplitudes and positions
Xs = ZZ(:, 1:3);
RE = ZZ(:, 4:3+Nsnap);
IM = ZZ(:, 4+Nsnap:end);
end

function [Xs, RE, IM] = clean(Xs, RE, IM, tolamp, tolpos)


idxr = mean(RE.^2 + IM.^2, 2) < tolamp;
IM(idxr, :) = [];
RE(idxr, :) = [];
Xs(idxr, :) = [];



dists = (Xs(:,1) - Xs(:,1)').^2 + (Xs(:,2) - Xs(:,2)').^2 + (Xs(:,3) - Xs(:,3)').^2 + 1000*eye(size(Xs, 1));
[mm, i1] = min(dists);
[m, i2] = min(mm);

i1 = i1(i2);
while m < tolpos    
    xxx = Xs(i1, :);
    Xs([i1 i2], :) = [];
    
    Xs = [Xs ; xxx];
    
    REz = RE(i1, :) + RE(i2, :);
    IMz = IM(i1, :) + IM(i2, :);

    RE([i1 i2], :) = [];
    IM([i1 i2], :) = [];
    
    RE = [RE ; REz];
    IM = [IM ; IMz];
    
    dists = (Xs(:,1) - Xs(:,1)').^2 + (Xs(:,2) - Xs(:,2)').^2 + (Xs(:,3) - Xs(:,3)').^2 + 1000*eye(size(Xs, 1));
    [mm, i1] = min(dists);
    [m, i2] = min(mm);
    i1 = i1(i2);
end

end

function [J, Jgrad] = obj_amplitudes_positions_multi(Xm, k, Data, lambda, x)

Ns = size(x, 1);

Xs = x(:, 1:3);

Nsnap = (size(x, 2)-3)/2;

xx = x(:, 4:end);

[D, Dx, Dy, Dz] = dictionary(Xm, Xs, k);

RE = xx(:, 1:Nsnap);

IM = xx(:, Nsnap+1:end);

A = RE + 1i * IM;
C = D * A;

Delta = C - Data;

J = 1/2 * norm(Delta, 'fro').^2 + lambda * sum(sqrt(sum(RE(:, 1:Nsnap).^2 + IM(:, 1:Nsnap).^2, 2)));

if nargout > 1
    JgradRE = zeros(Ns, Nsnap);
    JgradIM = zeros(Ns, Nsnap);
    Jgradpos = zeros(Ns, 3);

    for u = 1:Ns
        if norm(RE(u, :) + 1i*IM(u, :) , 2) > 0
            JgradRE(u, :) = real( sum(conj((D(:, u))).*Delta, 1)) + lambda * RE(u, :) / norm(RE(u, :) + 1i*IM(u, :), 2);
            JgradIM(u, :) = real( sum(conj(1i*(D(:, u))).*Delta, 1)) + lambda * IM(u, :) / norm(RE(u, :) + 1i*IM(u, :), 2);

        else
             JgradRE(u, :) = real( sum(conj((D(:, u))).*Delta, 1));
            JgradIM(u, :) = real( sum(conj(1i*(D(:, u))).*Delta, 1));
           
        end
       
        Jgradpos(u, 1) = real(trace(Delta' * Dx(:, u)*A(u, :)));
        Jgradpos(u, 2) = real(trace(Delta' * Dy(:, u)*A(u, :)));
        Jgradpos(u, 3) = real(trace(Delta' * Dz(:, u)*A(u, :)));

    end
    Jgrad = [Jgradpos JgradRE JgradIM];
    
end
end

function [J,Jgrad] = obj_amplitudes_multi(D, Data, lambda, x)

Nsnap = size(Data, 2);

Ns = size(x, 1)/2;

A = x(1:Ns, 1:Nsnap) + 1i * x(Ns+1:2*Ns, 1:Nsnap);
C = D * A;

Delta = C - Data;

J = 1/2 * norm(Delta, 'fro')^2 + lambda * sum(sqrt(sum(x.^2, 2)));


if nargout > 1
    Jgrad = zeros(size(x));
    for u = 1:Ns
        if norm(x(u, :) + 1i * x(u+Ns, :), 2) > 0
            Jgrad(u, :) = real( sum(conj(D(:, u)).*Delta, 1)) + lambda * x(u, :) / norm(x(u, :) + 1i * x(u+Ns, :), 2);
            Jgrad(u+Ns, :) = real( sum(conj(1i*D(:, u)).*Delta, 1)) + lambda * x(u+Ns, :) / norm(x(u, :) + 1i * x(u+Ns, :), 2);
        else
            Jgrad(u, :) = real( sum(conj((D(:, u))).*Delta, 1));
            Jgrad(u, :) = real( sum(conj((1i*D(:, u))).*Delta, 1));

        end
    end
end
end

function [nu] = nux_multi(Xm, k, Delta, Xnu)


    [d, gx, gy, gz] = dictionary(Xm, Xnu, k);
    d = d / norm(d);

    nu = - sum(abs(d'*Delta).^2);

end

