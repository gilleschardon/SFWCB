function [X, RE, IM, lambda, err] = sfw_multi_homotopy(Xsensors, k, Data, Xgrid, lambda0, c,errstop, Itermax, LX, UX)

%% SFW homotopy

% Xsensors microphones positions Mx3
% k  wavenumber
% Data acoustical data, MxS (S snapshots)
% Xgrid initialization grid Nx3
% lambda0 first lambda (high)
% c constant for new lambda (<1)
% iterlambdamax number of lambda
% Itermax max iter of each SFW
% LX UX domain


% X estimated positions
% RE IM amplitudes
% lambda
% err residiual


X = cell(1,1);
RE = cell(1, 1);
IM = cell(1, 1);


u = 1;
lambda(u) = lambda0;
tolpos = 0.0001;
tolamp = 0.0001;

[X{u}, RE{u}, IM{u}, nu, err(u)] = sfw_multi(Xsensors, k, Data, Xgrid, lambda(u), tolpos, tolamp, Itermax, LX, UX);
lambda(u+1) = lambda(u) * nu;
iter = 1;
while size(X{u}) < errstop
    iter = iter + 1;
    u = u+1;
[   X{u}, RE{u}, IM{u}, nu, err(u)] = sfw_multi(Xsensors, k, Data, Xgrid, lambda(u), tolpos, tolamp, Itermax, LX, UX, X{u-1}, RE{u-1}, IM{u-1});

lambda(u+1) = lambda(u) / c;
end

lambda = lambda(1:end-1);

end