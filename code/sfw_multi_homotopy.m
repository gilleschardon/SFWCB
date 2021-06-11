function [X, RE, IM, lambda, err] = sfw_multi_homotopy(Xsensors, k, Data, Xgrid, lambda0, c,errstop, Itermax, LX, UX)


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