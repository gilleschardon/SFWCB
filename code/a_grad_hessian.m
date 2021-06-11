function [a, a1, a2] = a_grad_hessian(PA, PX, k)

% like dictionary.m, but one source, and computes the hessian as well, for
% NOMP


global appels;
appels = appels + 1;
dx = PX(:,1).' - PA(:, 1);
dy = PX(:,2).' - PA(:, 2);
dz = PX(:,3).' - PA(:, 3);

d = sqrt(dx.^2 + dz.^2 + dy.^2);

g = exp(- 1i * k * d)./d;
g1 = (-1i * k - 1./d) .* g;
g2 = (-k^2 + 2i*k./d + 2./d.^2) .* g;

a = g;
a1 = [dx,dy,dz].*g1./d;

Cx = (dx./d).^2;
Cy = (dy./d).^2;
Cz = (dz./d).^2;

Cxy = dx.*dy./d.^2;
Cyz = dy.*dz./d.^2;
Cxz = dx.*dz./d.^2;



%a2_1 = [(1 - Cx)./d .*  g1 + Cx .*g2,  (1 - Cy)./d .*  g1 + Cy .*g2 , (1 - Cz)./d .*  g1 + Cz .*g2];
a2_1 = [(1 - Cx) (1 - Cy) (1 - Cz)]./d.*  g1 + [Cx Cy Cz] .*  g2;


a2_2 = [Cyz, Cxz, Cxy] .* (-g1./d +g2);

% x2 y2 z2 yz xz xy
a2 = [a2_1 , a2_2];

end
