function [D, Dx, Dy, Dz] = dictionary(PX, PS, k)

% Dictionary of sources and gradient

% PX positions of the sources
% PS positions of the array
% k wavenumber
% source dictionary D
% Dx Dy Dz gradients


dx = PX(:, 1) - PS(:, 1).';
dz = PX(:, 3) - PS(:, 3).';
dy = PX(:, 2) - PS(:, 2).';


d = sqrt(dx.^2 + dz.^2 + dy.^2);

D = exp(- 1i * k * d)./d;

if nargout > 1
    Dd = -1i * k * D - D./d;

    Dx = - Dd .* dx./d;
    Dy = - Dd .* dy./d;
    Dz = - Dd .* dz./d;
end

global appels;
appels = appels + size(D, 2);

end
