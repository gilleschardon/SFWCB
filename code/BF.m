function [BF] = BF(Y,XX, Pmic, k)


% beamforming

%% Initialisation
r=Y;

Dom = dictionary(Pmic, XX, k);

Domn = Dom ./ (sum(abs(Dom.^2), 1));

BF= sum(abs(Domn'* r).^2, 2) / size(Y, 2);
 end

