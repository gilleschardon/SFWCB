function [BF] = BF(Y,XX, Pmic, k)


%OMP this fonciton compute the inidce of the sources using the OMP method
%   Detailed explanation goes here

%% Initialisation
r=Y;

Dom = dictionary(Pmic, XX, k);

Domn = Dom ./ (sum(abs(Dom.^2), 1));

BF= sum(abs(Domn'* r).^2, 2) / size(Y, 2);
 end

