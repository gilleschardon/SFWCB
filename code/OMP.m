function [S, q_omp, BF] = OMP(Y,nbSources,XX, Pmic, k)


%OMP this fonciton compute the inidce of the sources using the OMP method
%   Detailed explanation goes here

%% Initialisation
r=Y;
S=zeros(nbSources,3);
DS=[];
Dom = dictionary(Pmic, XX, k);

Domn = Dom ./ sqrt(sum(abs(Dom.^2), 1));


for i = 1:nbSources
%% Computation of the correlations of the residual with each of the atoms of the dictionary,and identification of the maximum
    rho = sum(abs(Domn'* r).^2, 2);
    if i == 1
        BF = rho;
    end
    
%rho=sum(((Dom./sum(Dom.*conj(Dom),2))'*r).*(((Dom'*r)').'),2);
    [~,l_star]=max(rho);
    
%% Mise � jour du r�sidut
    DS=[DS,Dom(:,l_star)];
    P=DS*((DS'*DS)\(DS'));
    new_r=(eye(size(P))-P)*Y;
    r=new_r; 
%% Ajout de l'indice de la source k � l'ensemble des indices
    S(i,:)=XX(l_star,:);

end

q_omp = DS\Y;

end

