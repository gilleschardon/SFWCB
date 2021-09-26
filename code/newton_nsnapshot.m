function [S,q, err] = newton_nsnapshot(Y, nbSources,XX,Pmic,  tol, k)

%% newtonized OMP
% Y data
% nbSources iterations
% XX grid
% Pmic microphone positions
% tol tolerance for local optimizations
% k wavenumber

% S positions
% q amplitudes


if nargin < 5
    tol = 1e-6;
end

DD = 3;

[nbMic, ~] = size(Pmic);
[nbNoeuds, ~]=size(XX);

Pr=Y;

snapshots = size(Y, 2);
S=zeros(nbSources,3);
A=zeros(nbMic,nbSources);
q=zeros(nbSources,snapshots);

[Dom] = dictionary(Pmic, XX, k);

Domnorm = Dom ./ sqrt( sum(abs(Dom).^2, 1));



for s = 1:nbSources
    %%Grid based initial estimation
%     rho=0;
%     n_star=1;
%     for n=1:nbNoeuds
%         new_rho=norm(Domnorm(:,n)'*Pr)^2;
%         nugrid(n) = new_rho;
%         if rho<new_rho
%             n_star = n;
%             rho=new_rho;
%         end
%     end

    rhorho = sum(abs(Domnorm'*Pr).^2, 2);
    [~, n_star] = max(rhorho);

    %XX(n_star,:)
    
    S(s,:)=XX(n_star,:);

    q(s,:)=(Dom(:,n_star)'*Pr)/(norm(Dom(:,n_star),2)^2);
    %%Local Newton Optimization

    [dloc,gradloc,hessloc] = a_grad_hessian(Pmic, XX(n_star, :), k);

    [jacobian,hessian]=jac_hess(q(s,:),Pr,dloc,gradloc,hessloc);


    step = (hessian(1:DD, 1:DD)\jacobian(1:DD)).';



    S(s,1:DD)=S(s,1:DD)-step;

    A(:,s)=dictionary(Pmic,S(s, :), k);

    q(s,:)=(A(:,s)'*Pr)/(norm(A(:,s),2)^2);

    %S

    new_Pr=Y;
    for n=1:s
        new_Pr=new_Pr-A(:,n)*q(n,:);
    end
    %%Global cyclic feedback optimization
    l=0;
    while abs(norm(Pr,'fro')^2-norm(new_Pr,'fro')^2)>tol && l<100000
        for i=1:s
           [~,grad_a,hess_a]=a_grad_hessian(Pmic, S(i,:), k);
           interm_Pr=new_Pr+A(:,i)*q(i,:);
           interm_q=(A(:,i)'*interm_Pr)/(norm(A(:,i),2)^2);
           [jacobian,hessian]=jac_hess(interm_q,interm_Pr,A(:,i),grad_a,hess_a);
           step = (hessian(1:DD, 1:DD)\jacobian(1:DD)).';

           S(i,1:DD)=S(i,1:DD)-step;


           A(:,i)=dictionary(Pmic,S(i, :), k);
           q(i,:)=(A(:,i)'*interm_Pr)/(norm(A(:,i),2)^2);


           Pr=new_Pr;
           new_Pr=interm_Pr-A(:,i)*q(i,:);
        end
        l=l+1;

    end

    %%Orthogonal solution
    q(1:s,:)=A(:,1:s) \ Y;
    Pr=Y-A*q;
   % S
   err =  norm(Pr, 'fro');

end

end
