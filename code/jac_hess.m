function [jacobian,hessian] = jac_hess(q,Pr,A_r,A_jac,A_hess)
%q(i) = 1xsnapshots
%P_r = nxsnapshots
%A_r = nx1
%A_jac = nx3 dx dy dz
%A_hess = nx6 dxx dyy dzz dyz dxz dxy

%jacobian = 3x1
%hessian = 3x3

new_Pr=Pr - A_r * q;
norm_qi=abs(q.^2);%(q(i)')*q(i);

qPr = q*  (new_Pr)';

dL_dx = 2*real(qPr*A_jac(:, 1));
dL_dy = 2*real(qPr*A_jac(:, 2));
dL_dz = 2*real(qPr*A_jac(:, 3));

ddL_dxx= 2*real(qPr*A_hess(:,1)-(A_jac(:,1)')*A_jac(:,1)*sum(norm_qi));
ddL_dyy= 2*real(qPr*A_hess(:,2)-(A_jac(:,2)')*A_jac(:,2)*sum(norm_qi));
ddL_dzz= 2*real(qPr*A_hess(:,3)-(A_jac(:,3)')*A_jac(:,3)*sum(norm_qi));
ddL_dxy= 2*real(qPr*A_hess(:,6)-(A_jac(:,2)')*A_jac(:,1)*sum(norm_qi));
ddL_dxz= 2*real(qPr*A_hess(:,5)-(A_jac(:,3)')*A_jac(:,1)*sum(norm_qi));
ddL_dyz= 2*real(qPr*A_hess(:,4)-(A_jac(:,3)')*A_jac(:,2)*sum(norm_qi));



jacobian=[dL_dx dL_dy dL_dz].';
hessian = [ddL_dxx ddL_dxy ddL_dxz; ddL_dxy ddL_dyy ddL_dyz; ddL_dxz ddL_dyz ddL_dzz];

end

