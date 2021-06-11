%% Comparison of NOMP and SFW, experimental data

% please run FIG_realdata_del_lambda.m first

clear all


load ../data/damasdata150snaps2
close all

%% Parameters 

stepS = 0.1;
stepBF = 0.01;

% domaine
LBx = -2;
UBx = 1;
LBy = -1;
UBy = 0;
LBz = 4;
UBz = 5;
gap = 0.1;
xxS = (LBx+gap:stepS:UBx-gap)';
yyS = (LBy+gap:stepS:UBy-gap)';
zzS = (LBz+gap:stepS:UBz-gap)';

xxCMF = (LBx:stepS:UBx)';
yyCMF = (LBy:stepS:UBy)';
zzCMF = (LBz:stepS:UBz)';
xxBF = (LBx:stepBF:UBx)';
yyBF = (LBy:stepBF:UBy)';
[Xg, Yg, Zg] = meshgrid(xxS, yyS, zzS);
XXS = [Xg(:) Yg(:) Zg(:)];    


[Xg, Yg, Zg] = meshgrid(xxCMF, yyCMF, zzCMF);
XXCMF = [Xg(:) Yg(:) Zg(:)];  


[Xg, Yg, Zg] = meshgrid(xxBF, yyBF, 4.5);
XXBF = [Xg(:) Yg(:) Zg(:)];  






tol1 = 1;
tol2 = 0.1;
tol3 = 0.01;

nbSources = 4;
XX = XXCMF;

Y = data(:, 1:20);

tic 
[xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k);
toc
tic 
[qBF] = BF(Y,XXBF, Pmic, k);
toc
tic
[S_N1,q_N1] = newton_nsnapshot(data(:, 1:10),nbSources,XX,Pmic, tol1, k);
toc
tic
[S_N2,q_N2] = newton_nsnapshot(data(:, 1:10),nbSources,XX,Pmic, tol2, k);
toc
tic
[S_N3,q_N3] = newton_nsnapshot(data(:, 1:10),nbSources,XX,Pmic, tol3, k);
toc


tic
[XSFWm, RE, IM] = sfw_multi_greedy_complex(Pmic, k, data(:, 1:10), XXCMF, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
toc




ASFW = sqrt(RE.^2 + IM.^2);

load lambdapath
Xl = Xn{190};
REl = REn{190};
IMl = IMn{190};
 Al = sqrt(REl.^2 + IMl.^2);
 
Dl = dictionary(Pmic, Xl, k);
norms = sqrt(sum(abs(Dl).^2, 1));

Al = Al ./ norms';

Al2 = Dl\data(:, 1:10);


%%
ratio = 500;
figure

subplot(3, 1, 1)
imagesc(xxBF, yyBF, 10*log10(reshape(qBF, 101, 301)))
axis xy
axis image
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')

subplot(3, 1, 2)
scatter(XSFWm(:, 1), XSFWm(:, 2), mean(abs(ASFW.^2), 2)/ratio, 's', 'linewidth', 2)
hold on
scatter(S_N1(:, 1), S_N1(:, 2), mean(abs(q_N1.^2), 2)/ratio, 'x', 'linewidth', 2)
scatter(S_N2(:, 1), S_N2(:, 2), mean(abs(q_N2.^2), 2)/ratio, '+', 'linewidth', 2)
scatter(S_N3(:, 1), S_N3(:, 2), mean(abs(q_N3.^2), 2)/ratio, '*', 'linewidth', 2)

scatter(xomp(:, 1), xomp(:, 2), mean(abs(q_OMP.^2), 2)/ratio, 'o', 'linewidth', 2)

scatter(Xl(:, 1), Xl(:, 2), mean(abs(Al.^2), 2)/ratio, '^', 'linewidth', 2)

xlim([-2, 1])
ylim([-1, 0])

axis equal
xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])

subplot(3, 1, 3)

scatter(XSFWm(:, 1), XSFWm(:, 3), mean(abs(ASFW.^2), 2)/ratio, 's', 'linewidth', 2)
hold on
scatter(S_N1(:, 1), S_N1(:, 3), mean(abs(q_N1.^2), 2)/ratio, 'x', 'linewidth', 2)
scatter(S_N2(:, 1), S_N2(:, 3), mean(abs(q_N2.^2), 2)/ratio, '+', 'linewidth', 2)
scatter(S_N3(:, 1), S_N3(:, 3), mean(abs(q_N3.^2), 2)/ratio, '*', 'linewidth', 2)

scatter(xomp(:, 1), xomp(:, 3), mean(abs(q_OMP.^2), 2)/ratio, 'o', 'linewidth', 2)
scatter(Xl(:, 1), Xl(:, 3), mean(abs(Al.^2), 2)/ratio, '^', 'linewidth', 2)

xlabel('X (m)')
ylabel('Z (m)')
xlim([-2, 1])
ylim([4.35, 4.72])
legend('SFW greedy', 'NOMP \tau=1', 'NOMP \tau=0.1', 'NOMP \tau=0.01', 'OMP', 'SFW regul.', 'Interpreter', 'TeX')


% 10*log10(mean(abs(ASFW.^2), 2))
% 10*log10(mean(abs(q_N3.^2), 2))
% 10*log10(mean(abs(q_OMP.^2), 2))
% 10*log10(mean(abs(Al.^2), 2))
% 10*log10(mean(abs(Al2.^2), 2))

