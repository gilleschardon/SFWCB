%% Comparison of NOMP and SFW, experimental data

% please run FIG_realdata_det_lambda.m first

clear all


load data5
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

xxBF = (LBx:stepBF:UBx)';
yyBF = (LBy:stepBF:UBy)';
[Xg, Yg, Zg] = meshgrid(xxS, yyS, zzS);
XX = [Xg(:) Yg(:) Zg(:)];    

[Xg, Yg, Zg] = meshgrid(xxBF, yyBF, 4.5);
XXBF = [Xg(:) Yg(:) Zg(:)];  






tol1 = 1;
tol2 = 0.1;
tol3 = 0.01;

nbSources = 4;

Y = data(:, 1:10);

tic 
[xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k);
toc
tic 
[qBF] = BF(Y,XXBF, Pmic, k);
toc

tic
[S_N3,q_N3] = newton_nsnapshot(Y,nbSources,XX,Pmic, tol3, k);
toc


tic
[XSFWm, RE, IM] = sfw_multi_norm(Pmic, k, Y, XX, 0, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
toc


Dl2 = dictionary(Pmic, XSFWm, k);
norms = sqrt(sum(abs(Dl2).^2, 1));

ASFW = sqrt(RE.^2 + IM.^2);
ASFW = ASFW ./ norms.';
load lambdapath
Xl = Xn{190};
REl = REn{190};
IMl = IMn{190};
 Al = sqrt(REl.^2 + IMl.^2);
 
Dl = dictionary(Pmic, Xl, k);
norms = sqrt(sum(abs(Dl).^2, 1));

Al = Al ./ norms';

Al2 = Dl\data(:, 1:10);



D = dictionary(Pmic, XX, k);

norms = sqrt(sum(abs(D).^2, 1));
%%
MS=  200;
MSgt = 40;
lw = 2;
load loc1
ratio = 1000;
figure

C1= [0, 0.4470, 0.7410];
C2 = [0.8500, 0.3250, 0.0980];
C3 =           	[0.9290, 0.6940, 0.1250];
C4 = [0.4940, 0.1840, 0.5560];
C5 = [0.4660, 0.6740, 0.1880];
C6 = 	[0.3010, 0.7450, 0.9330];
C7 = [0.5 0.5 0.5];

subplot(1, 4, 1)
imagesc(xxBF, yyBF, 10*log10(reshape(qBF, 101, 301)))
axis xy
axis image
colormap(hot)
colorbar
xlabel('X (m)')
ylabel('Y (m)')




subplot(2, 3, 2)
scatter(XSFWm(:, 1), XSFWm(:, 2), MS, C1, '+', 'linewidth', lw)
hold on

scatter(Xl(:, 1), Xl(:, 2), MS, C2, 'x', 'linewidth', lw)
scatter(XSFWloc(:, 1),XSFWloc(:, 2), MSgt, C7, 'filled')

axis equal

xlim([-2, 1])
ylim([-1, 0])

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])

legend('SFW gr.', 'SFW p.')


subplot(2, 3, 5)
scatter(XSFWm(:, 1), XSFWm(:, 3), MS, C1, '+', 'linewidth', lw)
hold on

scatter(Xl(:, 1), Xl(:, 3), MS, C2, 'x', 'linewidth', lw)


scatter(XSFWloc(:, 1),XSFWloc(:, 3), MSgt, C7, 'filled')

axis equal

xlim([-2, 1])
ylim([4, 5])

xlabel('X (m)')
ylabel('Z (m)')
xlim([-2, 1])










subplot(2, 3, 3)

scatter(S_N3(:, 1), S_N3(:, 2), MS, C3, '+', 'linewidth', lw)
hold on
scatter(xomp(:, 1), xomp(:, 2), MS, C6, 'x', 'linewidth', lw)

scatter(XSFWloc(:, 1),XSFWloc(:, 2), MSgt, C7, 'filled')

legend('NOMP', 'OMP')

axis equal

xlim([-2, 1])
ylim([-1, 0])

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])


subplot(2, 3, 6)

scatter(S_N3(:, 1), S_N3(:, 3),MS, C3, '+', 'linewidth', lw)
hold on

scatter(xomp(:, 1), xomp(:, 3), MS, C6, 'x', 'linewidth', lw)


scatter(XSFWloc(:, 1),XSFWloc(:, 3), MSgt, C7, 'filled')


axis equal

xlim([-2, 1])
ylim([4, 5])

xlabel('X (m)')
ylabel('Z (m)')


