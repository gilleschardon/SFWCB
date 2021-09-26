%% DEMO

% application of greedy and penalized SFW to experimental data

clear all


load data5
close all

%% Parameters 

step = 0.1;

% domaine
LBx = -2;
UBx = 1;
LBy = -1;
UBy = 0;
LBz = 4;
UBz = 5;
gap = 0.1;
xx = (LBx+gap:step:UBx-gap)';
yy = (LBy+gap:step:UBy-gap)';
zz = (LBz+gap:step:UBz-gap)';


[Xg, Yg, Zg] = meshgrid(xx, yy, zz);
XX = [Xg(:) Yg(:) Zg(:)];    




Y = data(:, 1:10);

%% greedy

lambda = 0; % greedy, no lambda
nbSources = 4; % greedy, stopping at 4 iterations
tolpos = 0; % no merging of sources at identical positions
tolamp = 0; % no removal of sources with 0 amplitude
% Pmic positions of the microphones
% k wavenumber
% Y snapshots

tic
[XSFWm, RE, IM] = sfw_multi_norm(Pmic, k, Y, XX, lambda, tolpos, tolamp, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
toc

% amplitudes in RE and IM are given for the normalized dictionary, we
% compute the amplitudes in the original dictionary
D = dictionary(Pmic, XSFWm, k);
norms = sqrt(sum(abs(D).^2, 1));

ASFW = sqrt(RE.^2 + IM.^2);
ASFW = ASFW ./ norms.';

%% penalized

lambda = 650;
tolpos = 0.001;
tolamps = 0.001;
nbSources = 1000; % we stop when nu <=1, not when a fixed number of sources is found

tic
[Xl, REl, IMl] = sfw_multi_norm(Pmic, k, Y, XX, lambda, tolpos, tolamp, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
toc

% amplitudes in RE and IM are given for the normalized dictionary, we
% compute the amplitudes in the original dictionary
Dl = dictionary(Pmic, XSFWm, k);
norms = sqrt(sum(abs(Dl).^2, 1));

ASFWl = sqrt(REl.^2 + IMl.^2);
ASFWl = ASFWl ./ norms.';
%%
MS= 200;
MSgt = 40;
lw = 2;

ratio = 1000;
figure

subplot(2, 1, 1)
scatter(XSFWm(:, 1), XSFWm(:, 2), MS, '+', 'linewidth', lw)
hold on

scatter(Xl(:, 1), Xl(:, 2), MS, 'x', 'linewidth', lw)

axis equal

xlim([-2, 1])
ylim([-1, 0])

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])

legend('SFW gr.', 'SFW p.')


subplot(2, 1, 2)
scatter(XSFWm(:, 1), XSFWm(:, 3), MS, '+', 'linewidth', lw)
hold on

scatter(Xl(:, 1), Xl(:, 3), MS, 'x', 'linewidth', lw)


axis equal

xlim([-2, 1])
ylim([4, 5])

xlabel('X (m)')
ylabel('Z (m)')
xlim([-2, 1])

