%% Homotopy SFW, experimental data
clear all
%% Parameters

stepS = 0.05;
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

xxBF = (LBx:stepS:UBx)';
yyBF = (LBy:stepS:UBy)';
zzBF = (LBz:stepS:UBz)';

[Xg, Yg, Zg] = meshgrid(xxS, yyS, zzS);
XXS = [Xg(:) Yg(:) Zg(:)];    

[Xg, Yg, Zg] = meshgrid(xxBF, yyBF, zzBF);
XBF = [Xg(:) Yg(:) Zg(:)];

%% estimation of source 1 for normalization
load ../data/damasdata150snaps2s1
DD = dictionary(Pmic, XBF, k);

norms = sqrt(sum(abs(DD).^2, 1));
DD = DD ./ norms;

data10 = data(:, 1:10);

S = sum(abs(DD'*data10).^2, 2);

[z, idx] = max(S);

AA = sqrt(mean(abs(DD(:, idx)'*data10).^2));% / norms(idx);

%%

load ../data/damasdata150snaps2
close all


%%

max_source = 20;
alpha = 1.01;% new_lambda = lambda / alpha
lambda0 = 1e5;
itermax = 100; % max iterations of each SFW
nsnaps = 10;

% normalized dictionary

tic
[Xn, REn, IMn, lambdan, errn] = sfw_multi_homotopy_norm(Pmic, k, data(:, 1:nsnaps), XXS, lambda0, alpha, max_source, itermax, [LBx LBy LBz], [UBx UBy UBz]);
toc


% un-normalized
lambda0 = 1e5;
tic
[Xu, REu, IMu, lambdau, erru] = sfw_multi_homotopy(Pmic, k, data(:, 1:nsnaps), XXS, lambda0, alpha, max_source, itermax, [LBx LBy LBz], [UBx UBy UBz]);
toc


save lambdapath
 %%
 
 load lambdapath
 
 figure('Position', [100, 100, 600, 300])

 
 amps = cell(length(REn), 1);
 for u = 1:length(amps)
     ampsn{u} = sqrt(REn{u}.^2+IMn{u}.^2);
 end
 sourcesn = sources_tracking(lambdan,Xn,ampsn);
 
  % Amplitudes in function of lambda

 
hold on
for s = 1:length(sourcesn)
    plot(sourcesn{s}(:, 1), sqrt(mean(sourcesn{s}(:,5:end).^2, 2))/AA, 'linewidth', 2)
end
xlabel('\lambda', 'interpreter', 'TeX')
ylabel('RMS Amplitude (normalized)')
ylim([0,1])



plot(lambdan(190)*[1 1], [0 1000/AA], '--')


 
% Z in function of lambda

 figure('Position', [100, 100, 600, 300])
  subplot(1, 2, 1)

hold on
for s = 1:length(sourcesn)
    plot(sourcesn{s}(:, 1), sourcesn{s}(:, 4), 'linewidth', 2)
end
ylim([3.95, 5.05])
xlabel('\lambda', 'interpreter', 'TeX')
ylabel('Estimated Z coord.')
title('(A) Normalized dictionary')


 subplot(1, 2, 2)

hold on
for s = 1:length(sourcesu)
    plot(sourcesu{s}(:, 1), sourcesu{s}(:, 4), 'linewidth', 2)
end
ylim([3.95, 5.05])
xlabel('\lambda', 'interpreter', 'TeX')
ylabel('Estimated Z coord.')
title('(B) Unnormalized dictionary')

 set(gcf,'Renderer','Painter')

