%% MSE, number of snapshots

clear all

% load some data (mainly microphone positions)
load ../data/damasdemo
close all

%% Parameters of the simulation

nbSources=3;
snapshots=1:10;

step = 0.05;

snr = 10;

sample_size = 50;

nbPlots = length(snapshots);

% position errors (NOMP SFW OMP)
Errors_p=zeros(4,nbPlots, sample_size);
% amplitude err5ors
Errors_q=zeros(4,nbPlots, sample_size);

% running time
TN = zeros(nbPlots,1); % NOMP
TS = zeros(nbPlots, 1); % SFW
TO = zeros(nbPlots, 1); % OMP

tol = 1e-9; % NOMP stopping criterion

% domain
LBx = -1;
UBx = 1;
LBy = -1;
UBy = 1;
LBz = 3;
UBz = 5;

xx = (LBx:step:UBx)';
yy = (LBy:step:UBy)';
zz = (LBz:step:UBz)';

[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

XX = [Xg(:) Yg(:) Zg(:)];


stepM = 0.05;
xxM = (LBx:stepM:UBx)';
yyM = (LBy:stepM:UBy)';
zzM = (LBz:stepM:UBz)';
[Xgm, Ygm] = meshgrid(xxM, yyM);
[Xgm3, Ygm3, Zgm3] = meshgrid(xxM, yyM, zzM);

XXM = [Xgm3(:) Ygm3(:) Zgm3(:)];

% for the waitbar
ntotal = nbPlots * sample_size;
nn = 0;


for p = 1:nbPlots

%% Simulation
    
    for s = 1:sample_size
        waitbar(nn/ntotal);
        nn = nn + 1;
        
        
        % source positions
        XS=([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(nbSources, 3) - 0.5)*2;

        % amplitudes
        a = (randn(nbSources,snapshots(p)) + 1i * randn(nbSources,snapshots(p))) .* [1;2;4]/10/sqrt(2);

        % signals
        Dom = dictionary(Pmic, XS, k);
        Y0 = Dom * a;
        % noise
        noise = randn(size(Y0));
        noise = noise / norm(noise, 'fro') * norm(Y0, 'fro') * 10^(-snr/20);
        
        % noisy signal
        Y = Y0 + noise;

        % NOMP
        tic
        [S_N,q_N] = newton_nsnapshot(Y,nbSources,XX,Pmic, tol, k);
        [epN, eaN ] = compute_errors(S_N, XS, q_N, a);
        TN(p) = TN(p) + toc;
           
        % SFW
        tic
        [XSFW, RE, IM] = sfw_multi_norm(Pmic, k, Y, XX, 0, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
        TS(p) = TS(p) + toc;
        
        Dl = dictionary(Pmic, XSFW, k);
        norms = sqrt(sum(abs(Dl).^2, 1));

        Asfw = sqrt(RE.^2+IM.^2) ./ norms';

        
        [epS, eaS ] = compute_errors(XSFW, XS, Asfw, a);

            
        
        % OMP
        tic
        [xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k);
        [epO, eaO] = compute_errors(xomp, XS, q_OMP, a);
        TO(p) = TO(p) + toc;
        
        
epM = inf;
eaM = inf;
    if snapshots(p) > 3
        Xmusic = MUSIC_local(Y, nbSources,XXM,Pmic, k);
        Dmusic = dictionary(Pmic, Xmusic, k);
        [epM, eaM] = compute_errors_music(Xmusic, XS, Y, a, Dmusic);
    end


    
          Errors_p(:,p, s) = [epN, epS, epO, epM];
        Errors_q(:,p, s)= [eaN, eaS, eaO, eaM];
        
    end 
end

% MSE
Errors_p_m= mean(Errors_p, 3);
Errors_q_m= mean(Errors_q, 3);

save snap
%% Plot the points
load snap

figure('Position', [100, 100, 800, 300])
msize = 15;
subplot(1,2,1);
hold on
set(gca, 'YScale', 'log')
plot(snapshots,Errors_p_m(2,:),'-o','LineWidth',2, 'markersize', msize);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);
plot(snapshots,Errors_p_m(1,:),'--x','LineWidth',2, 'markersize', msize);
plot(snapshots,Errors_p_m(3,:),':s','LineWidth',2, 'markersize', msize);
plot(snapshots,Errors_p_m(4,:),':^','LineWidth',2, 'markersize', msize);


xlabel("Snapshots")
ylabel("MSE position (m^2) ")
yticks([0.0001, 0.001,0.01, 0.1,1])

subplot(1,2,2);
hold on
set(gca, 'YScale', 'log')
plot(snapshots,Errors_q_m(2,:),'-o','LineWidth',2, 'markersize', msize);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);
plot(snapshots,Errors_q_m(1,:),'--x','LineWidth',2, 'markersize', msize);
plot(snapshots,Errors_q_m(3,:),':s','LineWidth',2, 'markersize', msize);
plot(snapshots,Errors_q_m(4,:),':^','LineWidth',2, 'markersize', msize);

xlabel("Snapshots")
ylabel("MSE amplitudes (Pa^2)")
legend('SFW gr.','NOMP', 'OMP', 'MUSIC', 'Interpreter','tex')

