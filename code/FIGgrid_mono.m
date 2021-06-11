%% Influcence of the noise on precision

clear all


load ../data/damasdemo
close all
global appels;

%% Parameters of the simulation

sample_size = 50;

nbSources=3;
snapshots=1;

step = [0.02 0.05 0.1 0.2 0.5];
%step = [0.05 0.1 0.2 0.5];

snr = 20;

nbPlots = length(step);

% position errors (NOMP SFW OMP)
Errors_p=zeros(4,nbPlots, sample_size);
% amplitude errors
Errors_q=zeros(4,nbPlots, sample_size);

nbAppels=zeros(4,nbPlots);

% running time
TN1 = zeros(nbPlots,1); % NOMP
TN2 = zeros(nbPlots,1); % NOMP
TS = zeros(nbPlots, 1); % SFW
TO = zeros(nbPlots, 1); % OMP


% tolerances Newton
tol1 = 1e-7;
tol2 = 1e-9;


% domain
LBx = -1;
UBx = 1;
LBy = -1;
UBy = 1;
LBz = 3;
UBz = 5;



%%

for p = 1:nbPlots
    waitbar(p/nbPlots)

%% Simulation

    waitbar(p/nbPlots);
    % discretisation
xx = (LBx:step(p):UBx)';
yy = (LBy:step(p):UBy)';
zz = (LBz:step(p):UBz)';


[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

XX = [Xg(:) Yg(:) Zg(:)];


    
    for s = 1:sample_size
        
        % source positions
        XS=([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(nbSources, 3) - 0.5)*2;

        % amplitudes
        a = exp(1i*2*pi*rand(3, 1)).*[1;2;4]/10;

        % signal and noise
        Dom = dictionary(Pmic, XS, k);
        Y0 = Dom * a;

        
        noise = randn(size(Y0));
        noise = noise / norm(noise, 'fro') * norm(Y0, 'fro') * 10^(-snr/20);
        Y = Y0 + noise;

        appels = 0;

	% NOMP
        tic
        [S_N1,q_N1] = newton_nsnapshot(Y,nbSources,XX,Pmic, tol1, k);
        [epN1, eaN1 ] = compute_errors(S_N1, XS, q_N1, a);
        TN1(p) = TN1(p) + toc;
        nbAppels(1,p)=nbAppels(1,p)+appels;
        appels = 0;

        tic
        [S_N2,q_N2] = newton_nsnapshot(Y,nbSources,XX,Pmic, tol2, k);
        [epN2, eaN2 ] = compute_errors(S_N2, XS, q_N2, a);
        TN2(p) = TN2(p) + toc;  
        nbAppels(2,p)=nbAppels(2,p)+appels;
        appels = 0;
        tic

	% SFW  
        [XSFW, RE, IM] = sfw_multi_greedy_complex(Pmic, k, Y, XX, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
        TS(p) = TS(p) + toc;
        [epS, eaS ] = compute_errors(XSFW, XS, sqrt(RE.^2+IM.^2), a);
        nbAppels(3,p)=nbAppels(3,p)+appels;

        appels = 0;

	% OMP
        tic
        [xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k);
        [epO, eaO] = compute_errors(xomp, XS, q_OMP, a);
        TO(p) = TO(p) + toc;
        nbAppels(4,p)=nbAppels(4,p)+appels;

        Errors_p(:,p, s)=[epN1, epN2, epS, epO]
        Errors_q(:,p, s)=[eaN1, eaN2, eaS, eaO]
 end
end

Errors_p_m= mean(Errors_p, 3);
Errors_q_m= mean(Errors_q, 3);


save gridsize

%%

load gridsize
close all

figure('Position', [100, 100, 800, 300])

% Plot the points
subplot(1,2,1);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot(step,Errors_p_m(3,:),'-o','LineWidth',2, 'markersize', 10);
plot(step,Errors_p_m(1,:),'--+','LineWidth',2, 'markersize', 10);
plot(step,Errors_p_m(2,:),'-.x','LineWidth',2, 'markersize', 10);
plot(step,Errors_p_m(4,:),':s','LineWidth',2, 'markersize', 10);

ylim([9e-5, 1e4])
xlim([0.018, 0.6])

xlabel("step (m)")
ylabel("MSE position (m^2) ")
legend('SFW', sprintf('NOMP, \\tau = 10e%d ',log10(tol1)),sprintf('NOMP, \\tau = 10e%d ',log10(tol2)), 'OMP','Interpreter','tex')

xticks(step)
yticks([0.0001, 0.01, 1, 100, 10000])

subplot(1,2,2);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot(step,Errors_q_m(3,:),'-o','LineWidth',2, 'markersize', 10);
plot(step,Errors_q_m(1,:),'--+','LineWidth',2, 'markersize', 10);
plot(step,Errors_q_m(2,:),'-.x','LineWidth',2, 'markersize', 10);
plot(step,Errors_q_m(4,:),':s','LineWidth',2, 'markersize', 10);

xticks(step)

xlabel("step (m)")
ylabel("MSE amplitude (Pa^2)")
yticks([0.001, 0.01, 0.1,1, 10, 100])
figure('Position', [100, 100, 800, 300])

subplot(1,2,1);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot(step, TS/sample_size,'-o','LineWidth',2, 'markersize', 10);
plot(step,TN1/sample_size,'--+','LineWidth',2, 'markersize', 10);
plot(step,TN2/sample_size,'-.x','LineWidth',2, 'markersize', 10);
plot(step,TO/sample_size,':s','LineWidth',2, 'markersize', 10);

xlabel("step (m)")
ylabel("Time (s)")
xticks(step)
ylim([1e-3, 1e2])
xlim([0.018, 0.6])
yticks([0.001, 0.01, 0.1,1, 10,100])


subplot(1,2,2);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
plot(step, nbAppels(3, :)/sample_size,'-o','LineWidth',2, 'markersize', 10);
plot(step,nbAppels(1, :)/sample_size,'--+','LineWidth',2, 'markersize', 10);
plot(step,nbAppels(2, :)/sample_size,'-.x','LineWidth',2, 'markersize', 10);
plot(step,nbAppels(4, :)/sample_size,':s','LineWidth',2, 'markersize', 10);
legend('SFW', sprintf('NOMP, \\tau = 10e%d ',log10(tol1)),sprintf('NOMP, \\tau = 10e%d ',log10(tol2)), 'OMP','Interpreter','tex')
xticks(step)

xlabel("step (m)")
ylabel("#points")

yticks([100, 1000,10000, 100000,1000000])

xlim([0.018, 0.6])
