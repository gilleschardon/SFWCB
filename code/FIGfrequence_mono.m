%% Influcence of the noise on precision

load ../data/damasdemo
close all

%% Parameters of the simulation

sample_size = 50;

nbSources=3;
snapshots=1;
nbPlots = 9;

step = 0.05;

f = logspace(2.5, 4.5, nbPlots);
k = 2 * pi /340 * f;

snr = 20;

% position errors (NOMP SFW OMP)
Errors_p=zeros(3,nbPlots, sample_size);
% amplitude errors
Errors_q=zeros(3,nbPlots, sample_size);

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
        a = exp(1i*2*pi*rand(3, 1)).*[1;2;4]/10;

        % signal and noise
        Dom = dictionary(Pmic, XS, k(p));
        Y0 = Dom * a;

        noise = randn(size(Y0));
        noise = noise / norm(noise, 'fro') * norm(Y0, 'fro') * 10^(-snr/20);
        Y = Y0 + noise;

        % NOMP
        tic
        [S_N,q_N] = newton_nsnapshot(Y,nbSources,XX,Pmic, tol, k(p));        
        [epN, eaN ] = compute_errors(S_N, XS, q_N, a); 
        TN(p) = TN(p) + toc;
           
        % SFW
        tic
        [XSFW, RE, IM] = sfw_multi_greedy_complex(Pmic, k(p), Y, XX, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
        [epS, eaS] = compute_errors(XSFW, XS, sqrt(RE.^2+IM.^2), a);
        TS(p) = TS(p) + toc;

        % OMP
        tic
        [xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k(p));
        [epO, eaO] =compute_errors(xomp, XS, q_OMP, a);

        TO(p) = TO(p) + toc;

        Errors_p(:,p, s) = [epN, epS, epO];
        Errors_q(:,p, s)= [eaN, eaS, eaO];
    end
end

Errors_p_m= mean(Errors_p, 3);
Errors_q_m= mean(Errors_q, 3);

Eplast = sort(Errors_p, 3,'ascend');
Eqlast = sort(Errors_q, 3,'ascend');

Errors_p_m= mean(Eplast(:, :, 1:end-1), 3);
Errors_q_m= mean(Eqlast(:, :, 1:end-1), 3);

save freq

%% Plot the points

load freq

close all

figure('Position', [100, 100, 800, 300])

subplot(1,2,1);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

loglog(f,Errors_p_m(2,:),'-o','LineWidth',2, 'markersize', 10);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(f,Errors_p_m(1,:),'--x','LineWidth',2, 'markersize', 10);
loglog(f,Errors_p_m(3,:),':s','LineWidth',2, 'markersize', 10);
xticks([100,1000,10000])

xlabel("f (Hz)")
ylabel("MSE position (m^2) ")
yticks([0.00001, 0.0001,0.001, 0.01,0.1, 1, 10])

subplot(1,2,2);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

loglog(f,Errors_q_m(2,:),'-o','LineWidth',2, 'markersize', 10);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(f,Errors_q_m(1,:),'--x','LineWidth',2, 'markersize', 10);
loglog(f,Errors_q_m(3,:),':s','LineWidth',2, 'markersize', 10);
legend('SFWm','NOMP', 'OMP','Interpreter','tex')
xticks([100,1000,10000])
yticks([0.01, 0.1,1, 10,100]/100)


xlabel("f (Hz)")
ylabel("MSE amplitude (Pa^2)")

