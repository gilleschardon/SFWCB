%% MSE, frequency

load ../data/damasdemo
close all

%% Parameters of the simulation

sample_size = 50;

nbSources=3;
snapshots=1;
nbPlots = 11;

step = 0.05;

f = logspace(2.5, 5, nbPlots);

k = 2 * pi /340 * f;

snr = 20;

% position errors (NOMP SFW OMP)
Errors_X=zeros(3,nbPlots, sample_size);
Errors_Y=zeros(3,nbPlots, sample_size);
Errors_Z=zeros(3,nbPlots, sample_size);

% amplitude errors
Errors_q=zeros(3,nbPlots, sample_size);

% running time
TN = zeros(nbPlots,1); % NOMP
TS = zeros(nbPlots, 1); % SFW
TO = zeros(nbPlots, 1); % OMP

tol = 1e-9; % NOMP stopping criterion


% domain
LBx = -1;
UBx =  1;
LBy = -1;
UBy =  1;
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
        [epNX, epNY, epNZ, eaN ] = compute_errors_coherenceXYZ(S_N, XS, q_N, a, Pmic, k(p)); 
        TN(p) = TN(p) + toc;
           
        % SFW
        tic
        [XSFW, RE, IM] = sfw_multi_norm(Pmic, k(p), Y, XX, 0, 0, 0, nbSources, [LBx LBy LBz]-0.1, [UBx UBy UBz]+0.1);
        TS(p) = TS(p) + toc;

        Dest = dictionary(Pmic, XSFW, k(p));
        norms = sqrt(sum(abs(Dest).^2, 1))';
        
        RE = RE ./ norms;
        IM = IM ./ norms;

               
        [epSX, epSY, epSZ, eaS] = compute_errors_coherenceXYZ(XSFW, XS, RE + 1i*IM, a, Pmic, k(p));

       
        
        % OMP
        tic
        [xomp, q_OMP] = OMP(Y,nbSources,XX, Pmic, k(p));
        [epOX, epOY, epOZ, eaO] =compute_errors_coherenceXYZ(xomp, XS, q_OMP, a, Pmic, k(p));

        TO(p) = TO(p) + toc;

        Errors_X(:,p, s) = [epNX, epSX, epOX];
        Errors_Y(:,p, s) = [epNY, epSY, epOY];
        Errors_Z(:,p, s) = [epNZ, epSZ, epOZ];

        Errors_q(:,p, s)= [eaN, eaS, eaO];
        
        
        
    end
end

Errors_X_m= mean(Errors_X, 3);
Errors_Y_m= mean(Errors_Y, 3);
Errors_Z_m= mean(Errors_Z, 3);

Errors_p_m = Errors_X_m + Errors_Y_m + Errors_Z_m;

Errors_q_m= mean(Errors_q, 3);
save freq

%% Plot the points

load freq

close all
msize = 15;
figure('Position', [100, 100, 800, 300])

subplot(1,2,1);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

loglog(f,Errors_X_m(2,:) + Errors_Y_m(2,:)+ Errors_Z_m(2,:) ,'-o','LineWidth',2, 'markersize', msize);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(f,Errors_X_m(1,:)+ Errors_Y_m(1,:)+ Errors_Z_m(1,:),'--x','LineWidth',2, 'markersize', msize);
loglog(f,Errors_X_m(3,:)+ Errors_X_m(3,:)+ Errors_Z_m(3,:) ,':s','LineWidth',2, 'markersize', msize);


xticks([100,1000,10000])

xlabel("f (Hz)")
ylabel("MSE position (m^2) ")
yticks([0.00001, 0.0001,0.001, 0.01,0.1, 1, 10])

subplot(1,2,2);
hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

loglog(f,Errors_q_m(2,:),'-o','LineWidth',2, 'markersize', msize);
set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1);

loglog(f,Errors_q_m(1,:),'--x','LineWidth',2, 'markersize', msize);
loglog(f,Errors_q_m(3,:),':s','LineWidth',2, 'markersize', msize);

legend('SFW gr.','NOMP', 'OMP', 'Interpreter','tex')
xticks([100,1000,10000])
yticks([0.01, 0.1,1, 10,100]/100)


xlabel("f (Hz)")
ylabel("MSE amplitude (Pa^2)")

