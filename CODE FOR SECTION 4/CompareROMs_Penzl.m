% Script to reproduce plots 8(a-d) and Table 3.

%% construct data to use to find estimates of transfer function values
fprintf('Generating data...\n')
load Penzl_disc.mat
n_true = length(A);
rng(9876576)

T = 10000;
t_eval = 0:T;
U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);

red = 14;
freqs = logspace(-5,log10(.99*pi),10*red);
num = length(freqs);
z = exp(1i*freqs);

clear opts
opts.der_order = 1;
opts.num_windows = 40;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
opts.skip_condition_numbers = true;
opts.n = 900; 
%% Calulate Transfer Function Estimates
fprintf('Calculating frequency data\n')
tic
[Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc

%create vectors of input and putput closed under conjugation
Hz_WC = [Hz; conj(Hz)];
nstd_WC = [nstd_Hz(:,1); nstd_Hz(:,1)];
z_WC = [z.';conj(z.')];
[~,idx] = sort(z_WC,'ComparisonMethod','real');
Hz_WC = Hz_WC(idx,:);
nstd_WC = nstd_WC(idx);
z_WC = z_WC(idx,:);

%% Construct 3 ROMS from estimated data
fprintf('Constructing ROMs from learned data...\n')
% change order of interpolation points to interweve for Loewner
zz = [z.';conj(z.')];
HzHz = [Hz(:,1); conj(Hz(:,1))];
zz = [zz(1:2:end);zz(2:2:end)];
HzHz = [HzHz(1:2:end);HzHz(2:2:end)];
[~,idx2] = sort(zz(1:length(z)),'ComparisonMethod','real');
zz = zz([idx2;idx2+length(z)]);
HzHz = HzHz([idx2;idx2+length(z)]);
%Hermite Loewner
epsilon = 1e-8;
[Ap1,Bp1,Cp1,Ep1] = HermiteLoewner(z_WC,Hz_WC(:,1),Hz_WC(:,2),red,epsilon);

%Loewner
[Ap2,Bp2,Cp2,Ep2] = Loewner_sylvester(zz,HzHz(:,1),red,epsilon);

% Vector Fitting
eval_freqs = [freqs.';-freqs.'];
eval_freqs = eval_freqs(idx);

opts2.spy1=0; opts2.spy2=0; opts2.cmplx_ss = 0;
n_iter = 100;
tolVF = 1e-10;
% weight is 4-th root of standard deviation
weights = (1./sqrt(sqrt(nstd_WC)))';
initl_poles = .99*exp(1i*logspace(-3.5,pi,red/2));
initl_poles = [initl_poles, conj(initl_poles)];
initl_poles = sort(initl_poles,'ComparisonMethod','real');
count_VF = 0;
converged = false; diverged = false;
poles = initl_poles;
while ~converged && ~diverged
    count_VF = count_VF + 1;
    [SER,poles,rmserr,~,~]=...
        vectfit3_discrete(Hz_WC(:,1).',eval_freqs.',poles,weights,opts2);
    converged = rmserr < tolVF;
    diverged = count_VF > n_iter;
end
Avf = full(SER.A); Bvf = SER.B; Cvf = SER.C;
%% Construct 3 ROMS from true data
fprintf('Constructing ROMs from true data...\n')

%generate true data
I = eye(length(A));
Hp_func = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_func = @(s) C*((s*I-A)\B);
Hp_true = zeros(num*2,1);
H_interp_true = zeros(num*2,1);
HzHz_true = zeros(num*2,1);
parfor i = 1:num*2
    Hp_true(i) = Hp_func(z_WC(i));
    H_interp_true(i) = H_func(z_WC(i));
    HzHz_true(i) = H_func(zz(i));
end

%Hermite Loewner
[Ap1t,Bp1t,Cp1t,Ep1t] = HermiteLoewner(z_WC,H_interp_true,Hp_true,red,epsilon);
%Loewner
[Ap2t,Bp2t,Cp2t,Ep2t] = Loewner_sylvester(zz,HzHz_true,red,epsilon);
%Vector Fitting
count_VFt = 0;
converged = false; diverged = false;
weightst = ones(1,num*2); %dont weight
polest = initl_poles;
while ~converged && ~diverged
    count_VFt = count_VFt + 1;
    [SERt,polest,rmserrt,~,~]=...
        vectfit3_discrete(H_interp_true.',eval_freqs,polest,weightst,opts2);
    converged = rmserrt < tolVF;
    diverged = count_VFt > n_iter;
end
Avft = full(SERt.A); Bvft = SERt.B; Cvft = SERt.C;
%% Construct system objects and functions to check errors
fprintf('Constructing system objects...\n')
%system objects from approximated data
sysd_Low = ss(Ep2\Ap2,Ep2\Bp2,Cp2,0,1);
sysd_HerLow = ss(Ep1\Ap1,Ep1\Bp1,Cp1,0,1);
sysd_VF = ss(Avf,Bvf,Cvf,0,1);
%system objects from true data
sysd_Lowt = ss(Ep2t\Ap2t,Ep2t\Bp2t,Cp2t,0,1);
sysd_HerLowt = ss(Ep1t\Ap1t,Ep1t\Bp1t,Cp1t,0,1);
sysd_VFt = ss(Avft,Bvft,Cvft,0,1);
% Remove unstable part if present
if max(abs(eig(sysd_Low))) >= 1
    sysd_Low_Full = sysd_Low;
    [sysd_Low, Gus] = stabsep(sysd_Low);
    fprintf('Low was unstable!!!\n')
    fprintf('Stable Low dim: %d\n',length(sysd_Low.A))
end
if max(abs(eig(sysd_Lowt))) >= 1
    [sysd_Lowt, Gus] = stabsep(sysd_Lowt);
    fprintf('Low true was unstable!!!\n')
    fprintf('Stable Low true dim: %d\n',length(sysd_Lowt.A))
end
if max(abs(eig(sysd_HerLow))) >= 1
    sysd_HerLow_Full = sysd_HerLow;
    [sysd_HerLow, Gus] = stabsep(sysd_HerLow);
    fprintf('HerLow was unstable!!!\n')
    fprintf('Stable HerLow dim: %d\n',length(sysd_HerLow.A))
end
if max(abs(eig(sysd_HerLowt))) >= 1
    [sysd_HerLowt, Gus] = stabsep(sysd_HerLowt);
    fprintf('HerLow true was unstable!!!\n')
    fprintf('Stable HerLow true dim: %d\n',length(sysd_HerLowt.A))
end
%% Calculate values and errors on unit circle
fprintf('Calculating errors...\n')
% get points to plot frequency plot from sigma function
sysd = ss(A,B,C,D,1);
[~,freqs_plot] = sigma(sysd);
%Calulate frequency responses for ROMs from estimated data
H_Low = freqresp(sysd_Low,freqs_plot); H_Low = squeeze(H_Low);
H_HerLow = freqresp(sysd_HerLow_Full,freqs_plot); H_HerLow = squeeze(H_HerLow);
H_true = freqresp(sysd,freqs_plot); H_true = squeeze(H_true);
H_VF = freqresp(sysd_VF,freqs_plot); H_VF = squeeze(H_VF);

%Calulate frequency responses for ROMs from true data
H_Lowt = freqresp(sysd_Lowt,freqs_plot); H_Lowt = squeeze(H_Lowt);
H_HerLowt = freqresp(sysd_HerLowt,freqs_plot); H_HerLowt = squeeze(H_HerLowt);
H_VFt = freqresp(sysd_VFt,freqs_plot); H_VFt = squeeze(H_VFt);

err_Low = abs(H_Low - H_true)./abs(H_true);
err_HerLow = abs(H_HerLow - H_true)./abs(H_true);
err_VF = abs(H_VF - H_true)./abs(H_true);
%ROMs from true data
err_Lowt = abs(H_Lowt - H_true)./abs(H_true);
err_HerLowt = abs(H_HerLowt - H_true)./abs(H_true);
err_VFt = abs(H_VFt - H_true)./abs(H_true);
% Error in recovered frequency data
err_interp = norm(Hz_WC(:,1)-H_interp_true)./norm(H_interp_true);
err_interp_der = norm(Hz_WC(:,2)-Hp_true)./norm(Hp_true);
fprintf('2- norm relative error in recovered values:     %e\n',err_interp)
fprintf('2-norm relative error in recovered derivatives: %e\n',err_interp_der)
fprintf('Maximum pointwise realative error in recovered values:      %e\n',max(abs((Hz_WC(:,1)-H_interp_true)./H_interp_true)))
fprintf('Maximum pointwise realative error in recovered derivatives: %e\n',max(abs((Hz_WC(:,2)-Hp_true)./Hp_true)))
%% plot frequency response of ROMs made from estimated data
load ColorMat.mat

figure
loglog(freqs_plot,abs(H_true),'k','LineWidth',2)
hold on
loglog(freqs_plot, abs(H_VF),'-.','Color',ColorMat(1,:),'LineWidth',2)
loglog(freqs_plot, abs(H_HerLow),'--','Color',ColorMat(2,:),'LineWidth',2)
loglog(freqs_plot, abs(H_Low),':','Color',ColorMat(3,:),'LineWidth',2)
legend('$H$','$\hat H_{\texttt{VF}}$','$\hat H_{\texttt{LH}}$','$\hat H_{\texttt{L}}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 14;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
xlim([1e-5,pi])
ylim([10^(-1.5),1e2])
ylabel('$|\hat H_{\texttt X}(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
lgd = legend();
lgd.Location = 'northwest';
text(2e-5,.7e-1,'(a)','FontSize',30)


%% Plot errors of ROMs made from estiamted data
figure
loglog(freqs_plot, err_VF,'LineWidth',2)
hold on
loglog(freqs_plot, err_HerLow,'LineWidth',2)
loglog(freqs_plot, err_Low,'LineWidth',2)
legend('$\hat H_{\texttt{VF}}$','$\hat H_{\texttt{LH}}$','$\hat H_{\texttt{L}}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 14;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
xlim([1e-5,pi])
ylabel('$\epsilon_{rel}$','interpreter','latex','fontsize',20)
lgd = legend();
lgd.Location = 'northwest';
text(2e-5,5e-6,'(b)','FontSize',30)


%% Plot frequency response of ROMs made from true data
figure
loglog(freqs_plot,abs(H_true),'k','LineWidth',2)
hold on
loglog(freqs_plot, abs(H_VFt),'-.','Color',ColorMat(1,:),'LineWidth',2)
loglog(freqs_plot, abs(H_HerLowt),'--','Color',ColorMat(2,:),'LineWidth',2)
loglog(freqs_plot, abs(H_Lowt),':','Color',ColorMat(3,:),'LineWidth',2)
legend('$H$','$\tilde H_{\texttt{VF}}$','$\tilde H_{\texttt{LH}}$','$\tilde H_{\texttt{L}}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 14;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
xlim([1e-5,pi])
ylim([10^(-1.5),1e2])
ylabel('$|\tilde H_{\texttt X}(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
lgd = legend();
lgd.Location = 'northwest';
text(2e-5,.7e-1,'(c)','FontSize',30)


%% Plot errors of ROMs made from true data
figure
loglog(freqs_plot, err_VFt,'LineWidth',2)
hold on
loglog(freqs_plot, err_HerLowt,'LineWidth',2)
loglog(freqs_plot, err_Lowt,'LineWidth',2)
legend('$\tilde H_{\texttt{VF}}$','$\tilde H_{\texttt{LH}}$','$\tilde H_{\texttt{L}}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 14;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
xlim([1e-5,pi])
ylabel('$\epsilon_{rel}$','interpreter','latex','fontsize',20)
lgd = legend();
lgd.Location = 'northwest';
text(2e-5,3e-7,'(d)','FontSize',30)

%% System errors
fprintf('Calculating system errors...\n')
%H2 from approximate data
H2_Sysd = norm(sysd);
H2_Low = norm(sysd-sysd_Low)/H2_Sysd;
H2_HerLow = norm(sysd-sysd_HerLow)/H2_Sysd;
H2_VF = norm(sysd-sysd_VF)/H2_Sysd;

%H2 from true data
H2_Lowt = norm(sysd-sysd_Lowt)/H2_Sysd;
H2_HerLowt = norm(sysd-sysd_HerLowt)/H2_Sysd;
H2_VFt = norm(sysd-sysd_VFt)/H2_Sysd;

%H2 distance of ROM systems
H2_dist_Low = norm(sysd_Low-sysd_Lowt)/norm(sysd_Lowt);
H2_dist_HerLow = norm(sysd_HerLow-sysd_HerLowt)/norm(sysd_HerLowt);
H2_dist_VF = norm(sysd_VF-sysd_VFt)/norm(sysd_VFt);

% H_inf from approximate data
Hinf_Sysd = norm(sysd,'inf');
Hinf_Low = norm(sysd-sysd_Low,'inf')/Hinf_Sysd;
Hinf_HerLow = norm(sysd-sysd_HerLow,'inf')/Hinf_Sysd;
Hinf_VF = norm(sysd-sysd_VF,'inf')/Hinf_Sysd;


% H_inf from true data
Hinf_Lowt = norm(sysd-sysd_Lowt,'inf')/Hinf_Sysd;
Hinf_HerLowt = norm(sysd-sysd_HerLowt,'inf')/Hinf_Sysd;
Hinf_VFt = norm(sysd-sysd_VFt,'inf')/Hinf_Sysd;


%H_inf distance of ROM systems
Hinf_dist_Low = norm(sysd_Low-sysd_Lowt,'inf')/norm(sysd_Lowt,'inf');
Hinf_dist_HerLow = norm(sysd_HerLow-sysd_HerLowt,'inf')/norm(sysd_HerLowt,'inf');
Hinf_dist_VF = norm(sysd_VF-sysd_VFt,'inf')/norm(sysd_VFt,'inf');
%% Output Error norm results

fprintf('------ H2 FROM APX ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    H2_Low, H2_HerLow, H2_VF)
fprintf('------ H2 FROM TRUE ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    H2_Lowt, H2_HerLowt, H2_VFt)
fprintf('------ H2 ROM DISTANCES ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n\n',...
    H2_dist_Low, H2_dist_HerLow, H2_dist_VF)
fprintf('------ Hinf FROM APX ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    Hinf_Low, Hinf_HerLow, Hinf_VF)
fprintf('------ Hinf FROM TRUE ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    Hinf_Lowt, Hinf_HerLowt, Hinf_VFt)
fprintf('------ Hinf ROM DISTANCES ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    Hinf_dist_Low, Hinf_dist_HerLow, Hinf_dist_VF)