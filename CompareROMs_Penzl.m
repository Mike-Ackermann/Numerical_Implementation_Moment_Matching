%Plots errors in only Loewner, Hermite Loewner, approximations using
%Estimated TF values.  This is for a chosen, fixed r.

%% construct data to use to find estimates of transfer function values
fprintf('Generating data...\n')
load Penzl_disc.mat

fs = 1e4;%1e3?
Ts = 1/fs;
t_end = 1;

t_eval = 0:Ts:t_end;
T = length(t_eval);
%U = FiltNoise(fs,t_end);
%Y = runDTSys(A,B,C,D,U,t_eval);
load BestInput_Penzl.mat

red = 14;
freqs = logspace(-4,log10(.99*pi),10*red);
%freqs = logspace(-2.2924,log10(.99*pi),10*red);
%freqs = logspace(-1,log10(.99*pi),10*red);
num = length(freqs);
r = 1; % radius of points
z = r*exp(1i*freqs);

clear opts
opts.tol = 10^(-1);
opts.der_order = 1;
opts.num_est = 10;
opts.n = 300; %set n to 300 to follow section on estimating moments

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

%Hz_WC = Hz_WC_cheat;


%% Construct 3 ROMS from estimated data
fprintf('Constructing ROMs from learned data...\n')
nmax = red;
%Hermite Loewner
epsilon = 1e-8;
[Ap1,Bp1,Cp1,Ep1] = HermiteLoewner(z_WC,Hz_WC(:,1),Hz_WC(:,2),nmax,epsilon);
%Loewner
% change order of interpolation points to interweve
zz = [z.';conj(z.')];
HzHz = [Hz(:,1); conj(Hz(:,1))];
zz = [zz(1:2:end);zz(2:2:end)];
HzHz = [HzHz(1:2:end);HzHz(2:2:end)];
[~,idx2] = sort(zz(1:length(z)),'ComparisonMethod','real');
zz = zz([idx2;idx2+length(z)]);
HzHz = HzHz([idx2;idx2+length(z)]);

[Ap2,Bp2,Cp2,Ep2] = Loewner(zz,HzHz,nmax,epsilon);

% Vector Fitting
eval_freqs = [freqs.';-freqs.']/Ts;
eval_freqs = eval_freqs(idx);

opts2.spy1=0; opts2.spy2=0; opts2.cmplx_ss = 0;
n_iter = 100;
tolVF = 1e-10;
weights = (1./sqrt(sqrt(nstd_WC)))';
%weights = ones(1,num*2);
%initl_poles = 1*exp((1:nmax)*2*pi*1i/nmax); %set initial poles inside upper unit circle
initl_poles = .99*exp(1i*logspace(-3.5,pi,nmax/2));
initl_poles = [initl_poles, conj(initl_poles)];
initl_poles = sort(initl_poles,'ComparisonMethod','real');
%%%%%%%%%%%%%%%
%initl_poles = .99*exp(1i*logspace(-3.5,-1,floor(nmax/4)));
%initl_poles = [initl_poles, conj(initl_poles)];
%initl_poles = [initl_poles, linspace(.9,1,1+nmax/2)];
%initl_poles = sort(initl_poles,'ComparisonMethod','real');
%%%%%%%%%%%%%%%
count_VF = 0;
converged = false; diverged = false;
poles = initl_poles;
while ~converged && ~diverged
    count_VF = count_VF + 1;
    [SER,poles,rmserr,~,~]=...
        vectfit3_discrete(Hz_WC(:,1).',eval_freqs.',poles,weights,Ts,opts2);
    converged = rmserr < tolVF;
    diverged = count_VF > n_iter;
end
Avf = full(SER.A); Bvf = SER.B; Cvf = SER.C;
%% Construct 3 ROMS from true data
fprintf('Constructing ROMs from true data...\n')

%generate true data
sysd = ss(A,B,C,D,Ts);
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
[Ap1t,Bp1t,Cp1t,Ep1t] = HermiteLoewner(z_WC,H_interp_true,Hp_true,nmax,epsilon);
%Loewner
[Ap2t,Bp2t,Cp2t,Ep2t] = Loewner(zz,HzHz_true,nmax,epsilon);
%Vector Fitting
count_VFt = 0;
converged = false; diverged = false;
weightst = ones(1,num*2); %dont weight
polest = initl_poles;
while ~converged && ~diverged
    count_VFt = count_VFt + 1;
    [SERt,polest,rmserrt,~,~]=...
        vectfit3_discrete(H_interp_true.',eval_freqs,polest,weightst,Ts,opts2);
    converged = rmserrt < tolVF;
    diverged = count_VFt > n_iter;
end
Avft = full(SERt.A); Bvft = SERt.B; Cvft = SERt.C;
%% Construct system objects and functions to check errors
fprintf('Constructing system objects...\n')
%system objects from approximated data
sysd_Low = ss(Ep2\Ap2,Ep2\Bp2,Cp2,0,Ts);
sysd_HerLow = ss(Ep1\Ap1,Ep1\Bp1,Cp1,0,Ts);
sysd_VF = ss(Avf,Bvf,Cvf,0,Ts);
%system objects from true data
sysd_Lowt = ss(Ep2t\Ap2t,Ep2t\Bp2t,Cp2t,0,Ts);
sysd_HerLowt = ss(Ep1t\Ap1t,Ep1t\Bp1t,Cp1t,0,Ts);
sysd_VFt = ss(Avft,Bvft,Cvft,0,Ts);
% Remove unstable part if present
if max(abs(eig(sysd_Low))) >= 1
    sysd_Low_Full = sysd_Low;
    [sysd_Low, Gus] = stabsep(sysd_Low);
    fprintf('Low was unstable!!!\n')
    fprintf('Stable Low dim: %d\n',length(sysd_Low.A))
    %fprintf('Hinf unstable part: %e\n',norm(Gus,'inf'));
end
if max(abs(eig(sysd_Lowt))) >= 1
    [sysd_Lowt, Gus] = stabsep(sysd_Lowt);
    fprintf('Low true was unstable!!!\n')
    fprintf('Stable Low true dim: %d\n',length(sysd_Lowt.A))
    %fprintf('Hinf unstable part: %e\n',norm(Gus,'inf'));
end
if max(abs(eig(sysd_HerLow))) >= 1
    [sysd_HerLow, Gus] = stabsep(sysd_HerLow);
    fprintf('HerLow was unstable!!!\n')
    fprintf('Stable HerLow dim: %d\n',length(sysd_HerLow.A))
    %fprintf('Hinf unstable part: %e\n',norm(Gus,'inf'));
end
if max(abs(eig(sysd_HerLowt))) >= 1
    [sysd_HerLowt, Gus] = stabsep(sysd_HerLowt);
    fprintf('HerLow true was unstable!!!\n')
    fprintf('Stable HerLow true dim: %d\n',length(sysd_HerLowt.A))
    %fprintf('Hinf unstable part: %e\n',norm(Gus,'inf'));
end
%% Calculate values and errors on unit circle
fprintf('Calculating errors...\n')
%construct a random-log spaced distribution between log_min_freq and pi
%true frequencies interpolated at
%random frequencies for testing
freqs_used = freqs/Ts;
% get points to plot frequency plot from sigma function
[~,freqs_plot] = sigma(sysd);
%Calulate frequency responses at random points for ROMs from estimated data
H_Low = freqresp(sysd_Low,freqs_plot); H_Low = squeeze(H_Low);
H_HerLow = freqresp(sysd_HerLow,freqs_plot); H_HerLow = squeeze(H_HerLow);
H_true = freqresp(sysd,freqs_plot); H_true = squeeze(H_true);
%H_true_used = freqresp(sysd,freqs_used); H_true_used = squeeze(H_true_used);
H_true_used = H_interp_true;
H_VF = freqresp(sysd_VF,freqs_plot); H_VF = squeeze(H_VF);

%Calulate frequency responses at random points for ROMs from true data
H_Lowt = freqresp(sysd_Lowt,freqs_plot); H_Lowt = squeeze(H_Lowt);
H_HerLowt = freqresp(sysd_HerLowt,freqs_plot); H_HerLowt = squeeze(H_HerLowt);
H_VFt = freqresp(sysd_VFt,freqs_plot); H_VFt = squeeze(H_VFt);
%calculate errors.  Cant calculate error of TF estimates because these
%frequencies don't were not the interpolation frequencies.
err_Low = abs(H_Low - H_true)./abs(H_true);
err_HerLow = abs(H_HerLow - H_true)./abs(H_true);
err_VF = abs(H_VF - H_true)./abs(H_true);
%ROMs from true data
err_Lowt = abs(H_Lowt - H_true)./abs(H_true);
err_HerLowt = abs(H_HerLowt - H_true)./abs(H_true);
err_VFt = abs(H_VFt - H_true)./abs(H_true);
% Error in learned frequency data
err_interp = norm(Hz_WC(:,1)-H_true_used)./norm(H_true_used);
err_interp_der = norm(Hz_WC(:,2)-Hp_true)./norm(Hp_true);
fprintf('Error in learned transfer function values:      %e\n',err_interp)
fprintf('Error in learned transfer function derivatives: %e\n',err_interp_der)
%% plot bode plot from estimated data
%plot bode plot from estimated data
%need to make freqs be in [-log_min_freq, pi)

load ColorMat.mat

freqs_plt = freqs_plot*Ts;
figure
loglog(freqs_plt,abs(H_true),'k','LineWidth',2)
hold on
loglog(freqs_plt, abs(H_VF),'-.','Color',ColorMat(1,:),'LineWidth',2)
loglog(freqs_plt, abs(H_HerLow),'--','Color',ColorMat(2,:),'LineWidth',2)
loglog(freqs_plt, abs(H_Low),':','Color',ColorMat(3,:),'LineWidth',2)
legend('$H$','$\hat H_r^{VF}$','$\hat H_r^{HL}$','$\hat H_r^{L}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
xticks([1e-4,1e-3,1e-2,1e-1,1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
%set limits of plot
xlim([1e-4,1])
%labels
ylabel('$|\hat H_r^x(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
% xlabel('$r$','interpreter','latex','fontsize',30)
lgd = legend();
lgd.Location = 'northwest';

%% plot errors from estiamted data
%plot errors from estiamted data
figure
%loglog(freqs_plt_interp,err_interp,'LineWidth',2)
loglog(freqs_plt, err_VF,'LineWidth',2)
hold on
loglog(freqs_plt, err_HerLow,'LineWidth',2)
loglog(freqs_plt, err_Low,'LineWidth',2)
legend('$\hat H_r^{VF}$','$\hat H_r^{HL}$','$\hat H_r^{L}$',...
    'interpreter','latex')
%title('Relative error at random points on unit circle estiamted data')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
xticks([freqs_plt(1),1e-1,pi])
xticks([1e-4,1e-3,1e-2,1e-1,1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
%set limits of plot
xlabel('$\omega$','interpreter','latex','fontsize',25)
xlim([1e-4,pi])
%labels
ylabel('$|H(e^{\mathbf i \omega})-\hat H_r^x(e^{\mathbf i \omega})|/|H(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
% xlabel('$r$','interpreter','latex','fontsize',30)
lgd = legend();
lgd.Location = 'northwest';

%% plot bode plot from true data
%plot bode plot from estimated data
%need to make freqs be in [-log_min_freq, pi)
freqs_plt = freqs_plot*Ts;
figure
loglog(freqs_plt,abs(H_true),'k','LineWidth',2)
hold on
loglog(freqs_plt, abs(H_VFt),'-.','Color',ColorMat(1,:),'LineWidth',2)
loglog(freqs_plt, abs(H_HerLowt),'--','Color',ColorMat(2,:),'LineWidth',2)
loglog(freqs_plt, abs(H_Lowt),':','Color',ColorMat(3,:),'LineWidth',2)
legend('$H$','$\tilde H_r^{VF}$','$\tilde H_r^{HL}$','$\tilde H_r^{L}$',...
    'interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
xticks([1e-4,1e-3,1e-2,1e-1,1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
%set limits of plot
xlim([1e-4,pi])
%labels
ylabel('$|\tilde H_r^x(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
% xlabel('$r$','interpreter','latex','fontsize',30)
lgd = legend();
lgd.Location = 'northwest';

%% plot errors from true data
%plot errors from estiamted data
figure
%loglog(freqs_plt_interp,err_interp,'LineWidth',2)
loglog(freqs_plt, err_VFt,'LineWidth',2)
hold on
loglog(freqs_plt, err_HerLowt,'LineWidth',2)
loglog(freqs_plt, err_Lowt,'LineWidth',2)
legend('$\tilde H_r^{VF}$','$\tilde H_r^{HL}$','$\tilde H_r^{L}$',...
    'interpreter','latex')
%title('Relative error at random points on unit circle estiamted data')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
xticks([1e-4,1e-3,1e-2,1e-1,1])
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
%set limits of plot
xlim([1e-4,pi])
%labels
ylabel('$|H(e^{\mathbf i \omega})-\tilde H_r^x(e^{\mathbf i \omega})|/|H(e^{\mathbf i \omega})|$','interpreter','latex','fontsize',20)
% xlabel('$r$','interpreter','latex','fontsize',30)
lgd = legend();
lgd.Location = 'northwest';
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

% %H_inf from approximate data
% Hinf_Sysd = norm(sysd,'inf');
% Hinf_Low = norm(sysd-sysd_Low,'inf')/Hinf_Sysd;
% Hinf_HerLow = norm(sysd-sysd_HerLow,'inf')/Hinf_Sysd;
% Hinf_VF = norm(sysd-sysd_VF,'inf')/Hinf_Sysd;
% 
% 
% %H_inf from true data
% Hinf_Lowt = norm(sysd-sysd_Lowt,'inf')/Hinf_Sysd;
% Hinf_HerLowt = norm(sysd-sysd_HerLowt,'inf')/Hinf_Sysd;
% Hinf_VFt = norm(sysd-sysd_VFt,'inf')/Hinf_Sysd;

%H_inf distance of ROM systems
%Hinf_dist_HerLow = norm(sysd_HerLow-sysd_HerLowt,'inf')/norm(sysd_HerLowt,'inf');

%% Output Error norm results

fprintf('------ H2 FROM APX ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    H2_Low, H2_HerLow, H2_VF)
fprintf('------ H2 FROM TRUE ERRORS ------\n')
fprintf('Low: %e, HerLow: %e, VF: %e\n',...
    H2_Lowt, H2_HerLowt, H2_VFt)
% fprintf('------ Hinf FROM APX ERRORS ------\n')
% fprintf('Low: %e, HerLow: %e, VF: %e\n',...
%     Hinf_Low, Hinf_HerLow, Hinf_VF)
% fprintf('------ Hinf FROM TRUE ERRORS ------\n')
% fprintf('Low: %e, HerLow: %e, VF: %e\n',...
%     Hinf_Lowt, Hinf_HerLowt, Hinf_VFt)

%Matrix of H2 errors for converting to LaTex format
%exclude loewner because not in thesis
H2_err_mat = [H2_Low, H2_HerLow, H2_VF;...
              H2_Lowt, H2_HerLowt, H2_VFt;...
              H2_dist_Low H2_dist_HerLow, H2_dist_VF];

%Hinf_err_mat = [Hinf_HerLow, Hinf_VF, Hinf_AAA;
%                Hinf_HerLowt, Hinf_VFt, Hinf_AAAt;...
%                Hinf_dist_HerLow, Hinf_dist_VF, Hinf_dist_AAA];

%matrix2latex(H2_err_mat, 'H2_err_Penzl.txt')
%matrix2latex(Hinf_err_mat, 'Hinf_err_RandEx1_correct.txt')