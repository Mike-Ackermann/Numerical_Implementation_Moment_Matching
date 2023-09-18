%Constructs ROMs to the fully discretized heat model and compares 


%% construct data to use to find estimates of transfer function values
fprintf('Generating data...\n')
load heat-disc.mat
A = full(E\A); B = full(E\B); C = full(C);
D = 0;

n_true = length(A);

rng(29384234)
T = 1000;
t_eval = 0:T;
U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);
% load Reproduce_HeatModel.mat

red = 10;
num = 50*red;
freqs = logspace(-4,log10(.99*pi),num);

z = exp(1i*freqs);

clear opts
opts.der_order = 1;
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
opts.skip_condition_numbers = true;


%% Calulate Transfer Function Estimates
fprintf('Calculating frequency data\n')
tic
[Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc

%create vectors of input and putput closed under conjugation
Hz_WC = [Hz; conj(Hz)];
z_WC = [z.';conj(z.')];
[~,idx] = sort(z_WC,'ComparisonMethod','real');
Hz_WC = Hz_WC(idx,:);
z_WC = z_WC(idx,:);
%% Construct 3 ROMS from estimated data
fprintf('Constructing ROMs from learned data...\n')
%Hermite Loewner
epsilon = 1e-8;
[Ap1,Bp1,Cp1,Ep1] = HermiteLoewner(z_WC,Hz_WC(:,1),Hz_WC(:,2),red,epsilon);
%Loewner
% change order of interpolation points to interweve
zz = [z.';conj(z.')];
HzHz = [Hz(:,1); conj(Hz(:,1))];
zz = [zz(1:2:end);zz(2:2:end)];
HzHz = [HzHz(1:2:end);HzHz(2:2:end)];
[~,idx2] = sort(zz(1:length(z)),'ComparisonMethod','real');
zz = zz([idx2;idx2+length(z)]);
HzHz = HzHz([idx2;idx2+length(z)]);

[Ap2,Bp2,Cp2,Ep2] = Loewner_sylvester(zz,HzHz,red,epsilon);

% Vector Fitting
eval_freqs = [freqs.';-freqs.'];
eval_freqs = eval_freqs(idx);

opts2.spy1=0; opts2.spy2=0; opts2.cmplx_ss = 0;
n_iter = 100;
tolVF = 1e-5;
weights = ones(1,num*2); %dont weight
%initl_poles = 1*exp((1:red)*2*pi*1i/red); %set initial poles inside upper unit circle
initl_poles = exp(1i*logspace(-1,pi,red/2));
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
Ts = 1;
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
[Ap1t,Bp1t,Cp1t,Ep1t] = HermiteLoewner(z_WC,H_interp_true,Hp_true,red,epsilon);
%Loewner
[Ap2t,Bp2t,Cp2t,Ep2t] = Loewner_sylvester(zz,HzHz_true,red,epsilon);
%Vector Fitting
count_VFt = 0;
converged = false; diverged = false;
polest = initl_poles;
while ~converged && ~diverged
    count_VFt = count_VFt + 1;
    [SERt,polest,rmserrt,~,~]=...
        vectfit3_discrete(H_interp_true.',eval_freqs,polest,weights,opts2);
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
%Print errors in learned transfer function values
freqs_used = freqs/Ts;
err_interp = norm(Hz_WC(:,1)-H_interp_true)./norm(H_interp_true);
err_interp_der = norm(Hz_WC(:,2)-Hp_true)./norm(Hp_true);
fprintf('Error in learned transfer function values:      %e\n',err_interp)
fprintf('Error in learned transfer function derivatives: %e\n',err_interp_der)


%% Calculate trajectory for input to compare
t_eval_test = 0:(1e-3):1;
U_comp = sawtooth(2*pi*10*t_eval_test);
%U_comp = chirp(t_eval,1e-3,t_end,fs/10,'linear');

Y_true = runDTSys(A,B,C,D,U_comp,t_eval_test);
Y_low = runDTSys(Ep2\Ap2,Ep2\Bp2,Cp2,D,U_comp,t_eval_test);
Y_Hlow = runDTSys(Ep1\Ap1,Ep1\Bp1,Cp1,D,U_comp,t_eval_test);
Y_VF = runDTSys(Avf,Bvf,Cvf,D,U_comp,t_eval_test);

%% Plot output responses
%plot bode plot from estimated data
%need to make freqs be in [-log_min_freq, pi)

load ColorMat.mat

figure
plot(t_eval_test,Y_true,'k','LineWidth',2)
hold on
plot(t_eval_test, Y_VF,'-.','Color',ColorMat(1,:),'LineWidth',2)
plot(t_eval_test, Y_low,'--','Color',ColorMat(3,:),'LineWidth',2)
plot(t_eval_test, Y_Hlow,':','Color',ColorMat(2,:),'LineWidth',2)
legend('$Y$','$\hat Y_{\texttt{VF}}$','$\hat Y_{\texttt{L}}$','$\hat Y_{\texttt{LH}}$',...
    'interpreter','latex','Orientation','horizontal')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
xlabel('time (seconds)','interpreter','latex','fontsize',25)
ylabel('$\hat Y_x[t]$','interpreter','latex','fontsize',20)
yticks([-.04,-.02,0,.02,.04])
lgd = legend();
lgd.Location = 'north';

%% plot error trajectory
err_VF = Y_VF-Y_true;
err_Low = Y_low-Y_true;
err_Hlow = Y_Hlow-Y_true;

load ColorMat.mat
%plot errors from estiamted data
figure
plot(t_eval_test, err_VF,'LineWidth',2)
hold on
plot(t_eval_test, err_Low,'Color',ColorMat(3,:),'LineWidth',2)
plot(t_eval_test, err_Hlow,'Color',ColorMat(2,:),'LineWidth',2)
legend('$\hat Y_{\texttt{VF}}$','$\hat Y_{\texttt{L}}$','$\hat Y_{\texttt{LH}}$',...
    'interpreter','latex','Orientation','horizontal')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
xlabel('time (seconds)','interpreter','latex','fontsize',25)
%ylim([-2.5e-9,1.8e-9])
%yticks([-2,-1,0,1]*1e-9)
ylabel('$Y[t]-\hat Y_x[t]$','interpreter','latex','fontsize',20)
lgd = legend();
lgd.Location = 'south';

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

 %H_inf from approximate data
Hinf_Sysd = norm(sysd,'inf');
Hinf_Low = norm(sysd-sysd_Low,'inf')/Hinf_Sysd;
Hinf_HerLow = norm(sysd-sysd_HerLow,'inf')/Hinf_Sysd;
Hinf_VF = norm(sysd-sysd_VF,'inf')/Hinf_Sysd;
 
 
 %H_inf from true data
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