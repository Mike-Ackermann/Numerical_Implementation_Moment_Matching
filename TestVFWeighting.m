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
%freqs = [freqs 1.5 .99*pi];
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

% Calculate true frequency information

I = eye(length(A));
Hp_func = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_func = @(s) C*((s*I-A)\B);
H_interp_true = zeros(num*2,1);
rel_errs = zeros(num*2,1);
H_perturbed = zeros(num*2,1);
parfor i = 1:num*2
    H_interp_true(i) = H_func(z_WC(i));
    rel_errs(i) = abs(Hz_WC(i)-H_interp_true(i))/abs(H_interp_true(i));
end
H_perturbed(1:2:end) = H_interp_true(1:2:end) + H_interp_true(1:2:end).*rel_errs(1:2:end).*randn(num,1);
H_perturbed(2:2:end) = conj(H_perturbed(1:2:end));

%% Make VF Models
% Vector Fitting
eval_freqs = [freqs.';-freqs.']/Ts;
eval_freqs = eval_freqs(idx);

opts2.spy1=0; opts2.spy2=0; opts2.cmplx_ss = 0;
n_iter = 100;
tolVF = 1e-10;
%weights = (1./nstd_WC)';
weights = ones(1,num*2);
%initl_poles = 1*exp((1:nmax)*2*pi*1i/nmax); %set initial poles inside upper unit circle
initl_poles = exp(1i*logspace(-3,pi,red/2));
initl_poles = [initl_poles, conj(initl_poles)];
initl_poles = sort(initl_poles,'ComparisonMethod','real');
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

% Vector fitting from true data
count_VFt = 0;
converged = false; diverged = false;
polest = initl_poles;
while ~converged && ~diverged
    count_VFt = count_VFt + 1;
    [SERt,polest,rmserrt,~,~]=...
        vectfit3_discrete(H_interp_true.',eval_freqs,polest,weights,Ts,opts2);
    converged = rmserrt < tolVF;
    diverged = count_VFt > n_iter;
end
Avft = full(SERt.A); Bvft = SERt.B; Cvft = SERt.C;

% Vector fitting from perturbed data (no weighting)
count_VFp = 0;
converged = false; diverged = false;
polesp = initl_poles;
while ~converged && ~diverged
    count_VFp = count_VFp + 1;
    [SERp,polesp,rmserrp,~,~]=...
        vectfit3_discrete(H_perturbed.',eval_freqs,polesp,weights,Ts,opts2);
    converged = rmserrp < tolVF;
    diverged = count_VFp > n_iter;
end
Avfp = full(SERp.A); Bvfp = SERp.B; Cvfp = SERp.C;

% Vector fitting from perturbed data (with weighting)
count_VFp_w = 0;
converged = false; diverged = false;
weights2 = (1./rel_errs)';
polesp_w = initl_poles;
while ~converged && ~diverged
    count_VFp_w = count_VFp_w + 1;
    [SERp_w,polesp_w,rmserrp_w,~,~]=...
        vectfit3_discrete(H_perturbed.',eval_freqs,polesp_w,weights2,Ts,opts2);
    converged = rmserrp_w < tolVF;
    diverged = count_VFp_w > n_iter;
end
Avfp_w = full(SERp_w.A); Bvfp_w = SERp_w.B; Cvfp_w = SERp_w.C;

%Vector fitting from Discovered data with weighting
count_VF_w = 0;
converged = false; diverged = false;
poles_w = initl_poles;
while ~converged && ~diverged
    count_VF_w = count_VF_w + 1;
    [SER_w,poles_w,rmserr_w,~,~]=...
        vectfit3_discrete(Hz_WC(:,1).',eval_freqs,poles_w,weights2,Ts,opts2);
    converged = rmserr_w < tolVF;
    diverged = count_VF_w > n_iter;
end
Avf_w = full(SER_w.A); Bvf_w = SER_w.B; Cvf_w = SER_w.C;

%% Construct system objects and find H2 errors
fprintf('Finding H2 errors...\n')
sysd = ss(A,B,C,D,Ts);
sysd_VF = ss(Avf,Bvf,Cvf,0,Ts);
sysd_VFt = ss(Avft,Bvft,Cvft,0,Ts);
sysd_VFp = ss(Avfp,Bvfp,Cvfp,0,Ts);
sysd_VFp_w = ss(Avfp_w,Bvfp_w,Cvfp_w,0,Ts);
sysd_VF_w = ss(Avf_w,Bvf_w,Cvf_w,0,Ts);

H2_Sysd = norm(sysd);
H2_VF = norm(sysd-sysd_VF)/H2_Sysd;
H2_VFt = norm(sysd-sysd_VFt)/H2_Sysd;
H2_VFp = norm(sysd-sysd_VFp)/H2_Sysd;
H2_VFp_w = norm(sysd-sysd_VFp_w)/H2_Sysd;
H2_VF_w = norm(sysd-sysd_VF_w)/H2_Sysd;

fprintf('------ H2 ERRORS ------\n')
fprintf('Discovered: %e, Weighted Discovered: %e,\n Perturbed: %e, Weighted Perturbed: %e,\n True: %e\n',...
    H2_VF, H2_VF_w, H2_VFp,H2_VFp_w,H2_VFt)
