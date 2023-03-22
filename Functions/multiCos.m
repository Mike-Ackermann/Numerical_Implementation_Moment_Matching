function U = multiCos(t_eval,fs,freqs)
% generates a multicosine signal
% Inputs:
% t_eval - times to evaluate signal at
% fs - sampling frequency
% freqs - frequencies included in output signal

% Outputs:
% U - multi cosine signal

U = zeros(size(t_eval));
for k = 1:length(freqs)
    U = U + cos(t_eval*freqs(k)*fs);
end