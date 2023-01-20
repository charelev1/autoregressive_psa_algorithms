clc; clear; close all;

N=1024; % size of input signal
m=5; % order of the AR model 

% input signal parameters
k1 = N/3; % freq k (0 to N/2)
k2 = N/7; % freq k (0 to N/2)
A1 = 10; % Amplitude of sine wave
A2 = 10; % Amplitude of sine wave
N1 = 5; % Amplitude of added noise

% add input signal
w1 = 2*pi*k1/N; % normalized freq (0 to 1)
w2 = 2*pi*k2/N; % normalized freq (0 to 1)
x1 = A1*sin(w1*(0: N-1)');
x2 = A2*sin(w2*(0: N-1)');
noise = N1*randn([N, 1]); % gaussian noise
x = x1+x2+noise; % final input signal

% AR models
[aburg, sburg] = burg(x, m);
[aburg_fast, sburg_fast] = burg_fast(x, m);
[aburg_fft, sburg_fft] = burg_fast_fft(x, m);

% AR models (matlab)
[amat_burg, smat_burg] = arburg(x, m);
[amat_yule, smat_yule] = aryule(x, m);

%% PSD estimation methods
% FFT PSD method
Pfft = abs((fft(x))); 
% AR PSD methods
wout = 2*pi*(0:N-1)/N;
Pburg = ar_psd(sburg, aburg, m, wout);
Pmat_burg = ar_psd(smat_burg, amat_burg, m, wout);
Pmat_yule = ar_psd(smat_yule, amat_yule, m, wout);


%% Plot section
% Time domain plots
figure();
title('Time domain signals')
plot(x);
legend('True input');

% Freq domain plots
figure();
    total_figs=2;
    subplot(total_figs,1,1)
    plot(Pfft);
    title('FFT-based PSD');
    xlim([1, N/2]);
    subplot(total_figs,1,2)
    plot(Pburg);
    title('Burg-based PSD');
    xlim([1, N/2]);

% Freq domain plots (matlab)
figure();
    total_figs=2;
    subplot(total_figs,1,1)
    plot(Pmat_burg);
    title('Burg-based PSD');
    xlim([1, N/2]);
    subplot(total_figs,1,2)
    plot(Pmat_yule);
    title('Yule-based PSD');
    xlim([1, N/2]);
