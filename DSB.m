clear all;
close all;
[y, fs] = audioread('eric.wav');
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(Y));

% Plot the spectrum
plot(f, abs(Y)/length(Y));
title('Spectrum of m');
figure();

% Apply filter
bw = 4000;
filt = ones(size(Y));
filt(f > bw|f<-bw) = 0;
y_filter = Y .* filt;
plot(f, abs(y_filter)/length(y_filter));
title('filtered signal spectrum');
% Inverse transform
y_filtered_time = ifft(ifftshift(y_filter));

% Ensure data type and scaling
y_filtered_time = real(y_filtered_time); % Ensure real values
y_filtered_time = double(y_filtered_time); % Convert to double if necessary

%% Normalize if values are outside the range [-1, 1]
max_val = max(abs(y_filtered_time));
if max_val > 1
    y_filtered_time = y_filtered_time / max_val;
end
fc = 100000;
U = 0.5;
Am = max(y_filtered_time);
Ac = Am/U; %modulationindex = Am/Ac
new_fs = 5 * fc;

resampled_signal = resample(y_filtered_time, new_fs, fs);
t1 = linspace(0, length(resampled_signal) / new_fs, length(resampled_signal));
t1 = t1';
carrier = Ac .* cos(2 * pi * fc * t1);
DSB_SC = resampled_signal .* carrier;
DSB_TC = (1 + U * resampled_signal/Am ) .* carrier;
%% Plot DSB-SC spectrum
figure();
subplot(1,2,1)
DSB_SC_spectrum = fftshift(fft(DSB_SC));
f_DSB_SC =linspace(- new_fs/2,  new_fs/2, length(DSB_SC));
plot(f_DSB_SC, abs(DSB_SC_spectrum)/length(DSB_SC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-SC Modulated Signal Spectrum');
% Plot DSB-TC spectrum
subplot(1,2,2)
DSB_TC_spectrum = fftshift(fft(DSB_TC));
f_DSB_TC=linspace(-new_fs/2,new_fs/2,length(DSB_TC));
plot(f_DSB_TC, abs(DSB_TC_spectrum)/length(DSB_TC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-TC Modulated Signal Spectrum');

%% envelop for DSB-SC
envelope_DSB_SC = abs(hilbert(DSB_SC));
figure;
plot(t1, DSB_SC);
hold on;
plot(t1,-envelope_DSB_SC,'r-',t1, envelope_DSB_SC,'-r','Linewidth',1.5); 
hold off;
title('DSB_sc with envellope');
ylim([-2 2]);
xlim([2 2.5]);

%envelop for DSB-TC
envelope_DSB_TC = abs(hilbert(DSB_TC));
figure;
plot(t1, DSB_TC);
hold on;
plot(t1, -envelope_DSB_TC ,'r-', t1, envelope_DSB_TC, '-r', 'Linewidth', 1.5); 
hold off;
title('DSB_Tc with envellope');
ylim([-2 2]);
xlim([2 2.5]);


demod_DSB_SC=envelope_DSB_SC.*cos(2 * pi * fc * t1);
%% plot demodulation signal
figure;
subplot(2, 1, 1);
plot(t1, envelope_DSB_TC);
xlabel('Time');
ylabel('Amplitude');
title('Received Signal - DSB-TC');
subplot(2, 1, 2);
plot(t1, demod_DSB_SC);
xlabel('Time');
ylabel('Amplitude');
title('Received Signal - DSB-SC');
% resample the envelope DSB_SC to hear it
recived_sig_DSB_SC = resample(envelope_DSB_SC, fc, fs);
% resample the envelope DSB_TC to hear it
recived_sig_DSB_TC = resample(envelope_DSB_TC, fc, fs);
%% sound the three msgs
sounds={y_filtered_time,envelope_DSB_SC,envelope_DSB_TC};
for i = 1:length(sounds)
    sound(sounds{i}, fs);
    pause(10); 
end



%% Coherent detection
% when noise is added to DSB_SC with SNR = 0, 10, and 30
snr_values = [0, 10, 30];
for snr_dB = snr_values
    % Add noise to DSB_SC
    noisy_DSB_SC = awgn(DSB_SC, snr_dB);
    demodulated_noisy = noisy_DSB_SC .* cos(2*pi*fc*t1);
    f = new_fs / 2 * linspace(-1,1,length(DSB_SC));
    % Coherent detection with noise
    demodulated_noisy = double(real(demodulated_noisy));
    demodulated_noisy = fftshift(fft(demodulated_noisy));
    demodulated_noisy(f >= bw | f <= -bw) = 0;
    demodulated_noisy = ifft(ifftshift(demodulated_noisy));
    L = length(demodulated_noisy);
    % Plot received waveform and spectrum for each SNR value
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated_noisy);
    title(['Received Signal with SNR = ' num2str(snr_dB) ' dB - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2, 1, 2);
    plot(f, abs(fftshift(fft(demodulated_noisy))) / L);
    title(['Spectrum with SNR = ' num2str(snr_dB) ' dB']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    
    demodulated_noisy_sound = resample(abs(demodulated_noisy), fs, new_fs);
    sound(abs(demodulated_noisy_sound), fs);
    pause(10);
end
%% coherent detection with frequency error
fc_error=100100;
for snr_dB = snr_values
    % Add noise to DSB_SC
    noisy_DSB_SC = awgn(DSB_SC, snr_dB);
    demodulated_noisy = noisy_DSB_SC .* cos(2*pi*fc_error*t1);
    f = new_fs / 2 * linspace(-1,1,length(DSB_SC));
    % Coherent detection with noise
    demodulated_noisy = double(real(demodulated_noisy));
    demodulated_noisy = fftshift(fft(demodulated_noisy));
    demodulated_noisy(f >= bw | f <= -bw) = 0;
    demodulated_noisy = ifft(ifftshift(demodulated_noisy));
    L = length(demodulated_noisy);
    % Plot received waveform and spectrum for each SNR value
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated_noisy);
    title(['Received Signal with SNR (FE) = ' num2str(snr_dB) ' dB - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2, 1, 2);
    plot(f, abs(fftshift(fft(demodulated_noisy))) / L);
    title(['Spectrum with SNR (FE) = ' num2str(snr_dB) ' dB']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    
    demodulated_noisy_sound = resample(abs(demodulated_noisy), fs, new_fs);
    sound(abs(demodulated_noisy_sound), fs);
    pause(10);
end
%% phase error
for snr_dB = snr_values
    % Add noise to DSB_SC
    noisy_DSB_SC = awgn(DSB_SC, snr_dB);
    demodulated_noisy = noisy_DSB_SC .* cos(2*pi*fc*t1+20);
    f = new_fs / 2 * linspace(-1,1,length(DSB_SC));
    % Coherent detection with noise
    demodulated_noisy = double(real(demodulated_noisy));
    demodulated_noisy = fftshift(fft(demodulated_noisy));
    demodulated_noisy(f >= bw | f <= -bw) = 0;
    demodulated_noisy = ifft(ifftshift(demodulated_noisy));
    L = length(demodulated_noisy);
    % Plot received waveform and spectrum for each SNR value
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated_noisy);
    title(['Received Signal with SNR (PE) = ' num2str(snr_dB) ' dB - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2, 1, 2);
    plot(f, abs(fftshift(fft(demodulated_noisy))) / L);
    title(['Spectrum with SNR (PE) = ' num2str(snr_dB) ' dB']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    
    demodulated_noisy_sound = resample(abs(demodulated_noisy), fs, new_fs);
    sound(abs(demodulated_noisy_sound), fs);
    pause(10);
end

