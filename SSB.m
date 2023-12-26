[y, fs] = audioread('eric.wav');
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(Y));

% Plot the spectrum
figure;
subplot(2, 1, 1);
plot(f, abs(Y));
title('Spectrum of m');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Apply filter
bw = 4000;
filt = ones(size(Y));
filt(f > bw | f < -bw) = 0;
y_filter = Y .* filt;

% Plot filtered signal spectrum
subplot(2, 1, 2);
plot(f, abs(fftshift(y_filter)));
title('Filtered Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Inverse transform
y_filtered_time = ifft(ifftshift(y_filter));

% Ensure data type and scaling
y_filtered_time = real(y_filtered_time); % Ensure real values
y_filtered_time = double(y_filtered_time); % Convert to double if necessary

% Normalize if values are outside the range [-1, 1]
max_val = max(abs(y_filtered_time));
if max_val > 1
    y_filtered_time = y_filtered_time / max_val;
end

fc = 10e3;
U = 0.5;
Am = max(y_filtered_time);
Ac = Am/U; %modulationindex = Am/Ac
new_fs = 5 * fc;

resampled_signal = resample(y_filtered_time, new_fs, fs);
t1 = linspace(0, length(resampled_signal) / new_fs, length(resampled_signal));
t1 = t1';
carrier = Ac .* cos(2*pi*fc*t1);
DSB_SC = resampled_signal .* carrier;
DSB_TC = (1 + U * resampled_signal / Am) .* carrier;
DSB_SC_spectrum = fftshift(fft(DSB_SC));
f_DSB_SC = new_fs/2 * linspace(-1, 1, length(DSB_SC));

% Plot DSB-TC spectrum
figure;
subplot(1, 2, 1);
plot(f_DSB_SC, abs(fftshift(fft(DSB_TC))));
title('DSB-TC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Plot DSB-SC spectrum
subplot(1, 2, 2);
plot(f_DSB_SC, abs(DSB_SC_spectrum));
title('DSB-SC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 6. Coherent detection with no noise interference
demodulated_signal_ideal = DSB_SC .* carrier;

% Plot received waveform and spectrum
subplot(3, 3, 6);
t = linspace(0, length(demodulated_signal_ideal) / fs, length(demodulated_signal_ideal));
plot(t, demodulated_signal_ideal);
title('6. Received Signal (Coherent Detection) - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrum of received signal
Y_demodulated_ideal = fftshift(fft(demodulated_signal_ideal));
f_demodulated_ideal = linspace(-fs/2, fs/2, length(Y_demodulated_ideal));

subplot(3, 3, 9);
plot(f_demodulated_ideal, abs(Y_demodulated_ideal));
title('9. Spectrum of Received Signal (Coherent Detection)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 7. Repeat steps 5 and 6 with Butterworth filter
[b, a] = butter(4, fc / (fs / 2), 'low');
SSB_LSB_butter = filtfilt(b, a, DSB_SC);

% Plot the spectrum of SSB-LSB signal using Butterworth filter
Y_SSB_LSB_butter = fftshift(fft(SSB_LSB_butter));
f_SSB_LSB_butter = linspace(-fs/2, fs/2, length(Y_SSB_LSB_butter));

figure;
subplot(3, 3, 8);
plot(f_SSB_LSB_butter, abs(Y_SSB_LSB_butter));
title('8. Spectrum of SSB-LSB Signal (Butterworth Filter)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Coherent detection with Butterworth filter
demodulated_signal_butter = SSB_LSB_butter .* carrier;

% Plot received waveform and spectrum
subplot(3, 3, 7);
plot(t1, demodulated_signal_butter);
title('7. Received Signal (Coherent Detection with Butterworth) - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrum of received signal with Butterworth filter
Y_demodulated_butter = fftshift(fft(demodulated_signal_butter));
f_demodulated_butter = linspace(-fs/2, fs/2, length(Y_demodulated_butter));

subplot(3, 3, 9);
hold on;
plot(f_demodulated_butter, abs(Y_demodulated_butter), 'r');
legend('Ideal Filter', 'Butterworth Filter');
hold off;

% 8. For the ideal filter case, get the received signal again with noise
% when noise is added to SSB-SC with SNR = 0, 10, and 30
snr_values = [0, 10, 30];
for snr_dB = snr_values
    % Add noise to SSB-SC
    noisy_SSB_SC = awgn(SSB_LSB, snr_dB);
    
    % Coherent detection with noise
    demodulated_noisy_signal = noisy_SSB_SC .* carrier;

    % Plot received waveform and spectrum for each SNR value
    figure;
    subplot(2, 1, 1);
    plot(t, demodulated_noisy_signal);
    title(['Received Signal with SNR = ' num2str(snr_dB) ' dB - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2, 1, 2);
    plot(f_demodulated_ideal, abs(fftshift(fft(demodulated_noisy_signal))));
    title(['Spectrum with SNR = ' num2str(snr_dB) ' dB']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end

% 9. For the ideal filter case, generate a SSB-TC
SSB_TC = carrier + SSB_LSB;

% Envelope detection
envelope_SSB_TC = abs(hilbert(SSB_TC));

% Playback received message
sound(envelope_SSB_TC, fs);

% Plot received waveform
figure;
plot(t, SSB_TC, 'b', t, -envelope_SSB_TC, 'r', t, envelope_SSB_TC, 'r', 'Linewidth', 1.5);
title('Received Message with Envelope Detection');
xlabel('Time (s)');
ylabel('Amplitude');
legend('SSB-TC', 'Envelope');

% Spectrum of received signal
Y_SSB_TC = fftshift(fft(SSB_TC));
f_SSB_TC = linspace(-fs/2, fs/2, length(Y_SSB_TC));

figure;
plot(f_SSB_TC, abs(Y_SSB_TC));
title('Spectrum of Received SSB-TC Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');