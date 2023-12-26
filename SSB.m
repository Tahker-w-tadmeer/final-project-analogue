%% 1-2. Same as DSB

[y, fs] = audioread('eric.wav');
L = length(y);
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, L);

% Plot the spectrum
figure;
subplot(2, 1, 1);
plot(f, abs(Y) / L);
title('Original Spectrum of m');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Apply filter
bw = 4000;
Y(f >= bw | f <= -bw) = 0;

subplot(2, 1, 2);
plot(f, abs(Y) / L);
title('Filtered Spectrum of m');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% 3. Filtered signal in time domain (Inverse transform)
y_filtered_time = ifft(ifftshift(Y));

t1 = linspace(0, length(y_filtered_time) / fs, length(y_filtered_time));
t1 = t1';

y_filtered_time = real(double(y_filtered_time));

figure; plot(t1, y_filtered_time);
title('Filtered Signal of m Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

fc = 100000;
U = 0.5;
Am = max(y_filtered_time);
Ac = Am/U; % modulationindex = Am/Ac
new_fs = 5 * fc;


%% 4. Plot DSB-SC spectrum
message_for_sound = resample(y_filtered_time, new_fs, fs);
% sound(abs(message_for_sound), fs);

message = resample(y_filtered_time, new_fs, fs);
t1 = linspace(0, length(message) / new_fs, length(message));
t1 = t1';

L = length(message);
carrier = Ac .* cos(2*pi*fc*t1);
DSB_SC = message .* carrier;
DSB_SC_spectrum = fftshift(fft(DSB_SC));
f_DSB_SC = new_fs/2 * linspace(-1, 1, length(DSB_SC));

figure;
subplot(2, 1, 1);
plot(f_DSB_SC, abs(DSB_SC_spectrum) / L);
title('DSB-SC Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 5. Obtain SSB_LSB
SSB_LSB = DSB_SC;
L = length(SSB_LSB);
f = new_fs / 2 * linspace(-1, 1, L);
F = fftshift(fft(SSB_LSB));
F(f>=fc | f<=-fc) = 0;
SSB_LSB = ifft(ifftshift(F));

% Plot SSB_LSB spectrum
subplot(2, 1, 2);
plot(f, abs(F) / L);
title('SSB LSB Modulated Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%% 6. Coherent detection with no noise interference
demodulated_signal_ideal = SSB_LSB .* cos(2*pi*fc*t1);

demodulated_signal_ideal_sound = resample(abs(demodulated_signal_ideal), fs, new_fs);
sound(abs(demodulated_signal_ideal_sound), fs);

demodulated_signal_ideal = fftshift(fft(demodulated_signal_ideal));
demodulated_signal_ideal(f >= fc | f <= -fc) = 0;
demodulated_signal_ideal = ifft(ifftshift(demodulated_signal_ideal));

% Plot received waveform and spectrum
subplot(2, 1, 1);
plot(t1, demodulated_signal_ideal);
title('Received Signal (CD) in Time domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrum of received signal
L = length(demodulated_signal_ideal);
Y_demodulated_ideal = fftshift(fft(demodulated_signal_ideal));
f_demodulated_ideal = linspace(-new_fs/2, new_fs/2, L);

subplot(2, 1, 2);
plot(f_demodulated_ideal, abs(Y_demodulated_ideal) / L);
title('Spectrum of Received Signal (CD)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10000 10000]);
%% 7. Repeat steps 5 and 6 with Butterworth filter
[b, a] = butter(4, fc / (new_fs / 2), 'low');
SSB_LSB_butter = filtfilt(b, a, demodulated_signal_ideal);
L = length(SSB_LSB_butter);

% Plot the spectrum of SSB-LSB signal using Butterworth filter
Y_SSB_LSB_butter = fftshift(fft(SSB_LSB_butter));
f_SSB_LSB_butter = linspace(-fs/2, fs/2, length(Y_SSB_LSB_butter));

figure;
subplot(3, 1, 1);
plot(f_SSB_LSB_butter, abs(Y_SSB_LSB_butter) / L);
title('Spectrum of SSB-LSB Signal (Butterworth Filter)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Coherent detection with Butterworth filter
demodulated_signal_butter = SSB_LSB_butter .* carrier;
t = linspace(0, length(demodulated_signal_ideal) / fs, length(demodulated_signal_ideal));
% Plot received waveform and spectrum
subplot(3, 1, 2);
plot(t, demodulated_signal_butter);
title('Received Signal (Coherent Detection with Butterworth) - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrum of received signal with Butterworth filter
Y_demodulated_butter = fftshift(fft(demodulated_signal_butter));
f_demodulated_butter = linspace(-fs/2, fs/2, length(Y_demodulated_butter));

subplot(3, 1, 3);
plot(f_demodulated_butter, abs(Y_demodulated_butter) / L);
title('Spectrum of Received Signal (Coherent Detection with Butterworth)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% 8. For the ideal filter case, get the received signal again with noise
% when noise is added to SSB-SC with SNR = 0, 10, and 30
snr_values = [0, 10, 30];
for snr_dB = snr_values
    % Add noise to SSB-SC
    noisy_SSB_SC = awgn(SSB_LSB, snr_dB);
    demodulated_noisy = noisy_SSB_SC .* cos(2*pi*fc*t1);
    
    % Coherent detection with noise
    demodulated_noisy = double(real(demodulated_noisy));
    demodulated_noisy = fftshift(fft(demodulated_noisy));
    demodulated_noisy(f >= bw | f <= -bw) = 0;
    demodulated_noisy = ifft(ifftshift(demodulated_noisy));
    
    f = new_fs / 2 * linspace(-1,1,L);
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
end
%% 9. For the ideal filter case, generate a SSB-TC
SSB_TC = carrier + SSB_LSB;
L = length(SSB_TC);

% Spectrum of received signal
Y_SSB_TC = fftshift(fft(SSB_TC));
f_SSB_TC = linspace(-new_fs/2, new_fs/2, length(Y_SSB_TC));

figure;
plot(f_SSB_TC, abs(Y_SSB_TC) / L);
title('Spectrum of Received SSB-TC Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Envelope detection
envelope_SSB_TC = abs(hilbert(SSB_TC));

% Plot received waveform
figure;
plot(t1, SSB_TC, 'b', t1, -envelope_SSB_TC, 'r', t1, envelope_SSB_TC, 'r', 'Linewidth', 1.5);
title('Received Message with Envelope Detection');
xlabel('Time (s)');
ylabel('Amplitude');
legend('SSB-TC', 'Envelope');
ylim([-5 5]);
xlim([3 3.5]);

envelope_SSB_TC = resample(envelope_SSB_TC, fs, new_fs);
sound(abs(envelope_SSB_TC), fs);