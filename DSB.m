[y, fs] = audioread('eric.wav');
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(Y));

% Plot the spectrum
plot(f, abs(Y));
title('Spectrum of m');
figure();

% Apply filter
bw = 4000;
filt = ones(size(Y));
filt(f > bw|f<-bw) = 0;
y_filter = Y .* filt;
plot(f, y_filter);
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
DSB_TC = (1 + U * resampled_signal / Am) .* carrier;
DSB_SC_spectrum = fftshift(fft(DSB_SC));
f_DSB_SC = new_fs/2*linspace(-1, 1, length(DSB_SC));


%% Plot DSB-TC spectrum
figure();
subplot(1,2,1)
plot(f_DSB_SC, abs(DSB_SC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-TC Modulated Signal Spectrum');
DSB_TC_spectrum = fftshift(fft(DSB_TC));
f_DSB_TC=new_fs/2*linspace(-1,1,length(DSB_TC));

% Plot DSB-SC spectrum
subplot(1,2,2)
plot(f_DSB_TC, abs(DSB_SC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-SC Modulated Signal Spectrum');

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
%% plot demodi=ulation signal
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
snr = [0, 10, 30];
for i = 1:length(snr)
    noisy_signal = awgn(DSB_SC, snr(i));
    
    % Demodulation using coherent detection
    demodulated = noisy_signal .* cos(2 * pi * fc * t1);
    demodulated_FFT = fftshift(fft(demodulated));
    demodulated_FFT(f >= bw | f <= -bw) = 0;
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot and playback demodulated signal
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated); 
    title(['SNR demodulated signal in time domain = ' num2str(snr(i)) ' dB']);
    
    % Plot demodulated signal in frequency domain
    F = fftshift(fft(demodulated));
    f = new_fs / 2 * linspace(-1, 1, length(demodulated));
    subplot(2, 1, 2);
    plot(f, abs(F) / length(demodulated));
    title(['SNR demodulated signal in frequency domain = ' num2str(snr(i)) ' dB']);
    % Sound playback of the demodulated signal
    sound(abs(demodulated), new_fs); 
    pause(5); 
end
%% coherent detection with frequency error
fc_error=100100;
for i = 1:length(snr)
    noisy_signal = awgn(DSB_SC, snr(i));
    
    % Demodulation using coherent detection
    demodulated = noisy_signal .* cos(2 * pi * fc_error * t1);
    demodulated_FFT = fftshift(fft(demodulated));
    demodulated_FFT(f >= bw | f <= -bw) = 0;
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot and playback demodulated signal
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated); 
    title(['SNR demodulated signal in time domain = ' num2str(snr(i)) ' dB']);
    
    % Plot demodulated signal in frequency domain
    F = fftshift(fft(demodulated));
    f = new_fs / 2 * linspace(-1, 1, length(demodulated));
    subplot(2, 1, 2);
    plot(f, abs(F) / length(demodulated));
    title(['SNR demodulated signal in frequency domain = ' num2str(snr(i)) ' dB']);
    % Sound playback of the demodulated signal
    sound(abs(demodulated), new_fs); 
    pause(5); 
end
%% phase error
for i = 1:length(snr)
    noisy_signal = awgn(DSB_SC, snr(i));
    
    % Demodulation using coherent detection
    demodulated = noisy_signal .* cos(2 * pi * fc * t1 +pi/9);
    demodulated_FFT = fftshift(fft(demodulated));
    demodulated_FFT(f >= bw | f <= -bw) = 0;
    demodulated = ifft(ifftshift(demodulated_FFT));
    
    % Plot and playback demodulated signal
    figure;
    subplot(2, 1, 1);
    plot(t1, demodulated); 
    title(['SNR demodulated signal in time domain = ' num2str(snr(i)) ' dB']);
    
    % Plot demodulated signal in frequency domain
    F = fftshift(fft(demodulated));
    f = new_fs / 2 * linspace(-1, 1, length(demodulated));
    subplot(2, 1, 2);
    plot(f, abs(F) / length(demodulated));
    title(['SNR demodulated signal in frequency domain = ' num2str(snr(i)) ' dB']);
    % Sound playback of the demodulated signal
    sound(abs(demodulated), new_fs); 
    pause(5); 
end


