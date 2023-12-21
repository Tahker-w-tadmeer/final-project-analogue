[y, fs] = audioread('eric.wav');
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(Y));

% Plot the spectrum
plot(f, abs(Y));
title('Spectrum of m');
figure();

% Apply filter
filt = ones(size(Y));
filt(f > 4000) = 0;
y_filter = Y .* filt;
plot(f,y_filter);
title('filtered signal spectrum');
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

sound(y_filtered_time, fs);
Fc = 100000; 
new_fs = 5 * Fc; 
resampled_signal = resample(y_filtered_time, new_fs, fs);
t = 0:1/new_fs:(length(resampled_signal)-1)/new_fs;
carrier = cos(2*pi*Fc*t);
DSB_TC = resampled_signal .* (1 + 0.5 * cos(2*pi*Fc*t));
DSB_SC = resampled_signal .* cos(2*pi*Fc*t);
DSB_TC_spectrum = fftshift(fft(DSB_TC));
DSB_SC_spectrum = fftshift(fft(DSB_SC));
subplot(2, 1, 1);
plot(frequencies, abs(DSB_TC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-TC Modulated Signal Spectrum');

subplot(2, 1, 2);
plot(frequencies, abs(DSB_SC_spectrum));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DSB-SC Modulated Signal Spectrum');