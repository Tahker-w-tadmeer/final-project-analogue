[y, fs] = audioread('eric.wav');
Y = fftshift(fft(y));
f = linspace(-fs/2, fs/2, length(Y));

% Plot the spectrum
plot(f, abs(Y));
title('Spectrum of m');
figure();

% Apply filter
bw=4000;
filt = ones(size(Y));
filt(f > bw|f<-bw) = 0;
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
fc=10e3;
U=0.5;
Am=max(y_filtered_time);
Ac=Am/U;
new_fs=5*fc;
resampled_signal=resample(y_filtered_time,new_fs,fs);
t1=linspace(0,length(resampled_signal)/new_fs,length(resampled_signal));
t1=t1';
carrier=Ac.*cos(2*pi*fc*t1);
DSB_SC=resampled_signal.*carrier;
DSB_TC=(1+U*resampled_signal/Am).*carrier;
DSB_SC_spectrum = fftshift(fft(DSB_SC));
f_DSB_SC=new_fs/2*linspace(-1,1,length(DSB_SC));
% Plot DSB-TC spectrum
figure();
subplot(1,2,1)
plot(f_DSB_SC, abs(DSB_TC_spectrum));
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
%envelop detection for DSB_SC
envelope_DSB_SC = abs(hilbert(DSB_SC));
figure; 
plot(t1, DSB_SC);
hold on;
plot(t1,-envelope_DSB_SC,'r-',t1, envelope_DSB_SC,'-r','Linewidth',1.5); 
hold off;
title('DSB_sc with envellop ');
ylim([-2 2])
xlim([2 2.5])
