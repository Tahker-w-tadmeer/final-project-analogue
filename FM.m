% 1. Read the audio file and find the spectrum
[audio, Fs] = audioread('eric.wav');
N = length(audio);
t = (0:N-1) / Fs;

% Calculate the spectrum
audio_spectrum_shifted = fftshift(fft(audio));

% Plot the spectrum
figure;
plot(linspace(-Fs/2, Fs/2, N), abs(audio_spectrum_shifted));
title('Original Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 2. Ideal filter to remove frequencies above 4 KHz
cutoff_frequency = 4000;
filter = (abs(linspace(-Fs/2, Fs/2, N)) <= cutoff_frequency);
audio_spectrum_filtered = audio_spectrum_shifted .* filter';

% 3. Obtain filtered signal in time and frequency domain
audio_filtered =real(ifft(ifftshift(audio_spectrum_filtered)));

% Plot the filtered spectrum
figure;
plot(linspace(-Fs/2, Fs/2, N), abs(audio_spectrum_filtered));
title('Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Plot the filtered signal in time domain
figure;
plot(t, audio_filtered);
title('Filtered Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% 4. Sound the filtered audio signal
sound(abs(audio_filtered), Fs);

% 5. Generate NBFM signal
kf = 0.2/(2*pi*max(abs(cumsum(audio_filtered)))./Fs);
Ac = 1;
Fc = 100000; % Carrier frequency
resampled_audio = resample(audio_filtered, 5*Fc, Fs);
Fs_nbfm = 5 * Fc; % Sampling frequency for NBFM

%calculate time vector
tstart = 0;
tend = tstart + length(resampled_audio) / Fs_nbfm;
t1 = linspace(tstart, tend, length(resampled_audio));
t1 = t1';

%FM modulated signal
NBFM_signal = Ac * cos(2*pi*Fc*t1 + 2*pi*kf*cumsum(resampled_audio)./Fs_nbfm);

% Plot the NBFM spectrum
L = length(NBFM_signal);
NBFM_spectrum_shifted = real(fftshift(fft(NBFM_signal)));
f = Fs/2*linspace(-1,1,L);
figure;
plot(f, NBFM_spectrum_shifted/L)
title('NBFM Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% 6. Condition for NBFM: Frequency deviation (delta_f) should be much smaller than the message bandwidth (B)
B = cutoff_frequency;
delta_f = 75; % Frequency deviation
condition = delta_f < B;
disp(['Condition for NBFM: ', num2str(condition)]);

% 7. Demodulate NBFM signal
% Discriminator
dy = diff(NBFM_signal);
dy = [0; dy];

% envelope detector
demodulated_NBFM = abs(hilbert(dy)) -  mean(abs(hilbert(dy)));

% Plot the time-domain signal
figure;
plot(t1,demodulated_NBFM);
title('Demodulated NBFM Time-Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-2*10^-4 2*10^-4]);