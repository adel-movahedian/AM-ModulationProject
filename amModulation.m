clc;clear;

% Define m(t)
t = -10:0.01:10;
m = exp(-6*t).*(heaviside(t-4) - heaviside(t-8)) + exp(6*t).*(heaviside(-t-4) - heaviside(-t-8));

% Plotting in the time domain
figure;
plot(t, m, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Signal m(t) in Time Domain', 'Interpreter', 'latex');
grid on;

M = fftshift(fft(m));

% Define the frequency range
fs = 1/0.01; % Sampling frequency
f = (0:length(M)-1)*fs/length(M);
M_mag = abs(M);

% Find the index corresponding to the maximum magnitude of M_mag
[maxMag, maxIndex] = max(M_mag);
bandwidth = 2 * abs(f(maxIndex));
disp(['Signal Bandwidth: ', num2str(bandwidth), ' Hz']);

% Plot the magnitude of the signal as a function of frequency
figure;
plot(f, abs(M), 'r')
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Magnitude', 'Interpreter', 'latex')
title('Fourier transform of the signal in frequency domain', 'Interpreter', 'latex')
grid on;

% Calculate the maximum frequency component in the signal
f_max = bandwidth / 2;

% Determine the desired sampling rate
Fs_desired = 2 * f_max;

% Resample the signal to the desired sampling rate
t_resampled = -10:1/Fs_desired:10; % Resampled time axis
m_resampled = exp(-6*t_resampled).*(heaviside(t_resampled-4) - heaviside(t_resampled-8)) + exp(6*t_resampled).*(heaviside(-t_resampled-4) - heaviside(-t_resampled-8)); % Resampled signal

% Plot the resampled signal in the time domain
figure;
plot(t_resampled, m_resampled, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Resampled Signal m(t) in Time Domain', 'Interpreter', 'latex');
grid on;

% Define the cutoff frequency for the low-pass filter
cutoff_frequency = 0.2 * bandwidth; % Adjust the factor (0.2) as needed

% Design a Butterworth low-pass filter
filter_order = 4;
[b, a] = butter(filter_order, cutoff_frequency/(Fs_desired/2));

% Apply the filter to the resampled signal
filtered_signal = filter(b, a, M);

figure;
plot(t, filtered_signal, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Resaiooooooooooooooooin', 'Interpreter', 'latex');
grid on;

%% AM Modulation

Fc = 200; 
Fs = 2000; 
y = ammod(m,Fc,Fs); 

% Plot the signals
figure;
plot(t, y,color = 'red')
title('Message signal')
xlabel('Time (s)')
ylabel('Amplitude')
figure;
plot(t, abs(fftshift(fft(y))),color = 'blue')
title('AM modulated signal', 'Interpreter', 'latex')
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')

%% AM Demodulation
demodulated = coherentAMDemodulation(y , Fc , Fs , 200)
plot(t,demodulated,color = 'r')
title('AM demodulated signal', 'Interpreter', 'latex')
xlabel('Time (s)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
%% AM Envelope Demodulation
[envelopeDemodulatedUp , envelopeDemodulatedLow] = envelope_detector(y , Fs);
plot(t, envelopeDemodulatedUp,color = 'green');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Upper Envelope Component', 'Interpreter', 'latex');
figure
plot(t, envelopeDemodulatedLow,color = 'black');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Lower Envelope Component', 'Interpreter', 'latex');
%% my ussb and lssb
clc;

% Define m(t)
t = -10:0.01:10;
m = exp(-6*t).*(heaviside(t-4) - heaviside(t-8)) + exp(6*t).*(heaviside(-t-4) - heaviside(-t-8));

% Plotting in the time domain
figure;
plot(t, m, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Signal m(t) in Time Domain', 'Interpreter', 'latex');
grid on;

M = fftshift(fft(m));

% Define the frequency range
fs = 1/0.01; % Sampling frequency
f = (0:length(M)-1)*fs/length(M);
M_mag = abs(M);

% Find the index corresponding to the maximum magnitude of M_mag
[maxMag, maxIndex] = max(M_mag);
bandwidth = 2 * abs(f(maxIndex));
disp(['Signal Bandwidth: ', num2str(bandwidth), ' Hz']);

% Plot the magnitude of the signal as a function of frequency
figure;
plot(f, abs(M), 'r')
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Magnitude', 'Interpreter', 'latex')
title('Fourier transform of the signal in frequency domain', 'Interpreter', 'latex')
grid on;

% Calculate the maximum frequency component in the signal
f_max = bandwidth / 2;

% Determine the desired sampling rate
Fs_desired = 2 * f_max;

% Resample the signal to the desired sampling rate
t_resampled = -10:1/Fs_desired:10; % Resampled time axis
m_resampled = exp(-6*t_resampled).*(heaviside(t_resampled-4) - heaviside(t_resampled-8)) + exp(6*t_resampled).*(heaviside(-t_resampled-4) - heaviside(-t_resampled-8)); % Resampled signal

% Plot the resampled signal in the time domain
figure;
plot(t_resampled, m_resampled, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Resampled Signal m(t) in Time Domain', 'Interpreter', 'latex');
grid on;

% Define the cutoff frequency for the low-pass filter
cutoff_frequency = 0.2 * bandwidth; % Adjust the factor (0.2) as needed

% Design a Butterworth low-pass filter
filter_order = 4;
[b, a] = butter(filter_order, cutoff_frequency/(Fs_desired/2));

% Apply the filter to the resampled signal
filtered_signal = filter(b, a, m_resampled);

figure;
plot(t_resampled, filtered_signal, 'b');
xlabel('Time', 'Interpreter', 'latex');
ylabel('m(t)', 'Interpreter', 'latex');
title('Resampled Signal after Low-Pass Filtering', 'Interpreter', 'latex');
grid on;

%% USSB Modulation
Fc = 200; 
Fs = 2000; 

% Ensure the filtered signal is real-valued
filtered_signal = real(filtered_signal);

% Perform USSB modulation
ussb_modulated_signal = ammod(filtered_signal, Fc, Fs);

% Create time vector for USSB modulated signal
t_ussb = linspace(-10, 10, length(ussb_modulated_signal));

% Plot the USSB modulated signal
figure;
plot(t_ussb, ussb_modulated_signal, 'r');
title('USSB Modulated Signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
%% AM Demodulation for ussb
demodulated = coherentAMDemodulation(ussb_modulated_signal , Fc , Fs , 200);
plot(t_resampled,demodulated,color = 'r');
title('AM demodulated ussb signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
%% LSSB Modulation

lssb_modulated_signal = ssbmod(filtered_signal, Fc, Fs, 0);

% Plot the LSSB modulated signal
figure;
plot(t_resampled, lssb_modulated_signal, 'b');
title('LSSB Modulated Signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;

%% AM Demodulation for lssb
demodulated = coherentAMDemodulation(lssb_modulated_signal , Fc , Fs , 200);
plot(t_resampled,demodulated,color = 'r');
title('AM demodulated lssb signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;

%% transfer
clc
% Load the audio file
[y_audio, Fs_audio] = audioread('audio_2024-02-02_21-04-52.ogg');

% Define modulation parameters
Fc = 200; % Carrier Frequency
Fs = 2 * Fs_audio; % Modulated Sampling Frequency
cutoff_frequency = 4000; % Bandwidth

% AM modulation
y_modulated = ammod(y_audio, Fc, Fs);

% Generate time axis for audio signal
t_audio = (0:length(y_audio)-1) / Fs_audio;

% Generate time axis for modulated signal
t_modulated = (0:length(y_modulated)-1) / Fs;

% Plot audio signal in time domain
figure;
plot(t_audio, y_audio,'r');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Audio Signal in Time Domain', 'Interpreter', 'latex');
grid on;
% Plot modulated signal in time domain
figure;
plot(t_modulated, y_modulated,'b');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Modulated Signal in Time Domain', 'Interpreter', 'latex');
grid on;
% Compute frequencies present in the modulated signal
Y_modulated_fft = fftshift(fft(y_modulated));
f_modulated = (-length(Y_modulated_fft)/2:length(Y_modulated_fft)/2-1) * (Fs / length(Y_modulated_fft));

% Plot frequency spectrum of modulated signal
figure;
plot(f_modulated, abs(Y_modulated_fft),'g');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
title('Frequency Spectrum of Modulated Signal', 'Interpreter', 'latex');
grid on;
%%
clc
demodulated_coherent = coherentAMDemodulation2(y_modulated, Fc, Fs);
% Demodulate the modulated signal using Envelope Detection
demodulated_envelope = envelopeAMDemodulation2(y_modulated);
% Plot the original signal and the demodulated signals
figure;
subplot(3,1,1);
plot(t_audio, y_audio,'black');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Original Audio Signal', 'Interpreter', 'latex');
grid on;
subplot(3,1,2);
plot(t_audio, demodulated_coherent,'r');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Demodulated Signal (Coherent Detection)', 'Interpreter', 'latex');
grid on;
subplot(3,1,3);
plot(t_audio, demodulated_envelope, 'b');
hold on;
plot(t_audio, -demodulated_envelope, 'r');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Demodulated Signal (Envelope Detection)', 'Interpreter', 'latex');
legend('Positive Envelope', 'Negative Envelope', 'Interpreter', 'latex');
grid on;
%% play
clc
[y_audio, Fs_audio] = audioread('audio_2024-02-02_21-04-52.ogg');
sound(y_audio, Fs_audio); %original
pause(5);
sound(demodulated_coherent, Fs_audio); % coherent
pause(5);
sound(demodulated_envelope, Fs_audio); %envelope

%% functions

function y = heavisideplus(x)
    y = (x >= 0);
end

function demodulated_signal = coherentAMDemodulation(received_signal, carrier_frequency, fs, cutoff_frequency)
t = 0:1/fs:length(received_signal)/fs-1/fs;
local_oscillator = cos(2*pi*carrier_frequency*t);
demodulated_signal = received_signal .* local_oscillator;
[b, a] = butter(4, cutoff_frequency/(fs/2), 'low');
demodulated_signal = filter(b, a, demodulated_signal);
demodulated_signal = demodulated_signal / max(abs(demodulated_signal));
end

function [yupper, ylower] = envelope_detector(x, fs)
% yupper: upper envelope
% ylower: lower envelope
x = x - mean(x);
x_analytic = hilbert(x);
x_magnitude = abs(x_analytic);
fc = 50;
[b, a] = butter(2, fc / (fs / 2));
yupper = filter(b, a, x_magnitude);
ylower = -yupper;
end
% Demodulation using Coherent Detection
function demodulated_coherent = coherentAMDemodulation2(y_modulated, Fc, Fs)
    % Generate time axis
    t = (0:length(y_modulated)-1) / Fs;
    
    % Initialize demodulated signal
    demodulated_coherent = zeros(size(y_modulated));
    
    % Generate local oscillator
    local_oscillator = cos(2*pi*Fc*t);
    
    % Demodulate the signal element-wise
    for i = 1:length(y_modulated)
        demodulated_coherent(i) = y_modulated(i) * local_oscillator(i);
    end
end


% Envelope Detection
function demodulated_envelope = envelopeAMDemodulation2(y_modulated)
    % Convert to absolute value to get the envelope
    demodulated_envelope = abs(y_modulated);
end