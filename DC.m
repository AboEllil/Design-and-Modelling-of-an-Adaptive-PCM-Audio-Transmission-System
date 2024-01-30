% Reset environment
close all;
clc;

% User inputs
audio_file = input('Enter your Audio File Name (with extension): ', 's');
BW_in_kHz = input('Input signal bandwidth (kHz): ');
SNR_dB = input('Input SNR (dB): ');
A_law_compressor_parameter = input('A-law compressor parameter (A): ');

% Convert user inputs to appropriate units
BW = BW_in_kHz * 1000; % Convert to Hz

% Start of code
[x, fs] = audioread(audio_file); % Read audio file

% Calculate bit resolution based on Nyquist condition
l = floor(BW / fs);
if l < 1
    error('Error! Bandwidth must be greater than sampling frequency!');
end

% Calculate bit rate and bandwidth efficiency
Rb = fs * l;
BW_efficiency = (Rb / BW) * 100;

% Display calculated parameters
disp(['Bit Resolution (l): ' num2str(l)]);
disp(['Bit Rate (Rb): ' num2str(Rb)]);
disp(['Bandwidth Efficiency (%): ' num2str(BW_efficiency)]);

% Perform A-law compression
A_law_compressed_signal = compand(x, A_law_compressor_parameter, max(abs(x)), 'A/compressor');

% Calculate quantization parameters
max_amplitude = max(abs(A_law_compressed_signal));
number_of_levels = 2^l;
quantization_step = 2 * max_amplitude / number_of_levels;

% Quantize the A-law compressed signal
quantized_signal = A_law_compressed_signal;

% Create an empty matrix to store encoded data
encoded_data = zeros(length(A_law_compressed_signal), l);

% Loop through each sample and perform quantization and encoding
for i = 1:length(A_law_compressed_signal)
    quantization_level = floor((A_law_compressed_signal(i) + max_amplitude) / quantization_step);
    
    % Clamp quantization level to avoid exceeding maximum level
    if quantization_level == number_of_levels
        quantization_level = number_of_levels - 1;
    end
    
    % Quantize the sample and store it
    quantized_signal(i) = -max_amplitude + (quantization_level + 0.5) * quantization_step;
    
    % Convert quantization level to binary representation and store it
    encoded_data(i, :) = de2bi(quantization_level, l);
    
    % Apply bipolar encoding
    if all(encoded_data(i, :) == 0)
        encoded_data(i, :) = -1;
    end
end

% Calculate noise power
noise_power = 1 * 10^-6;

% Calculate pulse amplitude based on SNR
pulse_amplitude = sqrt(2 * noise_power * BW * 10^(SNR_dB / 10));

% Add scaled encoded data to AWGN
received_signal = pulse_amplitude * encoded_data + sqrt(noise_power * BW) * randn(size(encoded_data));

% Perform pulse detection
detected_data = received_signal > (pulse_amplitude / 2);

% Decode the received data
received_quantization_levels = bi2de(detected_data);

% Reconstruct the quantized signal
reconstructed_quantized_signal = -max_amplitude + (received_quantization_levels + 0.5) * quantization_step;

% Perform A-law expansion
reconstructed_signal = compand(reconstructed_quantized_signal, A_law_compressor_parameter, max(abs(reconstructed_quantized_signal)), 'A/expander');

% Calculate SQNR and NMSE
quantized_reconstructed_signal = compand(quantized_signal, A_law_compressor_parameter, max(abs(quantized_signal)), 'A/expander');
noise_power_quantization = mean((x - quantized_reconstructed_signal').^2);
signal_power = mean((x).^2);
SQNR = signal_power / noise_power_quantization;
NMSE = (sum((x - reconstructed_signal').^2)) / sum((x).^2);

% Display calculated metrics
disp(['NMSE: ' num2str(NMSE)]);
disp(['SQNR: ' num2str(SQNR)]);
disp(['SQNR (dB): ' num2str(10 * log10(SQNR))]);
 sound (reconstructed_signal, fs ) 
