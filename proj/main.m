% MIMO-OFDM with Applied Deep Learning
% Internship - UQTR 2019

%% Randomize State & Clear Figures
rng('shuffle');
clf('reset');

%% Implementations

% FFT length ~ Number of Subcarriers
n_subcarriers = 128; 

% Number of OFDM Symbols will be transferred
n_symbols = 1;

% CP
n_cps = 4;
cp = n_subcarriers / n_cps;

% Number of pilot carriers
pilot_distance = 8;
n_pilots = n_subcarriers / pilot_distance;

% Number of data carriers
n_data = n_subcarriers - n_pilots;

% Number of OFDM Symbols
n_ofdm_symbols = n_subcarriers + cp;

% Number of bits / OFDM Symbols
n_bit_per_symbol = 4;
QAM_Modulator = 2^n_bit_per_symbol;

% Carriers 
all_carriers = 1:1:n_subcarriers;
pilot_carriers = 1:pilot_distance:n_subcarriers;

% Signal-to-Noise Ratio
SNR = 25;

% Number of error symbols
n_error_symbols = 0;

% Signal Energy & QAM Normalization Factor
Es = 1;
A = sqrt(3/2/(QAM_Modulator-1)*Es);

%% TX-Side         
% Generate Pilot value
msg_pilot = 2*(randn(1,n_pilots) > 0) - 1;

% Generate input messages, values from 0 -> M-1 for 4-bits / number
msg_input = randi([0 QAM_Modulator-1], 1, n_data);
msg_data = qammod(msg_input, QAM_Modulator) * A;
Data = merge_pilot_into_msg(msg_data, msg_pilot);

% Convert Frequence-domain to Time-domain
time_domain_data = frq_to_time_converter(Data);

% Add CP
time_domain_data_with_cp = add_cp(time_domain_data, cp);

%% Channel Environment

% Generate 2-tap channel
% h = [(rand + 1i * rand) (rand + 1i * rand)/2];
h = [1+0j 0+0j 0.3+0.3j];

% True channel and its time-domain length
H = fft(h, n_subcarriers);
H_power_dB = 10*log10(abs(H.*conj(H)));
channel_length = length(h);

% Channel Path by Convolution
y_channel = conv(time_domain_data_with_cp, h);

% Add White Gaussian Noise
yt = awgn(y_channel, SNR, 'measured');

%% RX-Side
% Remove CP
y = yt(cp + 1 : n_ofdm_symbols);

% Convert to Frequency-domain
Y = time_to_frq_converter(y);

% Channel Estimation by MMSE
H_est = MMSE_CE(Y, msg_pilot, pilot_carriers, n_subcarriers, pilot_distance, h, SNR);

% Calculate Channel Estimation Power
H_est_power_dB = 10 * log10(abs(H_est.*conj(H_est)));

% DFT-Based Channel Estimation
h_est = ifft(H_est);
h_DFT = h_est(1:channel_length); 
H_DFT = fft(h_DFT, n_subcarriers);
H_DFT_power_dB = 10 * log10(abs(H_DFT.*conj(H_DFT)));


% Equalization
Y_eq = channel_equalizer(Y, H_est);

% Extract Data
extracted_Data = extract_data(Y_eq, n_pilots);

% Get Extracted Message
extracted_msg = qamdemod(extracted_Data/A, QAM_Modulator);

% Noise Calculation
noise = 0 + sum(extracted_msg~=msg_input, 'all');

%% Result

fprintf('Number of Error Symbol: %f\n', noise);
fprintf('BER: %f\n', noise / length(msg_input)); 
fprintf('MMSE : %6.4e\n', (H-H_est)*(H-H_est)');
fprintf('MMSE with DFT: %6.4e\n', (H-H_DFT)*(H-H_DFT)');

%% Plot Before & After Channel
figure(1);
plot(abs(time_domain_data), 'r-o','Markersize',4,'linewidth',1); hold on; grid on;
plot(abs(y), 'g-o','Markersize',4,'linewidth',1);
xlabel('Time');
ylabel('$x|t|$');
legend('True Channel','After Channel Env.','SouthEast');  set(gca,'fontsize',9)
signal_power = mean(abs(y_channel.^2));
sigma2 = signal_power * 10^(-SNR/10);
title(['Before & After Propagation over Chanel. RX Signal Power: ' num2str(signal_power) '. Noise Power: ' num2str(sigma2)]);

%% Plot Channel Estimation Power
method = 'MMSE';
figure(2), subplot(211), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
plot(H_est_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 n_subcarriers -10 10])
title(method); xlabel('Subcarrier Index'); ylabel('Power[dB]');
legend('True Channel',method,'SouthEast');  set(gca,'fontsize',9)
subplot(212), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
plot(H_DFT_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 n_subcarriers -10 10])
title([method ' with DFT']); xlabel('Subcarrier Index'); ylabel('Power[dB]');
legend('True Channel',[method ' with DFT'],'SouthEast');  set(gca,'fontsize',9)

%% Plot recevied value on Interpolation
figure(3), subplot(221), plot(Y,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
axis('equal'), set(gca,'fontsize',9), hold on,
title('Received data before equalization')
subplot(222), plot(Y_eq,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
axis('equal'), set(gca,'fontsize',9), hold on,
title('Received data after equalization');