% MIMO-OFDM with Applied Deep Learning
% Internship - UQTR 2019

%% Randomize State & Clear Figures
rng('shuffle');
clear, figure(1), clf, figure(2), clf, figure(3), clf

%% Implementations

% FFT length ~ Number of Subcarriers
n_subcarriers = 128; 

% Number of OFDM Symbols will be transferred
n_symbols = 100;

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
x = [];
X = [];
X_without_cp = [];

% Generate Pilot value
msg_pilot = 2*(randn(1,n_pilots) > 0) - 1;

for m=1:n_symbols
    % Generate input messages, values from 0 -> M-1 for 4-bits / number
    msg_input = randi([0 QAM_Modulator-1], 1, n_data);
    msg_data = qammod(msg_input, QAM_Modulator) * A;
    Data = merge_pilot_into_msg(msg_data, msg_pilot);

    % Convert Frequence-domain to Time-domain
    time_domain_data = frq_to_time_converter(Data);

    % Add CP
    time_domain_data_with_cp = add_cp(time_domain_data, cp);
    x = [x msg_input];
    X_without_cp = [X_without_cp time_domain_data];
    X = [X time_domain_data_with_cp];
end

%% Channel Environment
% 802.11 Channel Configuration
scale=1e-9;          % nano
Ts=50*scale;         % Sampling time
t_rms=25*scale;     % RMS delay spread
num_ch=1;       % Number of channels

% Calculate Power Delay Profile
PDP=IEEE802_11_model(t_rms,Ts);  

% Generate Channel
for k=1:length(PDP)
    h(:,k) = Ray_model(num_ch).'*sqrt(PDP(k));
    avg_pow_h(k)= mean(h(:,k).*conj(h(:,k)));
end
H = fft(h(1,:),n_subcarriers);

% True channel and its time-domain length
H_power_dB = 10*log10(abs(H.*conj(H)));
channel_length = length(h);

% Channel Path by Convolution
y_channel = conv(X, h);

% Add STO & CFO
nSTO = -2;
CFO = 0.15;

y_STO = add_STO(y_channel, -nSTO);
y_STO_CFO = add_CFO(y_STO,CFO,n_subcarriers);

% Add White Gaussian Noise
y_rx = awgn(y_STO_CFO, SNR, 'measured');

%% Channel Synchronization
% STO Estimation
com_delay = n_ofdm_symbols / 2;
[STO_cor,mag_cor]= STO_by_correlation(y_rx,n_subcarriers,cp,com_delay);
[STO_dif,mag_dif]= STO_by_difference(y_rx,n_subcarriers,cp,com_delay);

% CFO Estimation
Est_CFO_CP = CFO_CP(y_rx,n_subcarriers,cp);
Est_CFO_Moose = CFO_Moose(y_rx,n_subcarriers);
Est_CFO_Classen = CFO_Classen(y_rx,n_subcarriers,cp,msg_pilot);
       
%% RX-Side

extracted_msg_total = [];
y_without_cp = [];
start = 1;

% CFO compensation in Time-Domain
y_rx = add_CFO(y_rx,-Est_CFO_CP,n_subcarriers);

for m=1:n_symbols
	% Remove CP
    y = y_rx(cp + start : n_ofdm_symbols + start - 1);
    
    % Convert to Frequency-domain
    Y = time_to_frq_converter(y);
    
    % STO compensation in Frequency-Domain
    Y = add_CFO(Y,-STO_dif,n_subcarriers);

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
    y_without_cp = [y_without_cp y];
    extracted_msg_total = [extracted_msg_total extracted_msg];
    start = start + n_ofdm_symbols;
end

% Error Calculation
err = 0 + sum(extracted_msg_total~=x, 'all');

%% Result

fprintf('Number of Error Symbol: %f\n', err);
fprintf('BER: %f\n', err / length(x)); 
% fprintf('MMSE : %6.4e\n', (H-H_est)*(H-H_est)');
% fprintf('MMSE with DFT: %6.4e\n', (H-H_DFT)*(H-H_DFT)');

%% Plot Before & After Channel
figure(1);
plot(abs(X_without_cp), 'r-o','Markersize',4,'linewidth',1); hold on; grid on;
plot(abs(y_without_cp), 'g-o','Markersize',4,'linewidth',1);
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
