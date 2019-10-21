%% Randomize State & Clear Figures
rng('shuffle');
clear,clc;
% clear, figure(1), clf, figure(2), clf, figure(3), clf

%% Implementations

% FFT length ~ Number of Subcarriers
n_subcarriers = 128;

% Number of OFDM Symbols will be transferred
n_symbols = 500;

% CP
n_cps = 16;
cp = n_subcarriers / n_cps;

% Number of pilot carriers
pilot_distance = 16;
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

% Number of error symbols
n_error_symbols = 0;

% Signal Energy & QAM Normalization Factor
Es = 1;
A = sqrt(3/2/(QAM_Modulator-1)*Es);

% STO & CFO configurations
STOs = [0];
CFOs = [0];

%% Generate Dataset
n_times = 3;

% Generate Pilot value
Tx_pilot = 2*(randn(1,n_pilots) > 0) - 1;


% Calculate Power Delay Profile

% 802.11 Channel Configuration
scale=1e-6;          % nano
Ts=50*scale;         % Sampling time
t_rms=25*scale;     % RMS delay spread
num_ch=1;       % Number of channels
PDP=IEEE802_11_model(t_rms,Ts);

SNRs = [0:5:40];
n_snr = length(SNRs);
all_h = [];

for i=1:(n_snr*length(STOs)*length(CFOs)*n_times)
    for k=1:length(PDP)
        h(:,k) = Ray_model(num_ch).'*sqrt(PDP(k));
    end
    all_h = [all_h ; h];
end

for iSNR=1:length(SNRs)
    SNR = SNRs(iSNR);
    for iSTO=1:length(STOs)
        for iCFO=1:length(CFOs)
            %             BERs = [];
            parfor time=1:n_times
                %% Initialization
                % Get STO & CFO information
                nSTO = STOs(iSTO);
                CFO = CFOs(iCFO);
                
                fprintf('[%.0f/%.0f] SNR = %.0f - STO = %f - CFO = %f\n', time, n_times, SNR, nSTO, CFO);
                
                % temporary even if not having sto/cfo
                sto_mag_dif = [];
                STO_dif = nSTO;
                Est_CFO_CP = CFO;
                Est_CFO_Moose = CFO;
                Est_CFO_Classen = CFO;
                
                
                %% TX-Side
                Tx_data = [];
                Tx_modulated_data = [];
                Tx_CP = [];
                Tx_without_CP = [];
                
                for m=1:n_symbols
                    % Generate input messages
                    msg_input = randi([0 QAM_Modulator-1], 1, n_data);
                    msg_data = qammod(msg_input, QAM_Modulator) * A;
                    Data = merge_pilot_into_msg(msg_data, Tx_pilot);
                    
                    % Convert Frequence-domain to Time-domain
                    time_domain_data = frq_to_time_converter(Data);
                    
                    % Add CP
                    time_domain_data_with_cp = add_cp(time_domain_data, cp);
                    
                    % saving data
                    Tx_data = [Tx_data msg_input];
                    Tx_modulated_data = [Tx_modulated_data msg_data];
                    Tx_without_CP = [Tx_without_CP time_domain_data];
                    Tx_CP = [Tx_CP time_domain_data_with_cp];
                end
                
                %% Generate Channel
                h = all_h((iSNR-1)*n_times+time,:);
                H = fft(h(1,:),n_subcarriers);
                
                % True channel and its time-domain length
                
                H_power_dB = 10*log10(abs(H.*conj(H)));
                channel_length = length(h);
                % Channel Path by Convolution
                y_channel = conv(Tx_CP, h);
                
                if (nSTO ~= 0)
                    y_channel = add_STO(y_channel, -nSTO);
                end
                
                if (CFO ~= 0)
                    y_channel = add_CFO(y_channel,CFO,n_subcarriers);
                end
                
                % Add White Gaussian Noise
                Rx_CP = awgn(y_channel, SNR, 'measured');
                
                %% RX-Side
                % Channel Synchronization
                if (nSTO ~= 0)
                    % STO Estimation
                    com_delay = n_ofdm_symbols / 2;
                    %                     [STO_cor,sto_mag_cor]= STO_by_correlation(Rx_CP,n_subcarriers,cp,com_delay);
                    [STO_dif,sto_mag_dif]= STO_by_difference(Rx_CP,n_subcarriers,cp,com_delay);
                    %
                end
                
                if (CFO ~= 0)
                    % CFO Estimation
                    Est_CFO_CP = CFO_CP(Rx_CP,n_subcarriers,cp);
                    Est_CFO_Moose = CFO_Moose(Rx_CP,n_subcarriers);
                    Est_CFO_Classen = CFO_Classen(Rx_CP,n_subcarriers,cp,Tx_pilot);
                end
                
                Rx_without_CP = [];
                Rx_without_CP_frq = [];
                Rx_modulated_data = [];
                H_est = [];
                Rx_data = [];
                start = 1;
                
                % CFO compensation in Time-Domain
                if (CFO ~= 0)
                    Rx_CP = add_CFO(Rx_CP,-Est_CFO_CP,n_subcarriers);
                end
                
                for m=1:n_symbols
                    % Remove CP
                    y = Rx_CP(cp + start : n_ofdm_symbols + start - 1);
                    
                    % Convert to Frequency-domain
                    Y = time_to_frq_converter(y);
                    
                    if (nSTO ~= 0)
                        % STO compensation in Frequency-Domain
                        Y = add_CFO(Y,-STO_dif,n_subcarriers);
                    end
                    
                    % Channel Estimation by MMSE
                    H_MMSE_est = MMSE_CE(Y, Tx_pilot, pilot_carriers, n_subcarriers, pilot_distance, h, SNR);
                    
                    % Calculate Channel Estimation Power
                    % H_est_power_dB = 10 * log10(abs(H_MMSE_est.*conj(H_MMSE_est)));
                    
                    % DFT-Based Channel Estimation
                    %                     h_est = ifft(H_est);
                    %                     h_DFT = h_est(1:channel_length);
                    %                     H_DFT = fft(h_DFT, n_subcarriers);
                    %                     H_DFT_power_dB = 10 * log10(abs(H_DFT.*conj(H_DFT)));
                    
                    
                    % Equalization
                    Y_eq = channel_equalizer(Y, H_MMSE_est);
                    
                    % Extract Data
                    extracted_Data = extract_data(Y_eq, n_pilots);
                    
                    % Get Extracted Message
                    extracted_msg = qamdemod(extracted_Data/A, QAM_Modulator);
                    Rx_without_CP = [Rx_without_CP y];
                    Rx_without_CP_frq = [Rx_without_CP_frq ; Y];
                    H_est = [H_est ; H_MMSE_est];
                    Rx_modulated_data = [Rx_modulated_data extracted_Data];
                    Rx_data = [Rx_data extracted_msg];
                    start = start + n_ofdm_symbols;
                end
                
                % Error Calculation
                err = 0 + sum(de2bi(Rx_data(:))~=de2bi(Tx_data(:)), 'all');
                ber = err / (length(Tx_data)*n_bit_per_symbol);
                %% Result
                
                fprintf('BER: %f\n', ber);
                % fprintf('MMSE : %6.4e\n', (H-H_est)*(H-H_est)');
                % fprintf('MMSE with DFT: %6.4e\n', (H-H_DFT)*(H-H_DFT)');
                
                % %% Plot Before & After Channel
                % figure(1);
                % plot(abs(X_without_cp), 'r-o','Markersize',4,'linewidth',1); hold on; grid on;
                % plot(abs(y_without_cp), 'g-o','Markersize',4,'linewidth',1);
                % xlabel('Time');
                % ylabel('$x|t|$');
                % legend('True Channel','After Channel Env.','SouthEast');  set(gca,'fontsize',9)
                % signal_power = mean(abs(y_channel.^2));
                % sigma2 = signal_power * 10^(-SNR/10);
                % title(['Before & After Propagation over Chanel. RX Signal Power: ' num2str(signal_power) '. Noise Power: ' num2str(sigma2)]);
                %
                % %% Plot Channel Estimation Power
                % method = 'MMSE';
                % figure(2), subplot(211), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
                % plot(H_est_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 n_subcarriers -10 10])
                % title(method); xlabel('Subcarrier Index'); ylabel('Power[dB]');
                % legend('True Channel',method,'SouthEast');  set(gca,'fontsize',9)
                % subplot(212), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
                % plot(H_DFT_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 n_subcarriers -10 10])
                % title([method ' with DFT']); xlabel('Subcarrier Index'); ylabel('Power[dB]');
                % legend('True Channel',[method ' with DFT'],'SouthEast');  set(gca,'fontsize',9)
                %
                % %% Plot recevied value on Interpolation
                % figure(3), subplot(221), plot(Y,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
                % axis('equal'), set(gca,'fontsize',9), hold on,
                % title('Received data before equalization')
                % subplot(222), plot(Y_eq,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
                % axis('equal'), set(gca,'fontsize',9), hold on,
                % title('Received data after equalization');
                
                %% Save Dataset
                Tx_data_bit = de2bi(Tx_data(:));
                
                dir = ['data_3/' num2str(SNR)];
                if ~exist(dir, 'dir')
                    mkdir(dir)
                end
                
                dir = [dir '/' num2str(nSTO)];
                if ~exist(dir, 'dir')
                    mkdir(dir)
                end
                
                dir = [dir '/' num2str(CFO)];
                if ~exist(dir, 'dir')
                    mkdir(dir)
                end
                
                savetofile([dir '/' num2str(time) '.mat'],...
                    ber,...
                    n_bit_per_symbol,...
                    n_subcarriers,...
                    pilot_carriers,...
                    pilot_distance,...
                    Tx_pilot,...
                    Tx_data_bit,...
                    Tx_data,...
                    Tx_modulated_data,...
                    Tx_without_CP,...
                    Tx_CP,...
                    H,...
                    h,...
                    sto_mag_dif,...
                    STO_dif,...
                    Est_CFO_CP,...
                    Est_CFO_Moose,...
                    Est_CFO_Classen,...
                    H_est,...
                    Rx_CP,...
                    Rx_without_CP,...
                    Rx_without_CP_frq,...
                    Rx_modulated_data,...
                    Rx_data);
                %                 BERs = [BERs ber];
            end
            %             fprintf('======> Mean of BERs: %f\n',mean(BERs));
        end
    end
end