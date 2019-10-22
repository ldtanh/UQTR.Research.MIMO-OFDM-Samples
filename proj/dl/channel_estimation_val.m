clear,clc;

%% Configuration
n_subcarriers = 128;
pilot_distance = 16;
n_pilot = n_subcarriers / pilot_distance;

data_dir = '/home/anhldt/mimo-ofdm/proj/data_2/';

sto = 0;
cfo = 0;
times = 10;

Mean_BL_MSEs = [];
Mean_NN_MSEs = [];

%% Load Neural Networks
load('net/CE_adv_model_7.mat');

%% Calculation
for snr=0:5:30
    fprintf('------  SNR = %.0f  ----------\n',snr);    
    %% Generate Data   
    Baseline_MSEs = 0;
    NN_MSEs = 0;
    
    for time=1:times
        dataset_fn = [data_dir num2str(snr) '/' num2str(sto) '/' num2str(cfo) '/' num2str(time) '.mat'];
        fprintf('[%.0f/%.0f] Loading dataset [%s] and validating...\n', time, times, dataset_fn);
        load(dataset_fn);
        for k=1:length(H_est)
            X = zeros(1,n_subcarriers*2,1);
            
            % LS Est.
%             pilot_idx = 1:length(Tx_pilot);
%             ls_est = Rx_without_CP_frq(k,pilot_carriers(pilot_idx))./Tx_pilot(pilot_idx);
%             H_LS = interpolate(ls_est,pilot_carriers,n_subcarriers,'linear');
            
            h_ls_real = real(H_est(k,:));
            h_ls_imag = imag(H_est(k,:));
            
            X(1,1:n_subcarriers,1) = h_ls_real(:);
            X(1,n_subcarriers+1:n_subcarriers*2,1) = h_ls_imag(:);
            
            y_pred = net.predict(X);
            H_est_NN = complex(y_pred(1:n_subcarriers),y_pred(n_subcarriers+1:2*n_subcarriers));
            
            Baseline_MSEs = Baseline_MSEs + ((H-H_est(k,:))*(H-H_est(k,:))' / n_subcarriers);
            NN_MSEs = NN_MSEs + ((H-H_est_NN)*(H-H_est_NN)' / n_subcarriers);
        end
    end
    
    %% Validation
    
    fprintf('===========================\n');
    mean_BL = Baseline_MSEs / (length(H_est) * times * n_subcarriers);
    mean_NN = NN_MSEs / (length(H_est) * times * n_subcarriers);
    
    mse_BL = 10*log10(mean_BL);
    mse_NN = 10*log10(mean_NN);
    
    fprintf('Baseline: %f (dB)\n',mse_BL);
    fprintf('Deep-Learning: %f (dB)\n',mse_NN);
    
    Mean_BL_MSEs = [Mean_BL_MSEs mse_BL];
    Mean_NN_MSEs = [Mean_NN_MSEs mse_NN];
end

figure(1);
% plot([0:5:30], Mean_BL_MSEs, 'r','Markersize',4,'linewidth',1); hold on; grid on;
plot([0:5:30], Mean_NN_MSEs, 'm--o','Markersize',4,'linewidth',1);
% title('MSE per SNR'); xlabel('SNR'); ylabel('MSE');
% legend('Baseline','Deep Learning');  set(gca,'fontsize',9);