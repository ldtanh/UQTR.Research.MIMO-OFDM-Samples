clear,clc;

%% Configuration
n_bit_per_symbol = 4;
n_continous_subcarriers = 8;
inputSize = [1 n_continous_subcarriers*2 1];
outputSize = [1 n_continous_subcarriers*n_bit_per_symbol];

data_dir = '/home/anhldt/mimo-ofdm/proj/data_2/';

sto = 0;
cfo = 0;
times = 10;

Mean_BL_BERs = [];
Mean_NN_BERs = [];

%% Load Neural Networks
load('net/SD_5_0_0_model.mat');

%% Calculation
for snr=-10:5:30
    fprintf('------  SNR = %.0f  ----------\n',snr);
    Baseline_BERs = [];
    NN_BERs = [];
    for time=1:times
        fprintf('[%.0f/%.0f] STO = %.0f - CFO = %.2f\n', time, times, sto, cfo);
        %% Load dataset
        dataset_fn = [data_dir num2str(snr) '/' num2str(sto) '/' num2str(cfo) '/' num2str(time) '.mat'];
        load(dataset_fn);
        
        %% Generate Data
        n_total = length(Rx_modulated_data);
        
        x_real = real(Rx_modulated_data);
        x_imag = imag(Rx_modulated_data);
        
        n_dataset = n_total * 2 / n_continous_subcarriers;
        
        X = zeros(1,n_continous_subcarriers*2,1,n_dataset);
        Y = zeros(n_dataset,n_continous_subcarriers*n_bit_per_symbol);
        
        for i=1:n_dataset-1
            for k=1:n_continous_subcarriers
                idx = (i-1) * (n_continous_subcarriers / 2);
                
                X(1,(k-1)*2+1,1,i) = x_real(1,idx + k);
                X(1,(k-1)*2+2,1,i) = x_imag(1,idx + k);
                
                for j=1:n_bit_per_symbol
                    Y(i,(k-1)*n_bit_per_symbol + j) = Tx_data_bit(idx + k,j);
                end
            end
        end
        
        ratio = 0.25;
        n_train = round(n_dataset * ratio);
        n_test = n_dataset - n_train;
        
        X_train = X(:,:,:,1:n_train);
        Y_train = Y(1:n_train,:);
        
        X_test = X(:,:,:,n_train+1:n_dataset);
        Y_test = Y(n_train+1:n_dataset,:);
        
        % Cyclic for Testset
        idx = (n_dataset-1) * (n_continous_subcarriers / 2);
        for k=1:n_continous_subcarriers/2
            X_test(1,(k-1)*2+1,1,n_test) = x_real(1,idx + k);
            X_test(1,(k-1)*2+2,1,n_test) = x_imag(1,idx + k);
            
            for j=1:n_bit_per_symbol
                Y_test(n_test,(k-1)*n_bit_per_symbol + j) = Tx_data_bit(idx + k,j);
            end
        end
        
        idx = n_train * (n_continous_subcarriers / 2);
        for k=(n_continous_subcarriers/2+1):n_continous_subcarriers
            X_test(1,(k-1)*2+1,1,n_test) = x_real(1,idx+k-n_continous_subcarriers/2);
            X_test(1,(k-1)*2+2,1,n_test) = x_imag(1,idx+k-n_continous_subcarriers/2);
            
            for j=1:n_bit_per_symbol
                Y_test(n_test,(k-1)*n_bit_per_symbol + j) = Tx_data_bit(idx+k-n_continous_subcarriers/2,j);
            end
        end
        
        %% Validation
        Y_pred = zeros(1,n_test * n_bit_per_symbol * n_continous_subcarriers / 2);
        
        for i=1:(n_test-1)
            y_pred = net.predict(X_test(:,:,:,i));
            idx = (i-1)*length(y_pred)/2;
            Y_pred(idx + 1:idx + length(y_pred)) = Y_pred(idx + 1:idx + length(y_pred)) + y_pred(1,:);
        end
        
        % Cyclic
        y_pred = net.predict(X_test(:,:,:,n_test));
        idx = length(Y_pred) - 0.5 * length(y_pred);
        Y_pred(idx + 1 : length(Y_pred)) = Y_pred(idx + 1 : length(Y_pred)) + y_pred(1,1:0.5*length(y_pred));
        Y_pred(1 : 0.5 * length(y_pred)) = Y_pred(1 : 0.5 * length(y_pred)) + y_pred(1,0.5*length(y_pred)+1:length(y_pred));
        
        Y_pred = Y_pred / 2.0;
        Y_pred(Y_pred>0.5) = 1;
        Y_pred(Y_pred<=0.5) = 0;
        
        l = n_train*4 + 1;
        r = length(Tx_data);
        
        err_p = 0;
        dif = 0;
        
        tx_bit = de2bi(Tx_data(l:r));
        rx_bit = de2bi(Rx_data(l:r));
        
        for i=1:(r-l+1)
            err_p = err_p + sum(Y_pred(1,(i-1)*n_bit_per_symbol+1:i*n_bit_per_symbol)~=tx_bit(i,:),'all');
            dif = dif + sum(Y_pred(1,(i-1)*n_bit_per_symbol+1:i*n_bit_per_symbol)~=rx_bit(i,:),'all');
        end
        
        ber_p = err_p / length(Y_pred);
        fprintf('BER by DL-method: %f\n', ber_p);
        
        err_e = sum(tx_bit~=rx_bit, 'all');
        ber_e = err_e / ((r - l + 1)*n_bit_per_symbol);
        fprintf('BER by Baseline method: %f\n', ber_e);
        
        Baseline_BERs = [Baseline_BERs ber_e];
        NN_BERs = [NN_BERs ber_p];
    end
    
    fprintf('===========================\n');
    mean_BL = mean(Baseline_BERs);
    mean_NN = mean(NN_BERs);
    
    fprintf('Deep-Learning: %f\n',mean_NN);
    fprintf('Baseline: %f\n',mean_BL);
    
    Mean_BL_BERs = [Mean_BL_BERs mean_BL];
    Mean_NN_BERs = [Mean_NN_BERs mean_NN];
end

figure(2);
semilogy([-10:5:30], Mean_BL_BERs, 'r','Markersize',4,'linewidth',1); hold on; grid on;
semilogy([-10:5:30], Mean_NN_BERs, 'b--o','Markersize',4,'linewidth',1); axis([-10 30 1e-1 1]);
title('BER per SNR'); xlabel('SNR'); ylabel('BER');
legend('Baseline','Deep Learning');  set(gca,'fontsize',9);