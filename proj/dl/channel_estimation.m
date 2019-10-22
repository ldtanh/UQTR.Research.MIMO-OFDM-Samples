clear,clc;

%% Configuration
SNRs = [20:5:40];
sto = 0;
cfo = 0;

n_symbols = 5000;
n_subcarriers = 128;
pilot_distance = 16;
n_pilot = n_subcarriers / pilot_distance;

n_pilot = n_subcarriers / pilot_distance;
inputSize = [1 n_subcarriers*2 1];
outputSize = [1 n_subcarriers*2 1];

data_dir = '/home/anhldt/mimo-ofdm/proj/data_3/';

%% Generate Data
times = 100;
n_dataset = n_symbols * times * length(SNRs);

X = zeros(1,n_pilot*2,1,n_dataset);
Y = zeros(n_dataset,n_subcarriers*2);

idx = 0;
for i_snr=1:length(SNRs)
    snr = SNRs(i_snr);
    for time=1:times
        dataset_fn = [data_dir num2str(snr) '/' num2str(sto) '/' num2str(cfo) '/' num2str(time) '.mat'];
        fprintf('[%.0f][%.0f/%.0f] Loading dataset [%s]\n', snr, time, times, dataset_fn);
        load(dataset_fn);
        for k=1:length(H_est)
            idx = idx + 1;
            
%             % LS Est.
%             pilot_idx = 1:length(Tx_pilot);
%             ls_est = Rx_without_CP_frq(k,pilot_carriers(pilot_idx))./Tx_pilot(pilot_idx);
%             H_LS = interpolate(ls_est,pilot_carriers,n_subcarriers,'linear');
            
            h_ls_real = real(H_est(k,:));
            h_ls_imag = imag(H_est(k,:));
            
            X(1,1:n_subcarriers,1,idx) = h_ls_real(:);
            X(1,n_subcarriers+1:n_subcarriers*2,1,idx) = h_ls_imag(:);
            
            Y(idx,1:n_subcarriers) = real(H(1,:));
            Y(idx,n_subcarriers+1:n_subcarriers*2) = imag(H(1,:));
        end
    end
end

ratio = 1;
n_train = round(n_dataset * ratio);
n_test = n_dataset - n_train;

X_train = X(:,:,:,1:n_train);
Y_train = Y(1:n_train,:);

X_test = X(:,:,:,n_train+1:n_dataset);
Y_test = Y(n_train+1:n_dataset,:);

%% Neural Networks Implementation
layers = [
    imageInputLayer(inputSize,'Normalization','none')
    
    fullyConnectedLayer(n_subcarriers*2)
    fullyConnectedLayer(n_subcarriers*2)
        
    regressionLayer
    ];

miniBatchSize = 256;
validationFrequency = 200;

%     'ValidationFrequency',validationFrequency, ...
%     'ValidationData', {X_test Y_test},...

options = trainingOptions('adam', ...
    'MaxEpochs',500,...
    'InitialLearnRate',1e-3, ...
    'MiniBatchSize', miniBatchSize,...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.5,...
    'LearnRateDropPeriod',250,...
    'Shuffle', 'every-epoch');

net = trainNetwork(X_train, Y_train, layers, options);

%% Save Net
if ~exist('net/', 'dir')
    mkdir('net')
end

save(['net/CE_adv_model_8.mat'], 'layers', 'options', 'net');