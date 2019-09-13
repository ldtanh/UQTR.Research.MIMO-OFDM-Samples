%channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation|��������/������ֵ��LS/DFT�ŵ�����

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/
rng('shuffle');
clear all; close all; figure(1), clf, figure(2), clf
Nfft=128;  Ng=Nfft/8;  Nofdm=Nfft+Ng;  Nsym=1;
Nps=4; Np=Nfft/Nps; Nd=Nfft-Np; % Pilot spacing, Numbers of pilots and data per OFDM symbol|��Ƶ�����ÿ��OFDM���ŵĵ�Ƶ����ÿ��OFDM���ŵ���������
Nbps=4; M=2^Nbps; % Number of bits per (modulated) symbol|ÿ�����Ʒ��ŵı�����
% mod_object = modem.qammod('M',M, 'SymbolOrder','gray');
% demod_object = modem.qamdemod('M',M, 'SymbolOrder','gray');
Es=1; A=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor|�ź�������QAM��һ������
%fs = 10e6;  ts = 1/fs;  % Sampling frequency and Sampling period
SNRs = [30];  sq2=sqrt(2);
for i=1:length(SNRs)
   SNR = SNRs(i); 
   MSE = zeros(1,6); nose = 0; % ����nose����ͳ�ƴ����������Number_of_symbol_errors
   for nsym=1:Nsym
      Xp = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation|���ɵ�Ƶ���У�-1��+1���������
      %Data = ((2*(randn(1,Nd)>0)-1) + j*(2*(randn(1,Nd)>0)-1))/sq2; % QPSK modulation
      msgint=randi([0,M-1],1,Nfft-Np);    % bit generation|���ɱ���
      % Data = modulate(mod_object,msgint)*A;
      Data = qammod(msgint, M)*A;
      %Data = modulate(mod_object, msgint); Data = modnorm(Data,'avpow',1)*Data;   % normalization
      ip = 0;    pilot_loc = [];
      for k=1:Nfft % ��Ƶ����ض�λ�ü��뵼Ƶ������
         if mod(k,Nps)==1
           X(k) = Xp(floor(k/Nps)+1); pilot_loc = [pilot_loc k]; ip = ip+1;
          else        X(k) = Data(k-ip); % ipָʾ�˵�ǰOFDM�������Ѿ�����ĵ�Ƶ������
         end
      end
      x = ifft(X,Nfft);                            % IFFT
      xt = [x(Nfft-Ng+1:Nfft) x];                  % Add CP|��ѭ��ǰ׺
      h = [(randn+j*randn) (randn+j*randn)/2];     % generates a (2-tap) channel|����һ��2��ͷ�ŵ�
      H = fft(h,Nfft); channel_length = length(h); % True channel and its time-domain length|ʵ���ŵ������ĳ���
      H_power_dB = 10*log10(abs(H.*conj(H)));      % True channel power in dB|ʵ���ŵ��Ĺ���[dB]
      y_channel = conv(xt, h);                     % Channel path (convolution)|�ŵ�·���������
      sig_pow = mean(y_channel.*conj(y_channel));
      %y_aw(1,1:Nofdm) = y(1,1:Nofdm) + ...
      %   sqrt((10.^(-SNR/10))*sig_pow/2)*(randn(1,Nofdm)+j*randn(1,Nofdm)); % Add noise(AWGN)
      yt = awgn(y_channel,SNR,'measured');  
      y = yt(Ng+1:Nofdm);                   % Remove CP|ȥCP
      Y = fft(y);                           % FFT
      for m=1:3
         if m==1, H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear'); method='LS-linear'; % LS estimation with linear interpolation
          elseif m==2, H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline'); method='LS-spline'; % LS estimation with spline interpolation
          else  H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR); method='MMSE'; % MMSE estimation
         end
         H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
         h_est = ifft(H_est); h_DFT = h_est(1:channel_length); 
         H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
         H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
         if nsym==1
           figure(1), subplot(319+2*m), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
           plot(H_est_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 32 -6 10])
           title(method); xlabel('Subcarrier Index'); ylabel('Power[dB]');
           legend('True Channel',method,'SouthEast');  set(gca,'fontsize',9)
           subplot(320+2*m), plot(H_power_dB,'b','linewidth',1); grid on; hold on;
           plot(H_DFT_power_dB,'r:+','Markersize',4,'linewidth',1); axis([0 32 -6 10])
           title([method ' with DFT']); xlabel('Subcarrier Index'); ylabel('Power[dB]');
           legend('True Channel',[method ' with DFT'],'SouthEast');  set(gca,'fontsize',9)
         end
         MSE(m) = MSE(m) + (H-H_est)*(H-H_est)';
         MSE(m+3) = MSE(m+3) + (H-H_DFT)*(H-H_DFT)';
      end
      Y_eq = Y./H_est;
      if nsym>=Nsym-10
        figure(2), subplot(221), plot(Y,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
        axis('equal'), set(gca,'fontsize',9), hold on,
        subplot(222), plot(Y_eq,'.','Markersize',5), axis([-1.5 1.5 -1.5 1.5])
        axis('equal'), set(gca,'fontsize',9), hold on,       
      end
      ip = 0;
      for k=1:Nfft
         if mod(k,Nps)==1, ip=ip+1;  else  Data_extracted(k-ip)=Y_eq(k);  end
      end
      % msg_detected = demodulate(demod_object,Data_extracted/A);
      msg_detected = qamdemod(Data_extracted/A, M);
      nose = nose + sum(msg_detected~=msgint);
   end   
   MSEs(i,:) = MSE/(Nfft*Nsym);
end   
Number_of_symbol_errors=nose
figure(3), clf, semilogy(SNRs',MSEs(:,1),'-x', SNRs',MSEs(:,3),'-o')
legend('LS-linear','MMSE')
fprintf('MSE of LS-linear/LS-spline/MMSE Channel Estimation = %6.4e/%6.4e/%6.4e\n',MSEs(end,1:3));
fprintf('MSE of LS-linear/LS-spline/MMSE Channel Estimation with DFT = %6.4e/%6.4e/%6.4e\n',MSEs(end,4:6));
