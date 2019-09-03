function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
%function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,ts,SNR)
% MMSE channel estimation function|MMSE�ŵ�����
% Inputs:
%       Y         = Frequency-domain received signal|Ƶ������ź�
%       Xp        = Pilot signal|��Ƶ�ź�
%       pilot_loc = Pilot location|��Ƶλ��
%       Nfft      = FFT size|FFT��С
%       Nps       = Pilot spacing|��Ƶ���
%       h         = Channel impulse response|�ŵ�������Ӧ
%       ts        = Sampling time|����ʱ��
%       SNR       = Signal-to-Noise Ratio[dB]|�����dB
% output:
%      H_MMSE     = MMSE channel estimate|MMSE�ŵ�����

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

%H = fft(h,N);
snr = 10^(SNR*0.1);
Np=Nfft/Nps; % ��Ƶ����
k=1:Np; 
H_tilde = Y(1,pilot_loc(k))./Xp(k);  % LS estimate|LS�ŵ�����
k=0:length(h)-1; %k_ts = k*ts; 
hh = h*h'; 
tmp = h.*conj(h).*k; %tmp = h.*conj(h).*k_ts;
r = sum(tmp)/hh;    r2 = tmp*k.'/hh; %r2 = tmp*k_ts.'/hh;
tau_rms = sqrt(r2-r^2);     % rms delay|rmsʱ�ӣ���14ҳ��ʽ1.20
df = 1/Nfft;  %1/(ts*Nfft);
j2pi_tau_df = j*2*pi*tau_rms*df;
K1 = repmat([0:Nfft-1].',1,Np); % K1�Ĵ�СΪNfft*Np��ÿһ�о�Ϊ0:Nfft-1
K2 = repmat([0:Np-1],Nfft,1); % K2�Ĵ�СΪNfft*Np��ÿһ�о�Ϊ0:Np-1
rf = 1./(1+j2pi_tau_df*(K1-K2*Nps)); % ʽ6.16
K3 = repmat([0:Np-1].',1,Np); % K3�Ĵ�СΪNp*Np��ÿһ�о�Ϊ0:Np-1
K4 = repmat([0:Np-1],Np,1); % K4�Ĵ�СΪNp*Np��ÿһ�о�Ϊ0:Np-1
rf2 = 1./(1+j2pi_tau_df*Nps*(K3-K4)); % ʽ6.17
Rhp = rf;
Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr;
H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.');  % MMSE channel estimate|MMSE�ŵ����ƣ�ʽ6.14