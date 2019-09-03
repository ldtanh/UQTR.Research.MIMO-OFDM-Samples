% plot_IEEE80211_model.m

%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

clear, clf
scale=1e-9;          % nano
Ts=50*scale;         % Sampling time
t_rms=25*scale;     % RMS delay spread
num_ch=10000;       % Number of channels
N=128;                % FFT sizek
PDP=IEEE802_11_model(t_rms,Ts);  
for k=1:length(PDP)
    h(:,k) = Ray_model(num_ch).'*sqrt(PDP(k));
    avg_pow_h(k)= mean(h(:,k).*conj(h(:,k)));
end
H=fft(h(1,:),N);
subplot(221)
stem([0:length(PDP)-1],PDP,'ko'), hold on,
stem([0:length(PDP)-1],avg_pow_h,'k.');
xlabel('channel tap index, p'), ylabel('Average Channel Power[linear]');
title('IEEE 802.11 Model, \sigma_\tau=25ns, T_S=50ns');
legend('Ideal','Simulation');  axis([-1 7 0 1]);
subplot(222)
plot([-N/2+1:N/2]/N/Ts/10^6,10*log10(H.*conj(H)),'k-');
xlabel('Frequency[MHz]'), ylabel('Channel power[dB]')
title('Frequency response, \sigma_\tau=25ns, T_S=50ns');