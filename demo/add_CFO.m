function y_CFO=add_CFO(y,CFO,Nfft)
% To add an arbitrary frequency offset|ʩ��CFO
% Input: y    = Time-domain received signal|ʱ������ź�
%        dCFO = FFO (fractional CFO) + IFO (integral CFO)
%        Nfft = FFT size;|FFT��С

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

nn=0:length(y)-1; 
y_CFO = y.*exp(j*2*pi*CFO*nn/Nfft); % ��5.3��133ҳ��
% plot(real(y_CFO),imag(y_CFO),'.')