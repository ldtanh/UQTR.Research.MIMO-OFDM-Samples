function CFO_est=CFO_CP(y,Nfft,Ng)
% Time-domain CFO estimation based on CP (Cyclic Prefix)
% ����CP��ʱ��CFO����

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

nn=1:Ng; 
CFO_est = angle(y(nn+Nfft)*y(nn)')/(2*pi);  % Eq.(5.23){ԭ������Ϊ5.27�������ϲ���}|ʽ5.23��144ҳ��
% MSE = MSE + (CFO_est-CFO)*(CFO_est-CFO);