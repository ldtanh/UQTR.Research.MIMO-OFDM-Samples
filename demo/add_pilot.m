function xp=add_pilot(x,Nfft,Nps)
% CAZAC (Constant Amplitude Zero AutoCorrelation) sequence --> pilot
% CAZAC(�㶨�����������)���� --> ��Ƶ
% Nps : Pilot spacing|��Ƶ���

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

if nargin<3, 
    Nps=4; 
end
Np=Nfft/Nps; % Number of pilots|��Ƶ��
xp=x; % an OFDM signal including pilot signal|׼��������Ƶ�ź����ڵ�OFDM�ź�
for k=1:Np
   xp((k-1)*Nps+1)= exp(j*pi*(k-1)^2/Np);  % Pilot boosting with an even Np|ʽ7.17��185ҳ��
   %xp((k-1)*Nps+1)= exp(j*pi*(k-1)*k/Np);  % Pilot boosting with an odd Np
end
% CAZAC (Constant Amplitude Zero AutoCorrelation) sequence --> pilot