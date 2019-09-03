function H_LS = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt)
% LS channel estimation function|LS�ŵ�����
% Inputs:
%       Y         = Frequency-domain received signal|Ƶ������ź�
%       Xp        = Pilot signal|��Ƶ�ź�
%       pilot_loc = Pilot location|��Ƶλ��
%       N         = FFT size|FFT��С
%       Nps       = Pilot spacing|��Ƶ���
%       int_opt   = 'linear' or 'spline'|��ֵ����������/����������ֵ
% output:
%       H_LS      = LS channel etimate|LS�ŵ�����

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

Np=Nfft/Nps; % ��Ƶ����
k=1:Np; 
LS_est(k) = Y(pilot_loc(k))./Xp(k);  % LS channel estimation|LS�ŵ�����
if  lower(int_opt(1))=='l', % ȷ����ֵ����
    method='linear'; 
else
    method='spline';  
end
H_LS = interpolate(LS_est,pilot_loc,Nfft,method); % Linear/Spline interpolation|����/����������ֵ