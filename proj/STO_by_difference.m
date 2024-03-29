function [STO_est,Mag]=STO_by_difference(y,Nfft,Ng,com_delay)
% STO estimation by minimizing the difference between CP and rear part of OFDM symbol
% ͨ����С��CP��OFDM���ź󲿵Ĳ�ֵ��ʵ��STO�Ĺ���
% estimates STO by minimizing the difference between CP (cyclic prefix) 
%     and rear part of OFDM symbol
% Input:  y          = Received OFDM signal including CP|����CP��OFDM�����ź�
%          Ng         = Number of samples in CP (Guard Interval)|CP/GI�ڵĲ�����
%          com_delay = Common delay|����ʱ��
% Output: STO_est   = STO estimate|STO����
%           Mag        = Correlation function trajectory varying with time
%                   ��غ�����ʱ��켣

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

N_ofdm=Nfft+Ng; minimum=100; STO_est=0;
if nargin<4, com_delay = N_ofdm/2; end
for n=1:N_ofdm
   nn = n+com_delay+[0:Ng-1]; 
   %tmp1 = y(nn)-y(nn+Nfft); tmp = sum(abs(tmp1)); % Eq. (5.11) works for CFO=0
   %tmp1 = y(nn)-y(nn+Nfft); tmp = tmp1*tmp1'; % works for CFO=0
   %tmp = y(nn)*y(nn)'+y(nn+Nfft)*y(nn+Nfft)'-2*abs(y(nn))*abs(y(nn+Nfft))';
   tmp0 = abs(y(nn))-abs(y(nn+Nfft));
   Mag(n) = tmp0*tmp0'; % Eq.(5.11) is strong against CFO{ԭʼע��Ϊ5.12�������ϲ���}|ʽ5.11����5.10����CFO��
   %tmp0= abs(y(nn)-conj(y(nn+Nfft))); tmp= tmp0*tmp0.';
   %discrepancy = tmp - tmp2
   %Mag(n) = Mag(n) + tmp;
   if (Mag(n)<minimum)
     minimum=Mag(n);  STO_est = N_ofdm-com_delay -(n-1); 
   end
end
%if (STO_est>=Ng/2), STO_est= Ng/2-1;
% elseif (STO_est<-Ng/2), STO_est= -Ng/2;
%end