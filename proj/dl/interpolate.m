function H_interpolated = interpolate(H_est,pilot_loc,Nfft,method)
% Input:        H_est    = Channel estimate using pilot sequence|ʹ�õ�Ƶ���е��ŵ�����
%           pilot_loc    = location of pilot sequence|��Ƶ���е�λ��
%                Nfft    = FFT size|FFT��С
%              method    = 'linear'/'spline'|��ֵ����������/����������ֵ
% Output: H_interpolated = interpolated channel|��ֵ����ŵ�

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

if pilot_loc(1)>1 % �����һ����ŵ�
  slope = (H_est(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
  H_est = [H_est(1)-slope*(pilot_loc(1)-1)  H_est]; pilot_loc = [1 pilot_loc];
end
if pilot_loc(end)<Nfft % �������һ����ŵ�
  slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));  
  H_est = [H_est  H_est(end)+slope*(Nfft-pilot_loc(end))]; pilot_loc = [pilot_loc Nfft];
end
if lower(method(1))=='l', H_interpolated = interp1(pilot_loc,H_est,[1:Nfft]);   
 else      H_interpolated = interp1(pilot_loc,H_est,[1:Nfft],'spline');
end  