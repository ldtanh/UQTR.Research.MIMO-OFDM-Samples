function y=add_CP(x,Ncp)
% Add cyclic prefix|��ѭ��ǰ׺(CP)

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

y = [x(:,end-Ncp+1:end) x];