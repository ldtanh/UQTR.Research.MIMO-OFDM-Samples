function y=remove_GI(Ng,Lsym,NgType,ofdmSym)

% MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
% 2010 John Wiley & Sons (Asia) Pte Ltd

% http://www.wiley.com//legacy/wileychi/cho/

if Ng~=0
  if NgType==1  % cyclic prefix|ѭ��ǰ׺
    y=ofdmSym(Ng+1:Lsym);
   elseif NgType==2 % cyclic suffix|ѭ����׺
    y=ofdmSym(1:Lsym-Ng)+[ofdmSym(Lsym-Ng+1:Lsym) zeros(1,Lsym-2*Ng)];
  end
 else
  y=ofdmSym;
end
