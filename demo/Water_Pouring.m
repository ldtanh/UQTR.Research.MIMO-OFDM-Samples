function [Gamma]=Water_Pouring(Lamda,SNR,nT)

%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

Gamma = zeros(1,length(Lamda));
r=length(Lamda);  index=[1:r];
index_temp=index;
p=1;
while p<r
    irp=[1:r-p+1].'; 
    temp= sum(1./Lamda(index_temp(irp)));
    mu = nT/(r-p+1.)*(1+1/SNR*temp);
    Gamma(index_temp(irp))=mu-nT./(SNR*Lamda(index_temp(irp)));
    if min(Gamma(index_temp))<0
      i=find(Gamma==min(Gamma));  ii=find(index_temp==i);
      index_temp2=[index_temp([1:ii-1]) index_temp([ii+1:end])];
      clear index_temp;
      index_temp=index_temp2;
      p=p+1;
      clear Gamma;
     else
      p=r;
    end
end
Gamma_t=zeros(1,length(Lamda));
Gamma_t(index_temp)=Gamma(index_temp);
%clear Gamma;
Gamma=Gamma_t;