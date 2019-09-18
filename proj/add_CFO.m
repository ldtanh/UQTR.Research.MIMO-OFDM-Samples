function y_CFO=add_CFO(y,CFO,Nfft)
    nn=0:length(y)-1; 
    y_CFO = y.*exp(1i*2*pi*CFO*nn/Nfft);