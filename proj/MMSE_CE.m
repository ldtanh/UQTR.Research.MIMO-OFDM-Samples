function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
    snr = 10^(SNR*0.1);
    Np=Nfft/Nps; % µ¼ÆµÊýÁ¿
    k=1:Np; 
    H_tilde = Y(1,pilot_loc(k))./Xp(k);  % LS estimate|LSÐÅµÀ¹À¼Æ
    k=0:length(h)-1; %k_ts = k*ts; 
    hh = h*h'; 
    tmp = h.*conj(h).*k; %tmp = h.*conj(h).*k_ts;
    r = sum(tmp)/hh;    r2 = tmp*k.'/hh; %r2 = tmp*k_ts.'/hh;
    tau_rms = sqrt(r2-r^2);     % rms delay|rmsÊ±ÑÓ£¬Êé14Ò³£¬Ê½1.20
    df = 1/Nfft;  %1/(ts*Nfft);
    j2pi_tau_df = j*2*pi*tau_rms*df;
    K1 = repmat([0:Nfft-1].',1,Np); % K1µÄ´óÐ¡ÎªNfft*Np£¬Ã¿Ò»ÁÐ¾ùÎª0:Nfft-1
    K2 = repmat([0:Np-1],Nfft,1); % K2µÄ´óÐ¡ÎªNfft*Np£¬Ã¿Ò»ÐÐ¾ùÎª0:Np-1
    rf = 1./(1+j2pi_tau_df*(K1-K2*Nps)); % Ê½6.16
    K3 = repmat([0:Np-1].',1,Np); % K3µÄ´óÐ¡ÎªNp*Np£¬Ã¿Ò»ÁÐ¾ùÎª0:Np-1
    K4 = repmat([0:Np-1],Np,1); % K4µÄ´óÐ¡ÎªNp*Np£¬Ã¿Ò»ÐÐ¾ùÎª0:Np-1
    rf2 = 1./(1+j2pi_tau_df*Nps*(K3-K4)); % Ê½6.17
    Rhp = rf;
    Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr;
    H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.');  % MMSE channel estimate|MMSEÐÅµÀ¹À¼Æ£¬Ê½6.14