function idft = frq_to_time_converter(ofdm_data)
    idft = ifft(ofdm_data, length(ofdm_data));