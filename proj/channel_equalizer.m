function result = channel_equalizer(received_signal, channel_estimation)
    result = received_signal ./ channel_estimation;
    