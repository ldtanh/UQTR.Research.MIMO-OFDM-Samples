function data_with_cp = add_cp(data, cp)
    n_subcarriers = length(data);
    data_with_cp = [data(n_subcarriers - cp + 1:n_subcarriers) data];