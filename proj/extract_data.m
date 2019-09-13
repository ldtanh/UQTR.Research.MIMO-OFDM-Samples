function data = extract_data(equalized_data, n_pilots)
    n_subcarriers = length(equalized_data);
    n_data = n_subcarriers - n_pilots;
    data = zeros(1,n_data);
    passed_pilot = 0;
    for k=1:n_subcarriers
        if mod(k,fix(n_subcarriers / n_pilots)) == 1
            passed_pilot = passed_pilot + 1;
        else
            data(k - passed_pilot) = equalized_data(k);
        end
    end