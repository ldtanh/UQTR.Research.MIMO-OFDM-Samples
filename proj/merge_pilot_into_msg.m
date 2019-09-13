function X = merge_pilot_into_msg(data, pilot)
    % Pre-allocation
    n_pilot = length(pilot);
    n_subcarriers = length(data) + n_pilot;
    X = zeros(1,n_subcarriers);
    passed_pilot = 0;
    for k=1:n_subcarriers
        if mod(k,fix(n_subcarriers / n_pilot))==1 
            passed_pilot = passed_pilot + 1;
            X(k) = pilot(passed_pilot);
        else
            X(k) = data(k-passed_pilot);
        end
    end