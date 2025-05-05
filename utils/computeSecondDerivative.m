function second_derivative = computeSecondDerivative(signal, h)
    second_derivative = zeros(size(signal));
    for i = 4:length(signal)-3
        second_derivative(i) = (4*signal(i) + (signal(i+1) + signal(i-1)) ...
            - 2*(signal(i+2) + signal(i-2)) - (signal(i+3) + signal(i-3))) / (16*h^2);
    end

end
