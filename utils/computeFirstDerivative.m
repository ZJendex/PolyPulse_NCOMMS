function first_derivative = computeFirstDerivative(signal, h)
    first_derivative = zeros(size(signal));
    for i = 4:length(signal)-3
        first_derivative(i) = (signal(i+1) - signal(i-1)) / (2*h);
    end
end
