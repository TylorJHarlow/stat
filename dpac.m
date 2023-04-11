
% Debiased Phase-Amplitude Coupling
function output = dpac(low, high)
    output = squeeze(abs(mean(abs(high) .* (exp(1i*(angle(low))) - phasecluster(low)))));
end