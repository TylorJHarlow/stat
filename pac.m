function output = pac(low, high)

    output = squeeze(abs(mean(high.*exp(1i*low))));

end