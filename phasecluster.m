function pc = phasecluster(signal)

    pc = squeeze(mean(exp(1i*signal)));

end