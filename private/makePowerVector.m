function p = makePowerVector(parameters)
% Return a power vector such that the average power per section is
% parameters.aver_power
    L       = parameters.nSections;

    snr     = parameters.aver_power / parameters.sigmaw2;
    
    switch(parameters.method)
        case 'exponential'           
            I           = 1:L;
            capacity    = log2(1+snr)/2;
            p           = exp(-parameters.decay*capacity*I/L);
            
        case 'exponential_flat_tail'
            k = floor(parameters.cutoff*L);
            I = 1:k;
            capacity = log2(1+snr)/2;
            p           = exp(-parameters.decay*capacity*I/L);
            p(k+1:L)    = p(k);
            
        case 'flat' 
            p = ones(1,L)./L;
            
        otherwise
            fprintf('unknown power allocation method <%s>\n', method);
            exit;
    end
    
    p = L*parameters.aver_power*p/sum(p);
end


