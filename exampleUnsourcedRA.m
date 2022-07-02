function exampleUnsourcedRA
%EXAMPLEUNSOURCEDRA Calculate error rates of unsourced RA with SPARCs as a function of Eb/N0 
%   See Fig.8b in [1], Calculations may take ~2min per data point per iteration
    L       = 100;
    S       = 32;
    J       = 12;
    b       = 96;
    M       = 32;
    K_a     = 300;
    iter    = 10;
    
    % Specify fading
    fading.type         = 'no_fading';
    fading.lower_limit  = 0;
    fading.upper_limit  = 0;
    
    % Data_profile speficies the number of data bits per sections
    data_profile     = ones(S,1)*(3);
    data_profile(1)  = J;
    data_profile(32) = 0;
    data_profile(31) = 0;
    data_profile(30) = 0;
    
    
    assert(abs(sum(data_profile)-b)<eps);
    rate = sum(data_profile)/(S*L);
    
    Eb_N0_dB    = -5:1:5;
    P_dB        = Eb_N0_dB + 10*log10(rate);
    
    % Specify power allocation
    PA.nSections          = S;    
    PA.decay              = 15.5;     
    PA.cutoff             = 0.7; 
    PA.sigmaw2            = 1.0;
    PA.method             = 'exponential_flat_tail';
    
    for i=1:length(Eb_N0_dB)
        PA.aver_power   = 10.^(P_dB(i)/10);
        P               = makePowerVector(PA);
        [p_md,p_fa]     = unsourcedSPARC(L, S, J, K_a, M, data_profile, P, iter, fading);
        PE(i)           = p_md + p_fa;
    end
    
    semilogy(Eb_N0_dB,PE);
end

