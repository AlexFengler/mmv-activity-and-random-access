function out = exampleActivity
%EXAMPLEACTIVITY Calculate PE over M
% Fig.2 in [1]
    L       = 100;
    N       = 2000;
    K_a     = 300;
    P       = 1;
    iter    = 10;
    
    fading.type = 'no_fading';
    fading.lower_limit = 0;
    fading.upper_limit = 30;
    
    
    M_s = 50:50:400;
    for m = 1:length(M_s)
        M       = M_s(m);
        PE(m)   = activityDetectionPE(L,N,K_a,M,P,iter,'fading',fading);
    end
    semilogy(M_s,PE);
    grid on;
    xlabel('M');
    ylabel('P_{md}');
end

