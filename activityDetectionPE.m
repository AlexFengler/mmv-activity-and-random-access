function [p_md] = activityDetectionPE(L,N,K_a,M,P,iter,varargin)
%ACTIVITYDETECTIONPE Activity Detection with the relaxed ML algorithm
%   Compute the probability of missed detection in a fading MIMO-MAC when
%   the pilots are choosen at random and the activity is detected by
%   picking the positions of the K_a largest estimated coefficients.
%
%   L:      Pilot dimension
%   N:      total users
%   K_a:    active users
%   M:      receive antennas
%   P:      input power
%   iter:   number of simulation runs
%
%   optional inputs:
%   struct fading with the following fields
%   'type' = {'no_fading'(default),'uniform','exp','pathloss','shadowing_pathloss'}
%   'lower_limit'(default = 10),'upper_limit'(default = 30): limits for the uniform case
%   'iter_ml': No. of iterations for the ML algorithm default=10

    fading.type = 'uniform';
    fading.lower_limit = 10;
    fading.upper_limit = 30;
    
    iter_ml = 10;
    
    % Input parsing
    if nargin > 6
        names   = varargin(1:2:end);
        values  = varargin(2:2:end);
        for i=1:length(names)
            switch names{i}
                case 'iter_ml'
                    iter_ml = values{i};
                case 'fading'
                    fading = values{i};
            end
        end
    end
    
    
    % Create random pilot matrix
    A = (1 + 1i)/sqrt(2)*exp(2*pi*1i*rand(L,N));
    
    % Columns are normalized to L
    norms = sum(abs(A.^2));
    assert(all(abs(norms-L*ones(1,N))<1e-10));
    
    % Large-scale fading
    sampleLSF = @(N) LS_fading(N,fading.type, fading.lower_limit, fading.upper_limit);
    
    md = 0;
    
    for i=1:iter       
        % Create gamma
        gamma = sampleLSF(N);

        % Random activity
        sup = randperm(N,K_a);
        
        % Small-scale fading
        X           = zeros(N,M);
        X(sup,:)    = diag(sqrt(gamma(sup)))*(randn(K_a,M) + 1i*randn(K_a,M))/sqrt(2);
        X           = sqrt(P)*X;
        
        % AWGN
        Z = (randn(L,M) + 1i*randn(L,M))/sqrt(2);
        
        % MMV-MAC
        Y = A*X + Z;
        
        % Detection
        cov_m = Y*Y'/M;
        tic;
        g_hat = ML_coord_descent_round(cov_m, A, iter_ml, 1,[]);
        toc;
        
        % Count errors under picking the K_a largest entries
        md = md + countErrors(g_hat, sup, K_a)./K_a;
    end
    p_md = md/iter;
end

function md = countErrors(g_hat, sup, K_a)
    [~,sup_est] = maxk(g_hat,K_a);
    md      = K_a - numel(intersect(sup_est,sup));
end
