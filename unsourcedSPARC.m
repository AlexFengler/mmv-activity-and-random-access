function [p_md,p_fa] = unsourcedSPARC(L, S, J, K_a, M, data_profile, P, iter, fading)
%unsourcedSPARC Compute the error rate of the MU-SPARC coding scheme
%
%   L:              Blocklength per slot
%   S:              Number of slots
%   J:              log2 of the size of the sections
%   K_a:            active users
%   M:              receive antennas
%   data_profile:   vector of the number of data bits per section
%   P:              input power-per-user
%   iter:           number of simulation runs
%
%   struct fading with the following fields
%   'type' = {'no_fading'(default),'uniform','exp','pathloss','shadowing_pathloss'}
%   'lower_limit'(default = 10),'upper_limit'(default = 30): limits for the uniform case

% No. of iterations for the inner ML decoder
iter_ml         = 10;

N               = 2^J*S;
parity_profile  = J - data_profile;

% Coding Matrix, columns are normalized to 1 (entries scale like 1/sqrt(L))
A       = (1 + 1i)/sqrt(2)*exp(2*pi*1i*rand(L,2^J));
A       = A./sqrt(L);
norms   = sqrt(sum(abs(A.^2)));
assert(all(abs(norms-ones(1,2^J))<1e-10));

% Lookup table for outer code
LT  = createLT(2^max(data_profile));

p_md = 0;
p_fa = 0;

for i=1:iter
    % Vector of LSFCs
    lsf = LS_fading(K_a,fading.type, fading.lower_limit, fading.upper_limit);
    
    % random messages
    msg = generateMessages(data_profile, K_a);

    % C is the matrix of parity checks needed for the decoder
    [symbols, C] = outerEncoder( msg, 2^J, S, data_profile );

    % n times K_a matrix of transmitted messages
    TX  = innerEncoder(symbols, A, P);

    % MU-MIMO MAC-Channel, with Gaussian channel matrix,
    % User messages are attenuated by the LSFCs
    RX = gmacMMV(TX, M, 1, lsf);

    g_hat = innerDecoder(RX, A, S, P, J, iter_ml);
       
    % Hard decision
    supp        = (g_hat > 0.25*min(lsf));
    
    symbol_list = getSymbolList(reshape(supp,2^J,S));
    
    msg_list    = outerDecoder(symbol_list, C, parity_profile, 2^J, LT);

    [fn,fp]     = countErrors(msg_list, msg);
    p_md        = p_md + fn/K_a;
    p_fa        = p_fa + fp/length(msg_list);
end
p_md = p_md/iter;
p_fa = p_fa/iter;

end


%create lookuptable for sums of bits
function LT = createLT(N)
    LT = zeros(N,1);
    for i = 1:N
        LT(i) = mod(sum(dec2bin(i-1)-'0'),2);
    end
end

% Return: LxK matrix of symbols where each symbol has data_profile(l) bits
function msg_matrix = generateMessages(data_profile, K)
    L           = length(data_profile);
    msg_matrix  = zeros(L, K);
    for i = 1:L
        msg_matrix(i,:) = randi([0,2^data_profile(i) - 1], 1, K);
    end
end

% Superpose the columns of A indicated by the indices
% in symbols
function TX = innerEncoder(symbols, A, P)
    [S,K] = size(symbols);
    [L,~] = size(A);
    TX = zeros(S*L,K);
    for s = 1:S
        ind                    = symbols(s,:);
        TX((s-1)*L + 1:s*L,:)  = sqrt(L*P(s))*A(:,ind);
    end
end

% Apply the relaxed ML decoder to each section
function g = innerDecoder(RX, A, S, P, J, iter)
    [n,M] = size(RX);
    [L,N] = size(A);
    B     = 2^J;
    g     = zeros(N,1);
    
    for s=1:S
        x_l = (s-1)*B+1;
        x_r = s*B;
        y_l = (s-1)*L+1;
        y_r = s*L;
        
        Y       = RX(y_l:y_r,:);
        cov_m   = Y*Y'/M;
        tic;
        g(x_l:x_r) = ML_coord_descent_round(cov_m, A, iter, 1, [])./P(s)./L;
        toc;
    end
        
end

% Create list of position indices from a support vector
function symbols = getSymbolList(rpos_matrix)
    L = size(rpos_matrix,2);
    for i = 1:L
        % outer decoder wants symbols from 0:B-1;
        symbols{i} = find(rpos_matrix(:,i)==1)' - 1;
    end
end

% Count the messages which miss in the output list, and the messages
% which were not transmitted but appeared in the output list
function [errs, false_positives] = countErrors(est, sparc_matrix)
    K       = size(est,2);
    errs    = size(sparc_matrix,2);
    
    false_positives = 0;
    for i = 1:K
        codeword = est(:,i);
        [~,indx]=ismember(codeword',sparc_matrix','rows');
        if (indx~=0)
            sparc_matrix(:,indx) = [];
            errs = errs -1;
        else
            false_positives = false_positives + 1;
        end        
    end
end