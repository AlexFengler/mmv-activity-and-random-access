function RX = gmacMMV(TX, M, sigma2, lsf)
%GMACMMV Transmit multiple sparc messages over the GMAC
%   with M receive antennas
%   TX: message matrix where each column is a message from a user
%   RX: column i is the vector of received symbols of antenna i 
    K   = size(TX,2);
    n   = size(TX,1);

    % GMAC with M receive antennas
    mac_sum = zeros(n,M);
    for m = 1:M
        channel = (randn(1,K) + 1i*randn(1,K))/sqrt(2);
        channel = channel.*sqrt(lsf)'; % Large scale fading
        faded   = bsxfun(@times, TX, channel);
        mac_sum(:,m) = sum(faded,2);
    end

    Z  = sqrt(sigma2)*(randn(n,M) + 1i*randn(n,M))./sqrt(2);
    RX = mac_sum + Z;
end