function lam=ML_coord_descent_round(samp_cov,pilot_mat, nround, noise_var, gamma )

P = size(samp_cov,1);
N = size(pilot_mat,2);
% initialize estimate and covariance matrix
lam     = zeros(N,1);
lam_old = zeros(N,1);
cov_inv = eye(P)./noise_var;


% start optimization
for r=1:nround
    % generate a random permutation
    ind_perm=randperm(N);
    
    % Sort in descending order
    %[~,I] = sort(lam,'descend');
    %ind_perm = I';
    
    % make a nround of update according to the permutation
    for ind=ind_perm
        % the selected pilot vector
        pilot=pilot_mat(:,ind);
        
        p_vec = cov_inv*pilot;
        
        % set the step-size
        num1=real(p_vec'*(samp_cov*p_vec));
        
        num2=real(pilot'*p_vec);
        
        % find the optimal step size
        step=max( (num1 - num2)/ num2^2, -lam(ind));
        
        % box-constraint if LSF are known
        if (~isempty(gamma))
            step = min(step,gamma(ind)-lam(ind));
        end
        
        % update the lam vector
        lam(ind)=lam(ind)+step;
        
        if(step~=0)
            % update covariance inverse estimate
            cov_inv = cov_inv - step*(p_vec*p_vec')./(1 + step*num2);
        end
    end
    
    % Stopping criterion
    if(norm(lam-lam_old,1)<1e-4)
        break
    end
    lam_old = lam;
end


end