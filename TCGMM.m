function [C, pi, mu, sigma, A_kl,A_kl_init] = TCGMM(X, K, tol)
maxIter = 400;
iter = 1;
% Number of clusters
L = K; 
% Initialise GMM parameters
% Adjacency matrix
% A_kl = (1 * 10^(-3)) * rand(K);
A_kl = [0.5 0.5; 0.5 0.5];

% Normalise Adjacency matrix such that columns sum to 1
for i = 1 : K
    col = A_kl(:,i);
    col = col(col~=A_kl(i,i));
    A_kl(i,i) = 1 - sum(col);
end
A_kl_init = A_kl;
% A_kl = [0.99 0.005 0.005 0;...
%      0.002 0.99 0 0.005;...
%      0.008 0 0.99 0.005;...
%      0 0.005 0.005 0.99];
% Initialise initial means by running K-means plus plus algorithm
[~, mu] = kmeans(X, K,"Start","plus"); 

% Number of data points and dimensions, respectively
[T,D] = size(X); 
% Initialise covariance matrices
sigma = repmat(eye(D), [1 1 K]);
% Initialise responsibility matrix: p(z_k(t_i) = 1|x(t_i)
gamma = zeros(K,T);
% Initialise prior matrix / mixing coefficient matrix: p(z_k(t_i) =
% 1|x(t_i-1))
pi = zeros(K,T); 
% Set responsibility of first datapoint assigned to first cluster to 1:
% p(z_1(t_1) = 1|x(t_1)) = 1
gamma(1,1) = 1; 
% Expectation-Maximisation algorithm
diff = 1;
while diff > tol && iter <= maxIter
    a = mu;
    b = A_kl;
    for i = 2:T
        % Expectation - step: Compute responsibilities
        % pi_k(t_i) = Akl*gamma_k(t_i-1)
        % K x 1 = (K x K)*(K x 1)
        pi(:,i) = A_kl*gamma(:,(i - 1));
        % Compute responsibilities: p(z_k|x) =
        % pi_k*N(x|mu_k,sigma_k)/sum(pi_k*N(x|mu_k,sigma_k))
        % First, the numerator
        gamma(:,i) = pi(:,i) .* mvnpdf(X(i,:), mu(:,:), sigma(:,:,:) + eye(D)*1e-12) ;
        
        % Marginalise responsibilities, dividing by the denominator
        B = nonzeros([sum(gamma(:,i), 1) 1]);
        gamma(:,i) = gamma(:,i)./B(1);
    end
    % Maximisation-step: Update parameters
    % Sum rows
    Nk = sum(gamma, 2); 
    % Update means
    % mu_new = sum(gamma*X)/sum(gamma)
    % K x D = (K x T)(T x D)
    mu = (gamma*X) ./ Nk;
    M = L - 1;
    % Update Covariances
    for k = 1:K
        % Mean-centered data
        Xc = X - mu(k,:);
        % Update covariance matrices
        % sigma_new = sum(gamma*(X-mean(X))(X - mean(X))')/sum(gamma)
        % D x D = (D x T)*(T x T)*(T x D)
        nchunks = T/2;
        loops = T/nchunks;
        temp = 0;
%         for i = 1 : loops
%             temp = temp + (Xc((i - 1)*nchunks + 1 : i * nchunks,:))' * spdiags(gamma(k,(i - 1)*nchunks + 1 : i * nchunks)', 0, nchunks, nchunks) * (Xc((i - 1)*nchunks + 1 : i * nchunks,:));
%         end
        for i = 1 : loops
            temp = temp + (Xc((i - 1)*nchunks + 1 : i * nchunks,:))' * diag(gamma(k,(i - 1)*nchunks + 1 : i * nchunks)') * (Xc((i - 1)*nchunks + 1 : i * nchunks,:));
        end
        temp = temp/Nk(k);
        sigma(:,:,k) = temp;
    end

    for k = 1:K

        % Update adjacency matrix, Akl
        % Akl,new =
        % sum_T(N(x|mu_k,sigma_k)*Akl,old/sum_K(N(x|mu_k,sigma_k)*Akl,old)*1/(T
        % -1)
        denominator = ones(T-1,1)*1E-15;
        for l = 1:L
            if l == 1
                % Calculate the denominator
                for m = 0:M
                    denominator = denominator + mvnpdf(X(1:end-1,:), mu(l + m,:), sigma(:,:,l + m) + eye(D)*1e-12)*A_kl(l + m,k); 
                end
            end
          
            % Calculate the new matrix
            A_kl(l,k) = (1/(T -1))*sum((mvnpdf(X(1:end-1,:), mu(k,:), sigma(:,:,k) + eye(D)*1e-12)*A_kl(l,k))./...
                sum(denominator));
        end
    end
   A_kl = A_kl./(sum(A_kl,1));
   diff = sum(sum(abs(mu - a))) + sum(sum(abs(A_kl - b))); 
   iter = iter + 1;
end

[~,C] = max(gamma,[],1);
C = C';
pi = (Nk/T)';
end

