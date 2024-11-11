function [out] = AR_modes_func(phi,M,V,method)
rng(1)
warning ('off','all');


%% Define time region of interest
p.t = 10000; % Time span
t = linspace(1,p.t,p.t);
% Mode transition matrix
alpha = 0.01; % Variance change
eps = 0.01; % Mean change
stay = 1 - (alpha*eps + eps + alpha);
A = [stay alpha eps alpha*eps;...
     alpha stay alpha*eps eps;...
     eps alpha*eps stay alpha;...
     alpha*eps eps alpha stay];

C = cumsum(A,1);
mu = zeros(2,length(t));
sigma = zeros(2,length(t));

for i = 1:p.t
    if i == 1
        mode = 1;
        modes(i) = mode;
    else
        % Current mode
        r(i) = rand;
        mode = find(C(:,mode) >= r(i), 1);
        modes(i) = mode;
    end

    switch mode
        case {1} % Mean 1, Variance 1
            mu(:,i) = [300;150*10^3];
            sigma(:,i) = [0.08;0.1];

        case {2} % Mean 1, Variance 2 (Variance change, const mean) alpha
            mu(:,i) = [300;150*10^3];
            sigma(:,i) = [0.08 + V;0.1 + V];
        case {3} % Mean 2, Variance 1 (Mean change, const variance) epsilon
            mu(:,i) = [(1+M)*300;(1+M)*150*10^3];
            sigma(:,i) = [0.08;0.1];
        case {4} % Mean 2, Variance 2 (Mean change, variance change) alpha*epsilon
            mu(:,i) = [1.1*300;1.1*150*10^3];
            sigma(:,i) = [0.18;0.2];
    end
end
uu = zeros(2,length(t));
uu(:,1) = mu(:,1);
% phi = 0.75;
minimum = [150;50*10^3];
maximum = [450;250*10^3];
for j = 2 : p.t
    for k = 1 : 2
        uu(k,j) =  min(max(phi*uu(k,j-1) + (1-phi)*mu(k,j) + mu(k,j)*sigma(k,j)*sqrt(1-phi^2)*randn, minimum(k)), maximum(k));
    end
end
fields =  ["T", "P"];
for i = 1 : length(fields)
    u.(fields(i)) = griddedInterpolant(linspace(1, p.t, length(uu)), uu(i,:));
end
y_vec = [];
y_vec(1,:) = u.T(t);
y_vec(2,:) = u.P(t);


data = y_vec';
std_data = (data-mean(data))./std(data);
score = std_data;
C_T = modes';

%% GMM
if method == 1
J_best = -10^6;
K = 4;
Tol = 1e-10;

for i = 1 : 100
                rng(i)
                C_GMM = GMM(score(:,1:2),K,100);
              J_GMM = sum(diag(confusionmat(C_T, C_GMM)));
%                 J_GMM = sum_cluster_transitions(C_GMM);
%                 J_GMM = similarity(score(:,1:2),C_GMM);
                if J_GMM > J_best
                    J_best = J_GMM;
                    C_best = C_GMM;
                end
end

C_GMM = C_best;
[C_GMM,~] = ac_label(C_T,C_GMM,zeros(4,4));
out.C_GMM = C_GMM;
out.score = score;
out.C_T = C_T;


else
    


J_best = -10^6;
K = 4;
tol = 1e-10;

for j = 1 : 100
                rng(j)
                [C_TCGMM,~,~,~, A_kl,A_kl_init ] = TCGMM(score(:,1:2), K, tol);
%               J_TCGMM = sum(diag(confusionmat(C_T, C_TCGMM)));
%                 J_TCGMM = sum_cluster_transitions(C_TCGMM);
                J_TCGMM = similarity(score(:,1:2),C_TCGMM);
                if J_TCGMM > J_best
                    J_best = J_TCGMM;
                    C_best = C_TCGMM;
                    A_kl_best = A_kl;
                    A_kl_init_best = A_kl_init;
                end
end
C_TCGMM = C_best;
[C_TCGMM,A_kl] = ac_label(C_T,C_TCGMM,A_kl_best);

out.C_TCGMM = C_TCGMM;
out.score = score;
out.A_kl = A_kl;
out.C_T = C_T;
out.A_kl_init_best = A_kl_init_best;

end

end

