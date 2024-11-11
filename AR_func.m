function [out] = AR_func(alg,change)
%% Clustering of Auto-correlated data
%  Frikkie Olivier
rng(1)
%% Define time region of interest
p.t = 1000; % Time span
t = linspace(1,p.t,p.t);

%% Define parameters
if change ==1
    % Variance change, constant mean
    p.d0.Temperature = @(t) 300; % Mean temperature [K]
    p.d0.Pressure = @(t) 150*10^3; % Mean pressure [Pa]
    p.phi.Temperature = @(t) 0; % Auto-correlation
    p.phi.Pressure = @(t) 0; % Auto-correlation
    p.sigma.Temperature = @(t) 0.08 + 0.1*(t > 300)...
                                - 0.1*(t > 600); % Variance
    p.sigma.Pressure = @(t) 0.1 + 0.1*(t > 300)...
                              - 0.1*(t > 600); % Variance
else
    % Mean change, constant variance
    p.d0.Temperature = @(t) 300 + 0.1*300*(t > 300)...
                                - 0.1*300*(t > 600); % Mean temperature [K]
    p.d0.Pressure = @(t) 150*10^3 + 0.1*150*10^3*(t > 300)...
                              - 0.1*150*10^3*(t > 600); % Mean pressure [Pa]
    p.phi.Temperature = @(t) 0; % Auto-correlation
    p.phi.Pressure = @(t) 0; % Auto-correlation
    p.sigma.Temperature = @(t) 0.08;  % Variance
    p.sigma.Pressure = @(t) 0.1;  % Variance
end
    p.data_N = 2;

data_fields = ["Temperature", "Pressure"];
cf = data_fields;
for i = 1 : p.data_N
    dd.(cf(i)) = zeros(p.t, 1);
    dd.(cf(i))(1) = p.d0.(cf(i))(1);
end

for i = 2:p.t
    for j = 1:p.data_N
        dd.(cf(j))(i) = p.phi.(cf(j))(i)*dd.(cf(j))(i-1) + (1 - p.phi.(cf(j))(i))*p.d0.(cf(j))(i) + p.d0.(cf(j))(i)*p.sigma.(cf(j))(i)*sqrt(1 - p.phi.(cf(j))(i)^2)*randn; 
    end
end

y_vec = zeros(length(t),length(data_fields));
for i = 1:length(data_fields)
    d.(cf(i)) = griddedInterpolant(t, dd.(cf(i)));
    y_vec(:,i) = d.(cf(i)).Values;
end
data = y_vec;
std_data = (data-mean(data))./std(data);
[~, score, ~, ~, ~] = pca(std_data);
X = [t' score(:,1) score(:,2)];

%% Ground Truth
C_T = zeros(length(score),1);
C_T(1:300) = 1;
C_T(301:600) = 2;
C_T(601:1000) = 1;
%% Optimisation
switch alg
    case 1 %TCDBSCAN,125 loops
        tic
        J_best = -10^6;
        C_best = [];
        eps_0 = 0.2:0.2:1;
        minPts_0 = 12:2:20;
        tau_0 = 0.02:0.02:0.10;
        for i = 1 : length(eps_0)
            for j = 1 : length(minPts_0)
                for k = 1 : length(tau_0)
%                     x = fminsearch(@(x) -sum_cluster_transitions(TCDBSCAN(X,x(1),max(1,abs(round(x(2)))),x(3))),[eps_0(i);minPts_0(j);tau_0(k)])
                    x = fminsearch(@(x) -sum_silhouette(X(:,2:3), TCDBSCAN(X,x(1),max(1,abs(round(x(2)))),x(3))),[eps_0(i);minPts_0(j);tau_0(k)])

                    eps = x(1);
                    minPts = max(1,abs(round(x(2))));
                    tau = x(3);
                    C_TCD = TCDBSCAN(X,eps,minPts,tau);
%                     J_TCD = sum_cluster_transitions(C_TCD);
                    J_TCD = sum_silhouette(X(:,2:3), C_TCD);
         
                    if J_TCD > J_best
                        J_best = J_TCD;
                        C_best = C_TCD;
                    end
                end
            end
        end
        time = toc;
        out.time = time;
        C_TCD = C_best;
        out.C_TCD = C_TCD;
        out.C_T = C_T;
        out.score = score;
    
    case 2 %DBSCAN, 100 loops
        tic
        J_best = -10^6;
        C_best = [];
        eps_0 = 0.1:0.1:1;
        minPts_0 = 11:20;
        for i = 1 : length(eps_0)
            for j = 1 : length(minPts_0)
                    x = fminsearch(@(x) -sum_silhouette(X(:,2:3), dbscan(X(:,2:3),x(1),max(1,abs(round(x(2)))))),[eps_0(i);minPts_0(j)])
%                   x = fminsearch(@(x) -sum_cluster_transitions(dbscan(X(:,2:3),x(1),max(1,abs(round(x(2)))))),[eps_0(i);minPts_0(j)])
                    eps = x(1);
                    minPts = max(1,abs(round(x(2))));
                    C_D = dbscan(X(:,2:3),eps,minPts);
                    J_D = similarity(X(:,2:3), C_D);
%                     J_D = sum_cluster_transitions(C_D);
                    if J_D > J_best
                        J_best = J_D;
                        C_best = C_D;
                    end
             end
        end
        time = toc;
        out.time = time;
        C_D = C_best;
        out.C_D = C_D;
        out.C_T = C_T;
        out.score = score;
    case 3 % TCKMEANS , 100 initialisations
        tic
        J_best = -10^6;
        K = 2;
        Tol = 1e-5;
        for i = 1 : 100
            rng(i)
          tau_K = fminbnd(@(tau_K) -sum_cluster_transitions(TCKMEANS2(X, K, Tol, tau_K)), 0.2, 1);
%             tau_K = fminbnd(@(tau_K) -sum_silhouette(X(:,2:3), TCKMEANS2(X, K, Tol, tau_K)), 0.2, 1);

            C_TCK = TCKMEANS2(X, K, Tol, tau_K);
          J_TCK = sum_cluster_transitions(C_TCK);
%             J_TCK = sum_silhouette(X(:,2:3),C_TCK);
%           
            if J_TCK > J_best
                J_best = J_TCK;
                C_best = C_TCK;
            end
        end
        time = toc;
        out.time = time;
        C_TCK = C_best;
        C_TCK = ac_label(C_T,C_TCK,zeros(2,2));
        out.C_TCK = C_TCK;
        out.C_T = C_T;
        out.score = score;      
        
    case 4 % K-means, 100 initialisations
        tic
        K = 2;
        J_best = -10^6;

        for i = 1 : 100
            rng(i)
            C_K = kmeans(X(:,2:3),K,"Start","plus");
            J_K = sum_cluster_transitions(C_K);
%             J_K = sum_silhouette(X(:,2:3),C_K);
            if J_K > J_best
                J_best = J_K;
                C_best = C_K;
            end
        end
        time = toc;
        out.time = time;  
        C_K = C_best;
        C_K = ac_label(C_T,C_K,zeros(2,2));
        out.C_K = C_K;
        out.score = score;
        out.C_T = C_T;       
        
    case 5 % TCGMM, 100 intialisations
        tic
        J_best = -10^6;
        K = 2;       
        tol = 1e-10;
        for j = 1 : 100
            rng(j)
            [C_TCGMM,~,~,~, A_kl ] = TCGMM(X(:,2:3), K, tol);
          J_TCGMM = sum_cluster_transitions(C_TCGMM);
%             J_TCGMM = sum_silhouette(X(:,2:3),C_TCGMM);
            if J_TCGMM > J_best
                J_best = J_TCGMM;
                C_best = C_TCGMM;
                A_kl_best = A_kl;
            end
        end
        time = toc;
        out.time = time;
        C_TCGMM = C_best;
        C_TCGMM = ac_label(C_T,C_TCGMM,zeros(2,2));
        out.C_TCGMM = C_TCGMM;
        out.score = score;
        out.C_T = C_T;
    case 6 %GMM 100 intialisations
        tic
        J_best = -10^6;
        K = 2;
    
        for i = 1 : 100
            rng(i)
            C_GMM = GMM(X(:,2:3),K,100);
            J_GMM = sum_cluster_transitions(C_GMM);
%           J_GMM = similarity(X(:,2:3),C_GMM);
                if J_GMM > J_best
                    J_best = J_GMM;
                    C_best = C_GMM;
                end
        end
        time = toc;
        out.time = time;
        C_GMM = C_best;
        C_GMM = ac_label(C_T,C_GMM,zeros(2,2));
        out.C_GMM = C_GMM;
        out.score = score;
        out.C_T = C_T;
        
end
end
