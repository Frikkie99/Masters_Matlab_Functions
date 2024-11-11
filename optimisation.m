function [A_kl,C_best,J_TCGMM] = optimisation(score,K, PC, method, R)
tol = 1e-10;
C_best = [];
A_kl_best = [];

if method == 1
    for i = R
        rng(i)
        [C_TCGMM, ~, ~, ~, A_kl ] = TCGMM(score(:,1:PC), K, tol);
%         A_kl = A_kl./(sum(A_kl,1));
        J_TCGMM = sum_cluster_transitions(C_TCGMM);        
        
            
            C_best = C_TCGMM;
%             A_kl_best = A_kl;
      
   
    end
   else
    for i = R
        rng(i)
        [C_TCGMM, ~, ~, ~, A_kl ] = TCGMM(score(:,1:PC), K, tol);
%         A_kl = A_kl./(sum(A_kl,1));
        J_TCGMM = similarity(score(:,1:PC),C_TCGMM);
        
            
            C_best = C_TCGMM;
%             A_kl_best = A_kl;
       
   
    end
    
end
% data.A_kl = A_kl_best;
% data.C_TCGMM = C_best;
% data.J_TCGMM = J_TCGMM;
end

