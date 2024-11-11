function J = sum_cluster_transitions(C_org)
         J = 0;
         for i = 2 : length(C_org)
             if C_org(i) ~= C_org(i - 1)
                 J = J + 1;
             end
         end
         J = (-1)*J;
end