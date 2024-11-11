function [C_best,A_kl_best] = ac_label(C_T, C_org, A_kl_org)
labels = unique(C_org);
S = zeros(length(labels),length(labels));
for i = 1 : length(labels)
    for j = 1 : length(labels)
        S(i,j) = sum((C_org==i).*(C_T==j));
    end
end
comb = opt_comb(S);
C = zeros(size(C_org));
    for j = 1 : length(labels)
        C(C_org == comb(1,j)) = j;
    end
C_best = C;
new_labels = unique(C_best,"rows","stable");
j = 1;
A_kl_best = A_kl_org;
while j <= length(new_labels)
    if new_labels(j) ~= labels(j) && new_labels(j)~= 0
        A_kl_best(new_labels(j),new_labels(j)) = A_kl_org(labels(j),labels(j));
        A_kl_best(new_labels(j),labels(j)) = A_kl_org(labels(j),new_labels(j));
        A_kl_best(labels(j),new_labels(j)) = A_kl_org(new_labels(j),labels(j));
        A_kl_best(labels(j),labels(j)) = A_kl_org(new_labels(j),new_labels(j));

        new_labels(new_labels == labels(j)) = 0;
    end
    j = j + 1;
end
end