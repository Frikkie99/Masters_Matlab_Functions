function [comb] = opt_comb(S)

% Number of elements
n = size(S, 1);

% Objective function coefficients (to maximize sum)
f = -reshape(S, [], 1); % Reshape S into a column vector and negate for maximization

% Define linear constraints: each element can only be chosen once
% Aeq = zeros(n, n^2);
% beq = ones(n, 1);
% for i = 1:n
%     Aeq(i, (i-1)*n+1:i*n) = 1;
% end
% 
% % Additional constraints: ensure that each selected element comes from a different row and column
% Aeq = [Aeq; kron(eye(n), ones(1, n))]; % Each row must have only one selected element
% Aeq = [Aeq; kron(ones(1, n), eye(n))]; % Each column must have only one selected element
% beq = [beq; ones(2*n, 1)]; % Additional equality constraints

%TEST
Aeq = zeros(2*n,n^2);
beq = ones(2*n,1);
for i = 1:n
    Aeq(i, (i-1)*n+1:i*n) = 1;
    Aeq(n+1:end, (i-1)*n+1:i*n) = eye(n);
end



% Binary integer variables (0 or 1 indicating whether the element is chosen)
intcon = 1:n^2;

% Bounds on variables (0 or 1)
lb = zeros(n^2, 1);
ub = ones(n^2, 1);

% Solve the integer linear programming problem
[x, fval, exitflag] = intlinprog(f, intcon, [], [], Aeq, beq, lb, ub);

% Reshape the solution vector into a matrix
X = reshape(x, n, n);

% Extract the selected elements from the matrix S
selected_elements = S(logical(X));
comb = [];
for i = 1: length(selected_elements)
    comb(i) = find(1==X(:,i));
end
 end

