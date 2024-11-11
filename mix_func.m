function [out] = mix_func(phi,sig,sigC,method,PC)
rng(3);

 
% Time step
p.dt = 1;

% Time interval
t = 0 : p.dt : 60*24*2;


p.tN = length(t);

 
% Auto-regressive functions parameters
mu = 2; 
% phi = 0; sig = 0.1;

muC = 50; phiC = phi;

 
% Inputs
Fin1 = mu;

Fin2 = mu;

Ca1 = muC;

invControl = 1;

for i = 2:p.tN

    % Minimum value for Fin1 is 0.1
    Fin1(i) = max(phi*Fin1(i-1) + (1-phi)*mu + sig*mu*randn, 0.1);

    %Fin2(i) = max(phi*Fin2(i-1) + (1-phi)*mu + sig*randn, 0.1);

    % Minimum value for Ca1 is 1
    Ca1(i) = max(phiC*Ca1(i-1) + (1-phiC)*muC + sigC*muC*randn, 1);

    % 1% chance for inventory control to transistion to level control
    if rand < 0.01

        invControl(i) = invControl(i-1)*-1;

    else

        invControl(i) = invControl(i-1);

    end

end

u.Fin1 = griddedInterpolant(t, 2/4*Fin1);

%u.Fin2 = griddedInterpolant(t, 1/4*Fin2);

u.Ca1 = griddedInterpolant(t, Ca1);

 

options = optimoptions('fmincon', 'Display', 'off');

 

% System parameters

p.A = 10;

p.cv = 5;

p.K = 0.1;

 

% Control parameters

% Measurement interval
p.N = 5;

% Level setpoint
p.hSP = @(t) 1;

% Concentration setpoint
p.CaSP = 30;

 

penalty = @(uCV, u0) 1*sum( (uCV - [u0; uCV(1:end-1)]).^2 );

servoJ = @(h, hSP, uCV, u0) sum( (h - hSP).^2 ) + penalty(uCV, u0);

inventoryJ = @(h, hSP, uCV, u0) sum( (2 - tanh(50*(h-0.5)) + tanh(50*(h-1.5))).*(h-1).^2 ) + penalty(uCV, u0);

 

% Model

dhdt = @(t, h, Ca, u, p) [(u.Fin1(t) + min(1, max(0.1, p.K*(Ca - p.CaSP))) - u.fv(t)*p.cv*sqrt(h)) / p.A; (u.Fin1(t)*u.Ca1(t)    - u.fv(t)*p.cv*sqrt(h)*Ca) / p.A];

 

tic

% Start simulation

y.Fin1 = u.Fin1(0);

y.Ca1 = u.Ca1(0);

y.h = 1;

y.fv = 0.5;

y.Ca = 30;

 

uCV = y.fv*ones(p.N,1);

for i = 1:p.tN-1

    if invControl(i) == 1

        p.fJ = inventoryJ;

    else

        p.fJ = servoJ;

    end

    [y, uCV] = integrateStep(i, t, uCV, dhdt, u, y, p, options);

end

y.J(i+1) = 0;

y.fv(i+1) = y.fv(i);

 

toc

 

ym.t = t';

ym.invControl = invControl';

ym.Fin1 = y.Fin1' + 0.05*randn(p.tN,1);

ym.Fin2 = max(0, min(1, max(0.1, p.K*(y.Ca - p.CaSP)))' + 0.05*randn(p.tN,1));

ym.Ca1 = y.Ca1' + randn(p.tN,1);

ym.Fout = (y.fv.*p.cv.*sqrt(y.h))' + 0.05*randn(p.tN,1);

ym.h = y.h' + 0.025*randn(p.tN,1);

ym.Ca = y.Ca' + randn(p.tN,1);

 



 

X = [ym.Fin1 ym.Fin2 ym.Ca1 ym.Fout ym.h ym.Ca];
y_fields = ["Fin1", "Ca1", "h", "fv", "Ca", "J"];
% data = struct2vec(y,y_fields);

data = X;
data = data(121:end,:);
invControl = invControl(121:end);
std_data = (data-mean(data))./std(data);

[coeff, score, latent, ~, explained] = pca(std_data);
% plot(t(121:end),score(:,1),t(ym.invControl==1),ones(size(t(ym.invControl==1)))*0.9,'r*')
C_T = zeros(length(invControl),1);
C_T(invControl == 1) = 1;
C_T(invControl == -1) = 2;
if method == 1
    J_best = -10^6;
    K = 2;
    

for i = 1 : 20
                rng(i)
                C_GMM = GMM(score(:,1:PC),K,100);
%                 [C_GMM,~] = ac_label(C_T,C_GMM,zeros(2,2));
%               J_GMM = sum(diag(confusionmat(C_T, C_GMM)));
%                 J_GMM = sum_cluster_transitions(C_GMM);
                J_GMM = similarity(score(:,1:2),C_GMM);
                if J_GMM > J_best
                    J_best = J_GMM;
                    C_best = C_GMM;
                end
end

C_GMM = C_best;

out.C_GMM = C_GMM;
out.score = score;
out.C_T = C_T;
else
    J_best = -10^6;
    K = 2;
    tol = 1e-10;
for j = 1 : 20
                rng(j)
                [C_TCGMM,~,~,~, A_kl,A_init ] = TCGMM(score(:,1:PC), K, tol);
%               [C_TCGMM,A_kl] = ac_label(C_T,C_TCGMM,A_kl);  
%               J_TCGMM = sum(diag(confusionmat(C_T, C_TCGMM)));
%                 J_TCGMM = sum_cluster_transitions(C_TCGMM);
                J_TCGMM = similarity(score(:,1:2),C_TCGMM);
                if J_TCGMM > J_best
                    J_best = J_TCGMM;
                    C_best = C_TCGMM;
                    A_kl_best = A_kl;
                    A_init_best = A_init;
                end
end

C_TCGMM = C_best;
A_kl = A_kl_best;
A_init = A_init_best;
out.C_TCGMM = C_TCGMM;
out.A_kl = A_kl;
out.score = score;
out.C_T = C_T;
end


function J = objectiveMPC(uCV, t, dhdt, y, p)

    tHorizon = t : p.dt : t+(p.N-1)*p.dt;

    u.fv = griddedInterpolant(tHorizon, uCV, "previous");

    u.Fin1 = @(t) y.Fin1(end);

    %u.Fin2 = @(t) y.Fin2(end);

    u.Ca1  = @(t) y.Ca1(end);

 

    [~, x] = ode23(@(t, x) dhdt(t, x(1), x(2), u, p), tHorizon, [y.h(end); y.Ca(end)]);

   

    J = p.fJ(x(:,1), p.hSP(tHorizon)', uCV, y.fv(end));

end

 

function [y, uCV] = integrateStep(i, t, uCV, dhdt, u, y, p, options)

    [uCV, y.J(i)] = fmincon(@(uCV) objectiveMPC(uCV, t(i), dhdt, y, p), uCV, [], [], [], [], zeros(p.N,1), ones(p.N,1), [], options);

   

    u.fv = @(t) uCV(1);

    y.fv(i) = uCV(1);

    sol = ode45(@(t,x) dhdt(t,x(1),x(2),u,p), [t(i) t(i+1)], [y.h(end); y.Ca(end)]);

    y.h(i+1) =sol.y(1,end);

    y.Ca(i+1) = sol.y(2,end);

    y.Fin1(i+1) = u.Fin1(t(i+1));

    %y.Fin2(i+1) = u.Fin2(t(i+1));

    y.Ca1(i+1) = u.Ca1(t(i+1));

end
end

