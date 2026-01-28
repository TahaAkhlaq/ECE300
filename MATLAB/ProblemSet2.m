% Taha Akhlaq - Communication Theory Problem Set 2 Question 4 (A–D)
clc; clear; close all;

function color = pick_ball(r,b)
    color = 1 + (randi(r+b) > r);  % 1=Red if draw<=r, else 2=Blue
end

% Part A
% deterministic rule (1=Red, 2=Blue) -> [g_R, g_B]
rule = uint8([1, 2]);

% Parts B-D
N = 1e5;
% [r1 b1 r2 b2]
cases = [ ...          
     8  2  8  2;  % Case 1
     4  6  8  2;  % Case 2
     8  2  4  6;  % Case 3
     4  6  4  6]; % Case 4

for c = 1:size(cases,1)
    r1 = cases(c,1);  b1 = cases(c,2);
    r2 = cases(c,3);  b2 = cases(c,4);

    % Part B
    pX = [r1; b1] / (r1 + b1); % prior P(X=Red), P(X=Blue)
    T  = r2 + b2 + 1; % Urn II size after transfer
    pYgX = [ (r2+1)/T,  b2/T ; % P(Y|X)
             r2/T,    (b2+1)/T ];
    postR = (pYgX(:,1).*pX) / sum(pYgX(:,1).*pX); % P(X|Y=Red)
    postB = (pYgX(:,2).*pX) / sum(pYgX(:,2).*pX); % P(X|Y=Blue)

    [~,iR] = max(pYgX(:,1)); [~,iB] = max(pYgX(:,2)); rule_ML  = uint8([iR iB]); % ML
    [~,iR] = max(postR); [~,iB] = max(postB); rule_MAP = uint8([iR iB]); % MAP

    % theoretical probability of error for a rule
    perr = @(rule) ...
        pX(1)*((rule(1)~=1)*pYgX(1,1) + (rule(2)~=1)*pYgX(1,2)) + ...
        pX(2)*((rule(1)~=2)*pYgX(2,1) + (rule(2)~=2)*pYgX(2,2));
    Perr_ML_th  = perr(rule_ML);
    Perr_MAP_th = perr(rule_MAP);


    % Part C (time-efficient)
    Xi = randi(r1+b1, N, 1); % draw X from Urn I
    X  = uint8(1 + (Xi > r1)); % 1=Red, 2=Blue
    Ti = r2 + b2 + 1;  % Urn II size after transfer
    Yi = randi(Ti, N, 1); % Y from Urn II
    redCap = (r2 + 1) * (X==1) + (r2) * (X==2); 
    Y  = uint8(1 + (Yi > redCap)); % 1=Red, 2=Blue

    % rules and empirical error
    Xhat_ML  = rule_ML(1)*uint8(Y==1) + rule_ML(2)*uint8(Y==2);
    Xhat_MAP = rule_MAP(1)*uint8(Y==1) + rule_MAP(2)*uint8(Y==2);
    Perr_ML_emp  = mean(Xhat_ML  ~= X);
    Perr_MAP_emp = mean(Xhat_MAP ~= X);

    % Part D (priors, likelihoods, posteriors, rules, errors)
    fprintf('Case %d: r1=%d b1=%d | r2=%d b2=%d\n', c, r1,b1,r2,b2);
    
    fprintf('  Priors        [P(R), P(B)]            = [%.3f %.3f]\n', pX(1), pX(2));
    fprintf('  P(Y|X) rows=X=[Red;Blue], cols Y=[Red Blue]\n');
    fprintf('                [%.3f %.3f; %.3f %.3f]\n', pYgX(1,1),pYgX(1,2),pYgX(2,1),pYgX(2,2));
    fprintf('  Post X|Y=Red  [P(R), P(B)]            = [%.3f %.3f]\n', postR(1), postR(2));
    fprintf('  Post X|Y=Blue [P(R), P(B)]            = [%.3f %.3f]\n', postB(1), postB(2));
    
    fprintf('  ML  rule=[%u %u]  Perr(th)=%.5f  Perr(emp)=%.5f\n', ...
            rule_ML,  Perr_ML_th,  Perr_ML_emp);
    fprintf('  MAP rule=[%u %u]  Perr(th)=%.5f  Perr(emp)=%.5f\n', ...
            rule_MAP, Perr_MAP_th, Perr_MAP_emp);
    
    % comparison comment (theoretical)
    if Perr_MAP_th < Perr_ML_th
        fprintf('  MAP better than ML (uses priors when they are skewed)\n\n');
    elseif Perr_MAP_th > Perr_ML_th
        fprintf('  ML better than MAP (likelihood dominates here)\n\n');
    else
        fprintf('  MAP and ML tie (symmetric priors/likelihoods)\n\n');
    end
end

% explanation of target goal (c)
fprintf('Method: time-efficient -> vectorized randi draws and decisions (no per-trial loop) \n');
