% Taha Akhlaq - Communication Theory Problem Set 1 Question 4 (A-E)
clc; clear; close all;

%% Part B

% build constellations
S_psk8 = constellation_8psk(1); % 2 x 8
S_qam16 = constellation_16qam(1); % 2 x 16
S_orth4 = constellation_4orth(1); % 4 x 4

% check dmin
dm8 = min_distance(S_psk8);
dm16 = min_distance(S_qam16);
dm4 = min_distance(S_orth4);
fprintf('dmin:  8-PSK: %.6f  16-QAM: %.6f  4-orth: %.6f\n', dm8, dm16, dm4);

% compute
[Eb_dB_psk8,  eta_psk8,  Eb_psk8,  Es_psk8]  = constellation_metrics(S_psk8);
[Eb_dB_qam16, eta_qam16, Eb_qam16, Es_qam16] = constellation_metrics(S_qam16);
[Eb_dB_orth4, eta_orth4, Eb_orth4, Es_orth4] = constellation_metrics(S_orth4);

% output
fprintf('8-PSK:   Eb(dB/bit) = %+8.4f  eta = %.3f\n', Eb_dB_psk8,  eta_psk8);
fprintf('16-QAM:  Eb(dB/bit) = %+8.4f  eta = %.3f\n', Eb_dB_qam16, eta_qam16);
fprintf('4-ortho: Eb(dB/bit) = %+8.4f  eta = %.3f\n', Eb_dB_orth4, eta_orth4);

%% Part C

% table
Names  = ["8-PSK"; "16-QAM"; "4-orth"];
Nvec = [size(S_psk8,1);  size(S_qam16,1);  size(S_orth4,1)];
Mvec = [size(S_psk8,2);  size(S_qam16,2);  size(S_orth4,2)];
dmin_v = [dm8; dm16; dm4];
Es_v = [Es_psk8; Es_qam16; Es_orth4];
Eb_v = [Eb_psk8; Eb_qam16; Eb_orth4];
Eb_dB_v = [Eb_dB_psk8; Eb_dB_qam16; Eb_dB_orth4];
eta_v = [eta_psk8; eta_qam16; eta_orth4];

fprintf('\n');
T = table(Names, Nvec, Mvec, round(dmin_v,6), round(Es_v,6), round(Eb_v,6), round(Eb_dB_v,4), round(eta_v,3), ...
    'VariableNames', {'Constellation','N','M','dmin','Es','Eb','Eb_dB','eta'});
disp(T);

%% Part D

% max eta
[eta_max,i_eta_max] = max(eta_v);

fprintf('Most spectrally efficient: %s (eta = %.3f)\n',T.Constellation{i_eta_max},eta_max);

%% Part E

% power gap
[~,i_power_best] = min(Eb_dB_v);
[~,i_eta_min] = min(eta_v);
delta_dB = Eb_dB_v(i_eta_max) - Eb_dB_v(i_eta_min);

fprintf('Most power efficient: %s (Eb/bit = %.4f dB)\n',T.Constellation{i_power_best},Eb_dB_v(i_power_best));
fprintf('Eb/bit difference between most and least spectrally efficient: %.4f dB\n',delta_dB);

%% Part A

% compute Eb and spectral efficiency
function [Eb_dB, eta, Eb, Es] = constellation_metrics(S)

    % input checks
    if ~isnumeric(S) || ~isreal(S)
        error('S must be real numeric matrix');
    end

    [N, M] = size(S);
    if N < 1 || M < 2
        error('S must be N x M with N>=1 and M>=2');
    end

    k = log2(M);
    if ~isfinite(k) || k <= 0
        error('invalid M');
    end

    % energies
    Es = mean(sum(S.^2, 1));
    Eb = Es / k; % energy per bit
    if ~isfinite(Eb) || Eb <= 0
        error('nonpositive Eb');
    end

    % outputs
    Eb_dB = 10 * log10(Eb); % dB per bit
    eta = k / N; % bits per real dimension
end

%% Constellation Functions

% 8-PSK
function S = constellation_8psk(dmin)
    if nargin < 1
        dmin = 1;
    end
    M  = 8;
    A  = dmin / (2 * sin(pi / M));
    th = 2 * pi * (0:M-1) / M;
    S  = A * [cos(th); sin(th)];
end

% 16-QAM
function S = constellation_16qam(dmin)
    if nargin < 1
        dmin = 1;
    end
    base = [-3 -1 1 3];         
    s = dmin / 2;
    lv = s * base;
    [I, Q] = meshgrid(lv, lv);
    S = [I(:)'; Q(:)'];
end

% 4-orthogonal
function S = constellation_4orth(dmin)
    if nargin < 1
        dmin = 1;
    end
    a = dmin / sqrt(2);
    S = a * eye(4);
end

%% Distance Function

function d = min_distance(S)
    G  = S.' * S;
    g  = diag(G);
    D2 = g + g.' - 2 * G;
    D2(1:size(S,2)+1:end) = inf;
    d = sqrt(min(D2(:)));
end
