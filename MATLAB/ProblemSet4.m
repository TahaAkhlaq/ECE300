%Q3
% Parameters
Ac = 2; % carrier amplitude [V]
fc = 10; % carrier frequency [Hz]
fm = 2; % message frequency [Hz]
ka = 2; % amplitude sensitivity [V/V]
fs = 100; % sampling rate [Hz]
N = 2000; % number of samples
n = 0:N-1;
t = n/fs;

Am_vec = [0.40 0.45 0.60];
idx_plot = 1:400;

for k = 1:length(Am_vec)
    Am = Am_vec(k);

    % message and envelope term
    m    = Am * sin(2*pi*fm*t);
    env  = Ac * (1 + ka*m);

    % AM signal
    xAM  = env .* cos(2*pi*fc*t);

    figure;
    plot(t(idx_plot), xAM(idx_plot), 'LineWidth', 1);
    hold on;
    plot(t(idx_plot),  env(idx_plot), '--', 'LineWidth', 1);
    plot(t(idx_plot), -env(idx_plot), '--', 'LineWidth', 1);
    hold off;

    grid on;
    xlabel('t [s]');
    ylabel('x_{AM}(t)');
    title(sprintf('AM signal, A_m = %.2f V', Am_vec(k)));
    ylim([-5 5]);
end
