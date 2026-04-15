function R = gi_signal_enhancement_v4(opts)
% GI_SIGNAL_ENHANCEMENT_V4
% FINAL VERSION - Gaussian only
% MATLAB R2016a compatible
% Journal-ready figure export for Figure 4
%
% Usage:
%   R = gi_signal_enhancement_v4;
%   R = gi_signal_enhancement_v4(struct('show_plots',true));
%
% Output directory:
%   C:\Users\esray\Desktop\comparison_figures

if nargin < 1
    opts = struct();
end

clc;
close all;
rng(0);

% -------------------------------------------------------------------------
% 0) Synthetic signal generation
% -------------------------------------------------------------------------
N = 1024;
x_true = make_synthetic_signal(N);      % [0,1]
h = gauss_kernel(9, 2);
h = h(:);

Hx = conv_circ(x_true, h);

% Gaussian degradation only
sigma_input = 0.01;
b = clip01(Hx + sigma_input * randn(size(Hx)));

% -------------------------------------------------------------------------
% 1) Parameters
% -------------------------------------------------------------------------
P = struct();

% Algorithm parameters
P.maxit          = 500;
P.tau            = 0.06;
P.alpha          = 0.08;
P.stop_eps       = 1e-4;
P.stop_patience  = 8;
P.psnr_maxI      = 1;

% Ishikawa-type mixing
P.beta_a         = 8;
P.beta_b         = 60;
P.beta_min       = 0.03;
P.beta_max       = 0.15;

P.sigma_a        = 8;
P.sigma_b        = 60;
P.sigma_min      = 0.03;
P.sigma_max      = 0.15;

% Figure/export settings
P.show_plots       = true;
P.export_plots     = true;
P.output_dir       = 'C:\Users\esray\Desktop\comparison_figures';

P.fig_width_cm     = 18.0;
P.fig_height_cm    = 14.0;

P.font_name        = 'Arial';
P.font_size_axes   = 10;
P.font_size_label  = 11;
P.font_size_title  = 12;
P.font_size_legend = 9;
P.font_size_panel  = 13;

P.line_width_main  = 1.8;
P.line_width_minor = 1.2;
P.axes_line_width  = 1.0;
P.export_dpi       = 600;

% User overrides
fn = fieldnames(opts);
for k = 1:numel(fn)
    P.(fn{k}) = opts.(fn{k});
end

% Safe scaling
x_true = rescale01(double(x_true));
b      = rescale01(double(b));

% -------------------------------------------------------------------------
% 2) Initialization
% -------------------------------------------------------------------------
x = b;

PSNR = nan(P.maxit,1);
SSIM = nan(P.maxit,1);
J    = nan(P.maxit,1);

bestPSNR = -Inf;
bestx = x;
best_iter = 1;
stop_iter = P.maxit;

% -------------------------------------------------------------------------
% 3) Iterative process
% -------------------------------------------------------------------------
for n = 1:P.maxit

    beta_n  = clamp(P.beta_a  / (n + P.beta_b),  P.beta_min,  P.beta_max);
    sigma_n = clamp(P.sigma_a / (n + P.sigma_b), P.sigma_min, P.sigma_max);

    % C-operator
    Cx = grad_step_gaussian(x, b, h, P);
    y  = (1 - beta_n) * x + beta_n * Cx;

    % Backtracking for I-operator
    if n == 1 || ~isfinite(J(n-1))
        J_prev = inf;
    else
        J_prev = J(n-1);
    end

    local_tau = P.tau;
    max_bt = 8;
    bt = 0;

    while true
        Ptmp = P;
        Ptmp.tau = local_tau;

        Iy = grad_step_gaussian(y, b, h, Ptmp);
        x_cand = (1 - sigma_n) * x + sigma_n * Iy;

        J_cand = obj_val_gaussian(x_cand, b, h, P);
        if ~isfinite(J_cand)
            J_cand = inf;
        end

        if ~isfinite(J_prev) || (J_cand <= J_prev * (1 + 1e-5)) || bt >= max_bt
            x_new = x_cand;
            J(n) = J_cand;
            break;
        else
            local_tau = 0.5 * local_tau;
            bt = bt + 1;
        end
    end

    % Metrics
    PSNR(n) = psnr_safe(x_new, x_true, P.psnr_maxI);
    SSIM(n) = ssim1d(x_new, x_true, P.psnr_maxI);

    % Divergence protection
    if any(~isfinite([J(n), PSNR(n), SSIM(n)])) || J(n) > 1e12
        fprintf('[WARN] Divergence/overflow at iteration %d. Stopping.\n', n);
        x = x_new;
        stop_iter = n;
        break;
    end

    % Best PSNR
    if PSNR(n) > bestPSNR
        bestPSNR = PSNR(n);
        bestx = x_new;
        best_iter = n;
    end

    % Early stopping
    if n > P.stop_patience
        if (n - best_iter) >= P.stop_patience
            fprintf('Early stop at iteration %d (no PSNR improvement for %d iterations).\n', ...
                n, P.stop_patience);
            x = x_new;
            stop_iter = n;
            break;
        end

        Jw  = J(n-P.stop_patience+1:n);
        rel = abs(diff(Jw)) ./ max(abs(Jw(1:end-1)), eps);
        rel = rel(~isnan(rel) & ~isinf(rel));

        if ~isempty(rel) && mean(rel) < P.stop_eps
            fprintf('Early stop at iteration %d (flat objective).\n', n);
            x = x_new;
            stop_iter = n;
            break;
        end
    end

    x = x_new;
    stop_iter = n;
end

it_end = stop_iter;
PSNR = PSNR(1:it_end);
SSIM = SSIM(1:it_end);
J    = J(1:it_end);

if best_iter > numel(PSNR)
    best_iter = numel(PSNR);
end

% -------------------------------------------------------------------------
% 4) Figure creation and export
% -------------------------------------------------------------------------
saved_files = {};
caption_text = '';

if P.show_plots || P.export_plots
    [saved_files, caption_text] = create_plos_figure( ...
        x_true, b, bestx, J, PSNR, SSIM, best_iter, it_end, P);
end

% -------------------------------------------------------------------------
% 5) Output structure
% -------------------------------------------------------------------------
R = struct();
R.x_true       = x_true;
R.b            = b;
R.x            = bestx;
R.PSNR         = PSNR;
R.SSIM         = SSIM;
R.J            = J;
R.best_iter    = best_iter;
R.stop_iter    = it_end;
R.params       = P;
R.saved_files  = saved_files;
R.caption_text = caption_text;

fprintf('[GAUSSIAN] Ended at iteration %d | Best PSNR = %.2f dB (iteration %d) | Final SSIM = %.4f\n', ...
    it_end, bestPSNR, best_iter, SSIM(end));

end

% =========================================================================
% Figure creation
% =========================================================================
function [saved_files, caption_text] = create_plos_figure(x_true, b, x_best, J, PSNR, SSIM, best_iter, stop_iter, P)

saved_files = {};
iters = 1:numel(J);

i1 = 300;
i2 = 500;
i1 = max(1, min(i1, numel(x_true)));
i2 = max(i1, min(i2, numel(x_true)));

% Print-safe colors
c_true = [0.00 0.00 0.00];
c_obs  = [0.55 0.55 0.55];
c_rest = [0.00 0.45 0.74];
c_obj  = [0.49 0.18 0.56];
c_psnr = [0.85 0.33 0.10];
c_ssim = [0.00 0.60 0.50];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1 1 P.fig_width_cm P.fig_height_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 P.fig_width_cm P.fig_height_cm], ...
    'PaperSize', [P.fig_width_cm P.fig_height_cm], ...
    'Name', 'Figure 4 - Gaussian signal enhancement', ...
    'NumberTitle', 'off');

% ---------------- Panel A ----------------
ax1 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.08 0.60 0.39 0.25]);
hold(ax1, 'on');
plot(ax1, x_true, 'Color', c_true, 'LineWidth', P.line_width_main);
plot(ax1, b,      'Color', c_obs,  'LineWidth', P.line_width_minor);
plot(ax1, x_best, 'Color', c_rest, 'LineWidth', P.line_width_main);
hold(ax1, 'off');
set_common_axes_style(ax1, P);
xlim(ax1, [1 numel(x_true)]);
ylim(ax1, [-0.02 1.02]);
xlabel(ax1, 'Sample index', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
ylabel(ax1, 'Normalized amplitude', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
title(ax1, 'Signal restoration', 'FontName', P.font_name, 'FontSize', P.font_size_title, 'FontWeight', 'bold');

leg = legend(ax1, {'Ground truth', 'Blurred + noisy observation', 'Restored signal'});
set(leg, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'Box', 'off', ...
    'FontName', P.font_name, ...
    'FontSize', P.font_size_legend);

text(0.02, 0.05, sprintf('Best iter = %d', best_iter), ...
    'Parent', ax1, ...
    'Units', 'normalized', ...
    'FontName', P.font_name, ...
    'FontSize', P.font_size_legend, ...
    'BackgroundColor', [1 1 1], ...
    'EdgeColor', [0.80 0.80 0.80], ...
    'Margin', 3);

% ---------------- Panel B ----------------
ax2 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.55 0.60 0.37 0.25]);
hold(ax2, 'on');
plot(ax2, i1:i2, x_true(i1:i2), 'Color', c_true, 'LineWidth', P.line_width_main);
plot(ax2, i1:i2, b(i1:i2),      'Color', c_obs,  'LineWidth', P.line_width_minor);
plot(ax2, i1:i2, x_best(i1:i2), 'Color', c_rest, 'LineWidth', P.line_width_main);
hold(ax2, 'off');
set_common_axes_style(ax2, P);
xlim(ax2, [i1 i2]);
ylim(ax2, [-0.02 1.02]);
xlabel(ax2, 'Sample index', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
ylabel(ax2, 'Normalized amplitude', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
title(ax2, sprintf('Zoomed view (%d-%d)', i1, i2), ...
    'FontName', P.font_name, 'FontSize', P.font_size_title, 'FontWeight', 'bold');

% ---------------- Panel C ----------------
ax3 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.08 0.12 0.24 0.26]);
plot(ax3, iters, J, 'Color', c_obj, 'LineWidth', P.line_width_main);
set_common_axes_style(ax3, P);
xlim(ax3, [1 numel(J)]);
xlabel(ax3, 'Iteration', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
ylabel(ax3, 'Objective value J(x)', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
title(ax3, 'Objective convergence', 'FontName', P.font_name, 'FontSize', P.font_size_title, 'FontWeight', 'bold');

% ---------------- Panel D ----------------
ax4 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.39 0.12 0.24 0.26]);
plot(ax4, iters, PSNR, 'Color', c_psnr, 'LineWidth', P.line_width_main);
hold(ax4, 'on');
plot(ax4, best_iter, PSNR(best_iter), 'o', ...
    'MarkerSize', 7, ...
    'MarkerEdgeColor', c_psnr, ...
    'MarkerFaceColor', [1 1 1], ...
    'LineWidth', 1.2);
set_common_axes_style(ax4, P);
xlim(ax4, [1 numel(PSNR)]);
xlabel(ax4, 'Iteration', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
ylabel(ax4, 'PSNR (dB)', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
title(ax4, 'PSNR improvement', 'FontName', P.font_name, 'FontSize', P.font_size_title, 'FontWeight', 'bold');

yl4 = get(ax4, 'YLim');
plot(ax4, [stop_iter stop_iter], yl4, '--k', 'LineWidth', 1.0);

text(best_iter, PSNR(best_iter), sprintf('  %.2f dB', PSNR(best_iter)), ...
    'Parent', ax4, 'FontName', P.font_name, 'FontSize', P.font_size_legend);
text(stop_iter, yl4(2) - 0.03*(yl4(2)-yl4(1)), 'Early stop', ...
    'Parent', ax4, 'FontName', P.font_name, 'FontSize', P.font_size_legend, ...
    'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1]);
hold(ax4, 'off');

% ---------------- Panel E ----------------
ax5 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.69 0.12 0.23 0.26]);
plot(ax5, iters, SSIM, 'Color', c_ssim, 'LineWidth', P.line_width_main);
set_common_axes_style(ax5, P);
xlim(ax5, [1 numel(SSIM)]);
ylim(ax5, [max(0, min(SSIM)-0.03) min(1, max(SSIM)+0.03)]);
xlabel(ax5, 'Iteration', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
ylabel(ax5, 'SSIM', 'FontName', P.font_name, 'FontSize', P.font_size_label, 'FontWeight', 'bold');
title(ax5, 'SSIM evolution', 'FontName', P.font_name, 'FontSize', P.font_size_title, 'FontWeight', 'bold');

% Panel labels
annotation(fig, 'textbox', [0.03 0.84 0.03 0.03], 'String', 'A', ...
    'LineStyle', 'none', 'FontName', P.font_name, 'FontSize', P.font_size_panel, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.50 0.84 0.03 0.03], 'String', 'B', ...
    'LineStyle', 'none', 'FontName', P.font_name, 'FontSize', P.font_size_panel, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.03 0.36 0.03 0.03], 'String', 'C', ...
    'LineStyle', 'none', 'FontName', P.font_name, 'FontSize', P.font_size_panel, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.34 0.36 0.03 0.03], 'String', 'D', ...
    'LineStyle', 'none', 'FontName', P.font_name, 'FontSize', P.font_size_panel, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.64 0.36 0.03 0.03], 'String', 'E', ...
    'LineStyle', 'none', 'FontName', P.font_name, 'FontSize', P.font_size_panel, 'FontWeight', 'bold');

caption_text = build_caption(best_iter, stop_iter, PSNR(best_iter), SSIM(end));

if P.export_plots
    outdir = P.output_dir;
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end

    base_name = fullfile(outdir, 'Fig4_Gaussian_Main');

    % MATLAB R2016a safe export
    figure(fig);
    set(fig, 'PaperPositionMode', 'auto');

    print('-dtiff',  [base_name '.tif'], ['-r' num2str(P.export_dpi)]);
    print('-dpng',   [base_name '.png'], ['-r' num2str(P.export_dpi)]);
    print('-depsc2', '-painters', [base_name '.eps']);
    print('-dpdf',   '-painters', [base_name '.pdf']);

    fid = fopen(fullfile(outdir, 'Fig4_Gaussian_Caption.txt'), 'w');
    if fid ~= -1
        fprintf(fid, '%s\r\n', caption_text);
        fclose(fid);
    end

    saved_files = { ...
        [base_name '.tif'], ...
        [base_name '.png'], ...
        [base_name '.eps'], ...
        [base_name '.pdf'], ...
        fullfile(outdir, 'Fig4_Gaussian_Caption.txt')};
end

end

function set_common_axes_style(ax, P)
set(ax, ...
    'Box', 'off', ...
    'LineWidth', P.axes_line_width, ...
    'FontName', P.font_name, ...
    'FontSize', P.font_size_axes, ...
    'TickDir', 'out', ...
    'Layer', 'top', ...
    'XColor', [0 0 0], ...
    'YColor', [0 0 0], ...
    'Color', 'none');
grid(ax, 'on');
end

function txt = build_caption(best_iter, stop_iter, best_psnr, final_ssim)
txt = sprintf([ ...
    'Fig 4. Restoration and convergence behavior of the proposed two-step G-I-nonexpansive mapping for one-dimensional signal enhancement under Gaussian blur and additive Gaussian noise (sigma = 0.01). ' ...
    'Panel A compares the ground-truth signal, degraded observation, and restored signal. ' ...
    'Panel B presents a zoomed view of the restoration over samples 300-500. ' ...
    'Panel C shows the steady decrease of the objective value J(x). ' ...
    'Panel D illustrates the PSNR increase across iterations; the best performance is achieved at iteration %d with PSNR = %.2f dB, and the vertical dashed line marks the stopping iteration (%d). ' ...
    'Panel E shows the SSIM evolution, which remains at a high level throughout the iterative process. ' ...
    'The final SSIM is %.4f.'], ...
    best_iter, best_psnr, stop_iter, final_ssim);
end

% =========================================================================
% Core model
% =========================================================================
function y = grad_step_gaussian(x, b, h, P)
r = conv_circ(x, h) - b;
Ht_r = conv_circ(r, flipud(h(:)));
Lapx = lap1d(x);
y = x - P.tau * (Ht_r + P.alpha * Lapx);
y = clip01(y);
end

function val = obj_val_gaussian(x, b, h, P)
r = conv_circ(x, h) - b;
data = 0.5 * sum(r.^2);
reg  = 0.5 * P.alpha * sum((lap1d(x)).^2);
val  = data + reg;
end

% =========================================================================
% Helpers
% =========================================================================
function Lx = lap1d(x)
Lx = 2*x - circshift(x,1) - circshift(x,-1);
end

function y = conv_circ(x, h)
x = x(:);
h = h(:);
N = numel(x);
y = real(ifft(fft(x) .* fft(h, N)));
end

function k = gauss_kernel(w, s)
if mod(w,2) == 0
    w = w + 1;
end
m = (w-1)/2;
x = (-m:m)';
k = exp(-(x.^2)/(2*s^2));
k = k / sum(k);
end

function v = psnr_safe(x, y, maxI)
if nargin < 3
    maxI = 1;
end
x = double(x(:));
y = double(y(:));
mse = mean((x - y).^2);
if ~isfinite(mse) || mse <= eps
    v = 100;
else
    v = 10 * log10((maxI.^2) / mse);
end
end

function s = ssim1d(x, y, maxI)
if nargin < 3
    maxI = 1;
end

x = double(x(:));
y = double(y(:));

win = 11;
w = ones(win,1) / win;

mu_x = conv(x, w, 'same');
mu_y = conv(y, w, 'same');

x2 = x.^2;
y2 = y.^2;
xy = x .* y;

sigma_x2 = conv(x2, w, 'same') - mu_x.^2;
sigma_y2 = conv(y2, w, 'same') - mu_y.^2;
sigma_xy = conv(xy, w, 'same') - mu_x .* mu_y;

C1 = (0.01 * maxI)^2;
C2 = (0.03 * maxI)^2;

num = (2 * mu_x .* mu_y + C1) .* (2 * sigma_xy + C2);
den = (mu_x.^2 + mu_y.^2 + C1) .* (sigma_x2 + sigma_y2 + C2);

smap = num ./ (den + eps);

if numel(smap) > 2*(win-1)
    s = mean(smap(win:end-win+1));
else
    s = mean(smap);
end
end

function y = rescale01(x)
xmin = min(x(:));
xmax = max(x(:));
if xmax > xmin
    y = (x - xmin) / (xmax - xmin);
else
    y = zeros(size(x));
end
end

function y = clip01(x)
y = min(max(x,0),1);
end

function z = clamp(x, a, b)
z = min(max(x,a),b);
end

function s = make_synthetic_signal(N)
t = linspace(0,1,N)';

s = 0.4 + 0.3*sin(2*pi*3*t) + 0.2*exp(-150*(t-0.6).^2);
mask = (t > 0.7 & t < 0.8);
s(mask) = s(mask) + 0.3;

s = clip01(s);
end