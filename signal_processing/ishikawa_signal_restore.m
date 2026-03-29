function R = gi_signal_enhancement_v4(noise_type, opts)
% GI_SIGNAL_ENHANCEMENT_V4 (MATLAB R2016 compatible)
% Ishikawa-tipi 2-aşamalı şema + robust veri terimleri (Gaussian: L2,
% S&P: Huber), güvenli PSNR/SSIM, PSNR-tabanlı erken durdurma ve
% basit backtracking adım ayarı. Gaussian ve S&P için önerilen
% "son ayar"lar entegre edilmiştir.
%
% Kullanım:
%   R = gi_signal_enhancement_v4('gaussian');
%   R = gi_signal_enhancement_v4('sap');         % salt & pepper
%   R = gi_signal_enhancement_v4('gaussian', struct('tau',0.05)); % override

if nargin<1, noise_type = 'gaussian'; end
if nargin<2, opts = struct(); end

rng(0);   % <-- BURAYA
noise_type = lower(strtrim(noise_type));
assert(ismember(noise_type, {'gaussian','sap'}), 'noise_type must be gaussian or sap');

% -------------------------
% 0) Sentetik test sinyali
% -------------------------
N = 1024;
x_true = make_synthetic_signal(N);           % [0,1]
h = gauss_kernel(9, 2); h = h(:);            % normalize blur
Hx = conv_circ(x_true, h);

switch noise_type
    case 'gaussian'
        sigma_n = 0.01;
        b = clip01(Hx + sigma_n*randn(size(Hx)));
    case 'sap'
        d = 0.10; b = add_salt_pepper(Hx, d);
end

% -------------------------
% 1) Parametreler (önerilerin son hali)
% -------------------------
P.noise_type    = noise_type;
P.maxit         = 500;
P.psnr_maxI     = 1;
P.show_plots    = true;
P.stop_eps      = 1e-4;

% --- Ortak karışım parametreleri (daha muhafazakâr) ---
P.beta_a   = 8;  P.beta_b   = 60;  P.beta_min   = 0.03; P.beta_max   = 0.15;
P.sigma_a  = 8;  P.sigma_b  = 60;  P.sigma_min  = 0.03; P.sigma_max  = 0.15;
P.stop_patience = 8;   % PSNR iyileşmesi bekleme penceresi

% --- Gürültü modeline özel parametreler ---
switch noise_type
    case 'gaussian'    % (öneri: daha küçük adım + daha güçlü düzenleme)
        P.tau         = 0.06;
        P.alpha       = 0.08;
        P.huber_delta = 0.05;  % kullanılmaz ama tutarlı kalsın
        P.med_every   = 10;    % kullanılmaz
        P.med_win     = 3;
    case 'sap'         % (öneri: daha robust)
        P.tau         = 0.06;
        P.alpha       = 0.06;
        P.huber_delta = 0.03;
        P.med_every   = 5;     % daha sık medyan
        P.med_win     = 5;
end

% --- Kullanıcı overrides ---
fn = fieldnames(opts);
for k=1:numel(fn), P.(fn{k}) = opts.(fn{k}); end

% güvenli ölçek
x_true = rescale01(double(x_true));
b      = rescale01(double(b));

% -------------------------
% 2) Ayırmalar
% -------------------------
x = b;  % başlangıç
PSNR = nan(P.maxit,1); SSIM = nan(P.maxit,1); J = nan(P.maxit,1);
bestPSNR = -Inf; bestx = x; best_iter = 0;
nan_or_inf = @(z) any(~isfinite(z(:)));

% -------------------------
% 3) İterasyon döngüsü
% -------------------------
for n = 1:P.maxit
    beta_n  = clamp(P.beta_a /(n + P.beta_b),  P.beta_min,  P.beta_max);
    sigma_n = clamp(P.sigma_a/(n + P.sigma_b), P.sigma_min, P.sigma_max);

    % Operatör C
    Cx = grad_step(x, b, h, P);
    y  = (1 - beta_n)*x + beta_n*Cx;

    % S&P için ara medyan
    if strcmp(P.noise_type,'sap') && mod(n,P.med_every)==0
        if exist('medfilt1','file'), y = medfilt1(y, P.med_win, 'truncate');
        else,                         y = local_median(y, P.med_win);
        end
    end

    % --------- Backtracking ile I-operatorü ----------
    if n==1 || ~isfinite(J(n-1)), J_prev = inf; else, J_prev = J(n-1); end
    local_tau = P.tau;  max_bt = 8;  bt = 0;
    while true
        P_tmp = P; P_tmp.tau = local_tau;
        Iy = grad_step(y, b, h, P_tmp);
        x_cand = (1 - sigma_n)*x + sigma_n*Iy;
        J_cand = obj_val(x_cand, b, h, P);
        if ~isfinite(J_cand), J_cand = inf; end
        if ~isfinite(J_prev) || (J_cand <= J_prev*(1+1e-5)) || bt>=max_bt
            x_new = x_cand; J(n) = J_cand; break;
        else
            local_tau = local_tau * 0.5; bt = bt + 1;
        end
    end
    % --------------------------------------------------

    % Metrikler
    PSNR(n) = psnr_safe(x_new, x_true, P.psnr_maxI);
    SSIM(n) = ssim1d(x_new, x_true, P.psnr_maxI);

    % Koruma
    if nan_or_inf([J(n), PSNR(n), SSIM(n)]) || J(n) > 1e12
        fprintf('[WARN] Divergence/overflow at iter %d. Stopping.\n', n);
        x = x_new; break;
    end

    % En iyi PSNR
    if PSNR(n) > bestPSNR
        bestPSNR = PSNR(n); bestx = x_new; best_iter = n;
    end

    % ---- PSNR ve hedef düzlüğü tabanlı erken durdurma ----
    if n > P.stop_patience
        if (n - best_iter) >= P.stop_patience
            fprintf('Early stop at iter %d (no PSNR improvement %d iters).\n', n, P.stop_patience);
            x = x_new; break;
        end
        Jw  = J(n-P.stop_patience+1:n);
        rel = abs(diff(Jw)) ./ max(abs(Jw(1:end-1)), eps);
        rel = rel(~isnan(rel) & ~isinf(rel));
        if ~isempty(rel) && mean(rel) < P.stop_eps
            fprintf('Early stop at iter %d (flat objective).\n', n);
            x = x_new; break;
        end
    end
    % -------------------------------------------------------

    x = x_new;
end

it_end = n; PSNR = PSNR(1:it_end); SSIM = SSIM(1:it_end); J = J(1:it_end);

% -------------------------
% 4) Grafikler
% -------------------------
if P.show_plots
    figure('Name',['GI Enhancement - ', upper(P.noise_type)], 'Color','w');
    subplot(3,1,1); plot(1:it_end, J, 'LineWidth',1.2); grid on;
    title('Objective Function'); xlabel('Iteration'); ylabel('Value');
    subplot(3,1,2); plot(1:it_end, PSNR, 'LineWidth',1.2); grid on;
    ylabel('PSNR (dB)'); xlabel('Iteration'); title('PSNR vs Iteration');
    subplot(3,1,3); plot(1:it_end, SSIM, 'LineWidth',1.2); grid on;
    ylabel('SSIM'); xlabel('Iteration'); title('SSIM vs Iteration');
end

% -------------------------
% 5) Çıktı
% -------------------------
R = struct();
R.x_true    = x_true;
R.b         = b;
R.x         = bestx;
R.PSNR      = PSNR;
R.SSIM      = SSIM;
R.J         = J;
R.best_iter = best_iter;
R.params    = P;
R.noise_type= noise_type;

fprintf('[%s] Ended at iter %d, Best PSNR: %.2f dB (iter %d), SSIM(final)=%.4f\n',...
    upper(noise_type), it_end, bestPSNR, best_iter, SSIM(end));
end

% =========================
% ===== Yardımcılar =======
% =========================
function y = grad_step(x, b, h, P)
r = conv_circ(x, h) - b;
switch P.noise_type
    case 'gaussian', psi = r;                       % L2
    case 'sap',      psi = huber_grad(r,P.huber_delta); % robust
end
Ht_psi = conv_circ(psi, flipud(h(:)));
Lapx   = lap1d(x);
y = x - P.tau * (Ht_psi + P.alpha * Lapx);
y = clip01(y);
end

function val = obj_val(x, b, h, P)
r = conv_circ(x, h) - b;
switch P.noise_type
    case 'gaussian', data = 0.5*sum(r.^2);
    case 'sap',      data = sum(huber(r, P.huber_delta));
end
reg = 0.5*P.alpha * sum( (lap1d(x)).^2 );
val = data + reg;
end

function Lx = lap1d(x), Lx = 2*x - circshift(x,1) - circshift(x,-1); end

function y = conv_circ(x, h)
x = x(:); h = h(:); N = numel(x);
y = real(ifft( fft(x).*fft(h, N) ));
end

function k = gauss_kernel(w, s)
if mod(w,2)==0, w = w+1; end
m = (w-1)/2; x = (-m:m)';
k = exp(-(x.^2)/(2*s^2)); k = k/sum(k);
end

function y = add_salt_pepper(x, d)
N = numel(x); y = x; m = rand(N,1);
y(m < d/2) = 0; y(m > 1-d/2) = 1;
end

function r = huber(z, delta)
az = abs(z); r = zeros(size(z)); q = az <= delta;
r(q) = 0.5*z(q).^2; r(~q) = delta*az(~q) - 0.5*delta^2;
end

function g = huber_grad(z, delta)
az = abs(z); g = zeros(size(z)); q = az <= delta;
g(q) = z(q); g(~q) = delta .* sign(z(~q));
end

function v = psnr_safe(x, y, maxI)
if nargin<3, maxI = 1; end
x = double(x(:)); y = double(y(:));
mse = mean((x - y).^2);
if ~isfinite(mse) || mse <= eps, v = 100;
else, v = 10*log10((maxI.^2)/mse);
end
end

function s = ssim1d(x, y, maxI)
if nargin<3, maxI = 1; end
x = double(x(:)); y = double(y(:));
win = 11; w = ones(win,1)/win;
mu_x = conv(x,w,'same'); mu_y = conv(y,w,'same');
x2 = x.^2; y2 = y.^2; xy = x.*y;
sigma_x2 = conv(x2,w,'same') - mu_x.^2;
sigma_y2 = conv(y2,w,'same') - mu_y.^2;
sigma_xy = conv(xy,w,'same') - mu_x.*mu_y;
C1 = (0.01*maxI)^2; C2 = (0.03*maxI)^2;
num = (2*mu_x.*mu_y + C1).*(2*sigma_xy + C2);
den = (mu_x.^2 + mu_y.^2 + C1).*(sigma_x2 + sigma_y2 + C2);
smap = num ./ (den + eps);
s = mean(smap(win:end-win+1));
end

function y = rescale01(x)
xmin = min(x(:)); xmax = max(x(:));
if xmax>xmin, y = (x - xmin)/(xmax - xmin); else, y = zeros(size(x)); end
end

function y = clip01(x), y = min(max(x,0),1); end
function z = clamp(x, a, b), z = min(max(x,a),b); end

function y = local_median(x, w)
if mod(w,2)==0, w = w+1; end
r = (w-1)/2; N = numel(x); y = x;
for i=1:N, L=max(1,i-r); R=min(N,i+r); y(i)=median(x(L:R)); end
end

function s = make_synthetic_signal(N)
t = linspace(0,1,N)';
s = 0.4 + 0.3*sin(2*pi*3*t) + 0.2*exp(-150*(t-0.6).^2);
s(t>0.7 & t<0.8) = s(t>0.7 & t<0.8) + 0.3;
s = clip01(s);
end
