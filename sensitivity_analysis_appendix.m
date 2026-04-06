function sensitivity_analysis_appendix()

clc; clear; close all;
rng(0);

%% ================= PATH =================
inputFolder  = 'C:\Users\esray\Desktop\all_images';
outputFolder = 'C:\Users\esray\Desktop\comparison_figures';

imageFiles  = {'04.png','10.png'}; % Starfish, Boats
imageLabels = {'Starfish','Boats'};

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ================= PARAMETERS =================
PSF = fspecial('motion', 15, 0);
noiseVar = 0.001;

rhoVals   = [0.1 0.2 0.4 0.6 0.8 1.0];
sigmaVals = [0.1 0.2 0.3 0.4 0.5 0.6];

tau_data   = 0.8;
tau_reg    = 0.15;
lambda_reg = 0.01;
maxIter    = 80;

rho_fixed   = 0.6;
sigma_fixed = 0.3;

%% ================= FIGURE =================
hFig = figure('Color','w','Position',[100 100 1400 700]);

plotIdx = 1;

for n = 1:2

    %% ================= LOAD IMAGE =================
    I = im2double(imread(fullfile(inputFolder, imageFiles{n})));
    if size(I,3) == 3
        I = rgb2gray(I);
    end
    I = mat2gray(I);

    %% ================= DEGRADATION =================
    blurred = imfilter(I, PSF, 'circular', 'conv');
    noisy   = imnoise(blurred, 'gaussian', 0, noiseVar);

    %% ================= RHO ANALYSIS =================
    PSNR_rho = zeros(size(rhoVals));
    SSIM_rho = zeros(size(rhoVals));

    for i = 1:length(rhoVals)
        out = proposed_ishikawa_restore(noisy, PSF, ...
            rhoVals(i), sigma_fixed, tau_data, tau_reg, lambda_reg, maxIter);

        out = min(max(out,0),1);

        PSNR_rho(i) = compute_psnr_manual(out, I);
        SSIM_rho(i) = compute_ssim_safe(out, I);
    end

    subplot(2,2,plotIdx);
    [ax,h1,h2] = plotyy(rhoVals, PSNR_rho, rhoVals, SSIM_rho);
    plotIdx = plotIdx + 1;

    set(h1, 'LineStyle','-',  'Marker','o', 'LineWidth',1.8, 'MarkerSize',7);
    set(h2, 'LineStyle','--', 'Marker','s', 'LineWidth',1.8, 'MarkerSize',7);

    % Colors
    set(h1, 'Color', [0 0.4470 0.7410]);         % blue for PSNR
    set(h2, 'Color', [0.8500 0.3250 0.0980]);    % orange for SSIM

    % Axis labels
    xlabel('\rho', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel(ax(1), 'PSNR (dB)', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ylabel(ax(2), 'SSIM', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    % Title
    title([imageLabels{n}, ' - \rho'], ...
        'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    % Axis style
    set(ax(1), 'FontSize', 11, 'LineWidth', 0.8, ...
        'FontName', 'Times New Roman', ...
        'Color', 'w', ...
        'YColor', [0 0.4470 0.7410]);
    set(ax(2), 'FontSize', 11, 'LineWidth', 0.8, ...
        'FontName', 'Times New Roman', ...
        'Color', 'none', ...
        'YColor', [0.8500 0.3250 0.0980]);

    grid(ax(1), 'on');
    box(ax(1), 'on');

    % Axis limits to match paper-style appearance
    if strcmp(imageLabels{n}, 'Starfish')
        ylim(ax(1), [20.5 22.5]);
        ylim(ax(2), [0.56 0.64]);
    elseif strcmp(imageLabels{n}, 'Boats')
        ylim(ax(1), [23.4 24.4]);
        ylim(ax(2), [0.52 0.62]);
    end

    legend([h1 h2], {'PSNR','SSIM'}, ...
        'Location', 'northwest', ...
        'FontSize', 10, ...
        'FontName', 'Times New Roman');

    %% ================= SIGMA ANALYSIS =================
    PSNR_sigma = zeros(size(sigmaVals));
    SSIM_sigma = zeros(size(sigmaVals));

    for i = 1:length(sigmaVals)
        out = proposed_ishikawa_restore(noisy, PSF, ...
            rho_fixed, sigmaVals(i), tau_data, tau_reg, lambda_reg, maxIter);

        out = min(max(out,0),1);

        PSNR_sigma(i) = compute_psnr_manual(out, I);
        SSIM_sigma(i) = compute_ssim_safe(out, I);
    end

    subplot(2,2,plotIdx);
    [ax,h3,h4] = plotyy(sigmaVals, PSNR_sigma, sigmaVals, SSIM_sigma);
    plotIdx = plotIdx + 1;

    set(h3, 'LineStyle','-',  'Marker','o', 'LineWidth',1.8, 'MarkerSize',7);
    set(h4, 'LineStyle','--', 'Marker','s', 'LineWidth',1.8, 'MarkerSize',7);

    % Colors
    set(h3, 'Color', [0 0.4470 0.7410]);         % blue for PSNR
    set(h4, 'Color', [0.8500 0.3250 0.0980]);    % orange for SSIM

    % Axis labels
    xlabel('\sigma_I', 'FontSize', 12, 'FontName', 'Times New Roman');
    ylabel(ax(1), 'PSNR (dB)', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    ylabel(ax(2), 'SSIM', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    % Title
    title([imageLabels{n}, ' - \sigma_I'], ...
        'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    % Axis style
    set(ax(1), 'FontSize', 11, 'LineWidth', 0.8, ...
        'FontName', 'Times New Roman', ...
        'Color', 'w', ...
        'YColor', [0 0.4470 0.7410]);
    set(ax(2), 'FontSize', 11, 'LineWidth', 0.8, ...
        'FontName', 'Times New Roman', ...
        'Color', 'none', ...
        'YColor', [0.8500 0.3250 0.0980]);

    grid(ax(1), 'on');
    box(ax(1), 'on');

    if strcmp(imageLabels{n}, 'Starfish')
        ylim(ax(1), [21.0 23.0]);
        ylim(ax(2), [0.50 0.70]);
    elseif strcmp(imageLabels{n}, 'Boats')
        ylim(ax(1), [23.2 24.6]);
        ylim(ax(2), [0.53 0.60]);
    end

    legend([h3 h4], {'PSNR','SSIM'}, ...
        'Location', 'northwest', ...
        'FontSize', 10, ...
        'FontName', 'Times New Roman');
end

%% ================= FIGURE TITLE =================
annotation('textbox', [0 0.94 1 0.04], ...
    'String', 'Fig. S4. Sensitivity analysis for additional test images.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

%% ================= SAVE =================
set(hFig, 'PaperPositionMode', 'auto');

print(hFig, fullfile(outputFolder, 'Figure_S4.tif'), '-dtiff', '-r600');
print(hFig, fullfile(outputFolder, 'Figure_S4.png'), '-dpng', '-r300');

fprintf('\nFigure kaydedildi:\n');
fprintf('%s\n', fullfile(outputFolder, 'Figure_S4.tif'));
fprintf('%s\n', fullfile(outputFolder, 'Figure_S4.png'));

end