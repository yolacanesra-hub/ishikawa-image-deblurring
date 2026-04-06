function ishikawa_main_selected_clean_plos_FINAL()

clc; clear; close all;
rng(0);

%% ================= PATH =================
inputFolder  = 'C:\Users\esray\Desktop\all_images';
outputFolder = 'C:\Users\esray\Desktop\comparison_figures';

if ~exist(inputFolder, 'dir')
    error('Girdi klasoru bulunamadi: %s', inputFolder);
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ================= IMAGE LIST =================
imageFiles = { ...
    '04.png', ...
    '06.png', ...
    '09.png', ...
    '10.png', ...
    '11.png', ...
    '12.png'};

imageLabels = { ...
    'Starfish', ...
    'Plane', ...
    'Woman', ...
    'Boats', ...
    'Pirate', ...
    'Couple'};

numImages = length(imageFiles);

fprintf('Toplam %d goruntu islenecek.\n', numImages);

%% ================= PARAMETERS =================
PSF = fspecial('motion', 15, 0);
noiseVar = 0.001;

rho = 0.6;
sigmaI = 0.3;
tau_data = 0.8;
tau_reg = 0.15;
lambda_reg = 0.01;
maxIter = 80;

tv_lambda = 0.005;
tv_rho    = 1.0;
tv_iter   = 150;

%% ================= STORAGE =================
PSNR_all = nan(numImages,5);
SSIM_all = nan(numImages,5);

GT_all       = cell(numImages,1);
Noisy_all    = cell(numImages,1);
Proposed_all = cell(numImages,1);
TV_all       = cell(numImages,1);
Wiener_all   = cell(numImages,1);
LR_all       = cell(numImages,1);
FISTA_all    = cell(numImages,1);

%% ================= MAIN LOOP =================
for n = 1:numImages

    fname = imageFiles{n};
    fpath = fullfile(inputFolder, fname);

    if ~exist(fpath, 'file')
        fprintf('Eksik dosya: %s\n', fname);
        continue;
    end

    I = im2double(imread(fpath));

    if size(I,3) == 3
        I = rgb2gray(I);
    end

    I = mat2gray(I);

    fprintf('Isleniyor: %s (%s)\n', fname, imageLabels{n});

    %% Degradation
    blurred = imfilter(I, PSF, 'circular', 'conv');
    noisy   = imnoise(blurred, 'gaussian', 0, noiseVar);

    %% Methods
    proposed = proposed_ishikawa_restore(noisy, PSF, ...
        rho, sigmaI, tau_data, tau_reg, lambda_reg, maxIter);

    wiener = deconvwnr(noisy, PSF, noiseVar);
    lr     = deconvlucy(noisy, PSF, 20);
    tv     = tv_admm_deblur(noisy, PSF, tv_lambda, tv_rho, tv_iter);
    fista  = fista_deblur(noisy, PSF, 0.005, 0.8, 300);

    %% Clamp
    proposed = min(max(proposed, 0), 1);
    wiener   = min(max(wiener,   0), 1);
    lr       = min(max(lr,       0), 1);
    tv       = min(max(tv,       0), 1);
    fista    = min(max(fista,    0), 1);

    %% Metrics
    PSNR_all(n,:) = [ ...
        compute_psnr_manual(proposed, I), ...
        compute_psnr_manual(wiener,   I), ...
        compute_psnr_manual(lr,       I), ...
        compute_psnr_manual(tv,       I), ...
        compute_psnr_manual(fista,    I)];

    SSIM_all(n,:) = [ ...
        compute_ssim_safe(proposed, I), ...
        compute_ssim_safe(wiener,   I), ...
        compute_ssim_safe(lr,       I), ...
        compute_ssim_safe(tv,       I), ...
        compute_ssim_safe(fista,    I)];

    %% Save outputs for combined figures
    GT_all{n}       = I;
    Noisy_all{n}    = noisy;
    Proposed_all{n} = proposed;
    TV_all{n}       = tv;
    Wiener_all{n}   = wiener;
    LR_all{n}       = lr;
    FISTA_all{n}    = fista;
end

%% ================= TABLES =================
validRows = ~isnan(PSNR_all(:,1));

T_PSNR = table( ...
    imageLabels(validRows)', ...
    PSNR_all(validRows,1), ...
    PSNR_all(validRows,2), ...
    PSNR_all(validRows,3), ...
    PSNR_all(validRows,4), ...
    PSNR_all(validRows,5), ...
    'VariableNames', {'Image','Proposed','Wiener','LR','TV','FISTA'});

writetable(T_PSNR, fullfile(outputFolder, 'PSNR_Table.xlsx'));

T_SSIM = table( ...
    imageLabels(validRows)', ...
    SSIM_all(validRows,1), ...
    SSIM_all(validRows,2), ...
    SSIM_all(validRows,3), ...
    SSIM_all(validRows,4), ...
    SSIM_all(validRows,5), ...
    'VariableNames', {'Image','Proposed','Wiener','LR','TV','FISTA'});

writetable(T_SSIM, fullfile(outputFolder, 'SSIM_Table.xlsx'));

%% ================= BAR GRAPHS =================
hBar1 = figure('Color','w','Position',[100 100 900 500]);
hb1 = bar(PSNR_all(validRows,:), 'LineWidth', 1.0);
set(gca, 'XTick', 1:sum(validRows), ...
    'XTickLabel', imageLabels(validRows), ...
    'XTickLabelRotation', 45, ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'FontName', 'Times New Roman');
legend(hb1, {'Proposed','Wiener','LR','TV','FISTA'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'FontName','Times New Roman');
title('PSNR Comparison','FontSize',13,'FontWeight','bold','FontName','Times New Roman');
ylabel('dB','FontSize',12,'FontName','Times New Roman');
grid on; box on;
set(hBar1,'PaperPositionMode','auto');
print(hBar1, fullfile(outputFolder, 'PSNR_BarGraph_300dpi.png'), '-dpng', '-r300');
close(hBar1);

hBar2 = figure('Color','w','Position',[100 100 900 500]);
hb2 = bar(SSIM_all(validRows,:), 'LineWidth', 1.0);
set(gca, 'XTick', 1:sum(validRows), ...
    'XTickLabel', imageLabels(validRows), ...
    'XTickLabelRotation', 45, ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'FontName', 'Times New Roman');
legend(hb2, {'Proposed','Wiener','LR','TV','FISTA'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'FontName','Times New Roman');
title('SSIM Comparison','FontSize',13,'FontWeight','bold','FontName','Times New Roman');
ylabel('SSIM','FontSize',12,'FontName','Times New Roman');
grid on; box on;
set(hBar2,'PaperPositionMode','auto');
print(hBar2, fullfile(outputFolder, 'SSIM_BarGraph_300dpi.png'), '-dpng', '-r300');
close(hBar2);

%% ================= FIGURE 2 =================
selectedIdx = [2, 1, 4]; % Plane, Starfish, Boats
rowTags = {'(A) Plane', '(B) Starfish', '(C) Boats'};

zoomRects = {
    [120 90 60 60], ...
    [90 110 60 60], ...
    [110 80 60 60]};

savePathFig2 = fullfile(outputFolder, 'Figure2_NEW.png');

create_figure2_publishable( ...
    GT_all, Noisy_all, Proposed_all, TV_all, FISTA_all, LR_all, Wiener_all, ...
    PSNR_all, SSIM_all, selectedIdx, rowTags, zoomRects, savePathFig2);

fprintf('Yeni Figure 2 kaydedildi: %s\n', savePathFig2);

%% ================= APPENDIX FIGURE A2 =================
selectedIdx_A2 = [3, 5, 6]; % Woman, Pirate, Couple
rowTags_A2 = {'(A) Woman', '(B) Pirate', '(C) Couple'};

zoomRects_A2 = { ...
    [90 80 60 60], ...
    [100 90 60 60], ...
    [110 85 60 60]};

savePathA2 = fullfile(outputFolder, 'Appendix_Figure_A2.png');

create_appendix_figureA2( ...
    GT_all, Noisy_all, Proposed_all, TV_all, FISTA_all, LR_all, Wiener_all, ...
    PSNR_all, SSIM_all, selectedIdx_A2, rowTags_A2, zoomRects_A2, savePathA2);

fprintf('Appendix Figure A2 kaydedildi: %s\n', savePathA2);

%% ================= FIGURE 3 (3D SENSITIVITY) =================
rhoVals   = [0.2 0.4 0.6 0.8 1.0];
sigmaVals = [0.1 0.2 0.3 0.4 0.5];

planeIdx = 2; % Plane
I_plane  = GT_all{planeIdx};

blurred_plane = imfilter(I_plane, PSF, 'circular', 'conv');
noisy_plane   = imnoise(blurred_plane, 'gaussian', 0, noiseVar);

sensPSNR = zeros(length(rhoVals), length(sigmaVals));
sensSSIM = zeros(length(rhoVals), length(sigmaVals));

for i = 1:length(rhoVals)
    for j = 1:length(sigmaVals)
        out = proposed_ishikawa_restore(noisy_plane, PSF, ...
            rhoVals(i), sigmaVals(j), tau_data, tau_reg, lambda_reg, maxIter);

        out = min(max(out,0),1);

        sensPSNR(i,j) = compute_psnr_manual(out, I_plane);
        sensSSIM(i,j) = compute_ssim_safe(out, I_plane);
    end
end

[SIG, RHO] = meshgrid(sigmaVals, rhoVals);

hFig3 = figure('Color','w','Position',[100 100 1400 600]);

subplot(1,2,1);
surf(SIG, RHO, sensPSNR);
shading interp;
colormap(parula);
xlabel('\sigma_I','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
ylabel('\rho','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
zlabel('PSNR (dB)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
title('(A) PSNR','FontSize',13,'FontWeight','bold','FontName','Times New Roman');
set(gca,'FontSize',11,'LineWidth',1,'FontName','Times New Roman');
grid on; box on;
view(135,30);
cb = colorbar;
ylabel(cb, 'PSNR (dB)', 'FontSize', 11, 'FontName', 'Times New Roman');
hold on;
zChosen1 = compute_psnr_manual( ...
    min(max(proposed_ishikawa_restore(noisy_plane, PSF, ...
    rho, sigmaI, tau_data, tau_reg, lambda_reg, maxIter),0),1), I_plane);
plot3(sigmaI, rho, zChosen1, 'kp', 'MarkerSize', 10, ...
      'MarkerFaceColor', 'y', 'LineWidth', 1.5);
hold off;

subplot(1,2,2);
surf(SIG, RHO, sensSSIM);
shading interp;
colormap(parula);
xlabel('\sigma_I','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
ylabel('\rho','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
zlabel('SSIM','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
title('(B) SSIM','FontSize',13,'FontWeight','bold','FontName','Times New Roman');
set(gca,'FontSize',11,'LineWidth',1,'FontName','Times New Roman');
grid on; box on;
view(135,30);
cb = colorbar;
ylabel(cb, 'SSIM', 'FontSize', 11, 'FontName', 'Times New Roman');
hold on;
zChosen2 = compute_ssim_safe( ...
    min(max(proposed_ishikawa_restore(noisy_plane, PSF, ...
    rho, sigmaI, tau_data, tau_reg, lambda_reg, maxIter),0),1), I_plane);
plot3(sigmaI, rho, zChosen2, 'kp', 'MarkerSize', 10, ...
      'MarkerFaceColor', 'y', 'LineWidth', 1.5);
hold off;

set(hFig3,'PaperPositionMode','auto');
print(hFig3, fullfile(outputFolder, 'Figure3.png'), '-dpng', '-r300');
close(hFig3);

%% ================= SENSITIVITY TABLES =================
for n = 1:numImages

    fname = imageFiles{n};
    fpath = fullfile(inputFolder, fname);

    if ~exist(fpath, 'file')
        fprintf('Sensitivity icin atlandi: %s\n', fname);
        continue;
    end

    I_sens = im2double(imread(fpath));
    if size(I_sens,3) == 3
        I_sens = rgb2gray(I_sens);
    end
    I_sens = mat2gray(I_sens);

    blurred_sens = imfilter(I_sens, PSF, 'circular', 'conv');
    noisy_sens   = imnoise(blurred_sens, 'gaussian', 0, noiseVar);

    sensPSNR_local = zeros(length(rhoVals), length(sigmaVals));
    sensSSIM_local = zeros(length(rhoVals), length(sigmaVals));

    for i = 1:length(rhoVals)
        for j = 1:length(sigmaVals)
            out = proposed_ishikawa_restore(noisy_sens, PSF, ...
                rhoVals(i), sigmaVals(j), tau_data, tau_reg, lambda_reg, maxIter);

            out = min(max(out,0),1);

            sensPSNR_local(i,j) = compute_psnr_manual(out, I_sens);
            sensSSIM_local(i,j) = compute_ssim_safe(out, I_sens);
        end
    end

    rowNames = cell(length(rhoVals),1);
    for i = 1:length(rhoVals)
        rowNames{i} = sprintf('rho_%d', round(rhoVals(i)*10));
    end

    colNames = cell(1,length(sigmaVals));
    for j = 1:length(sigmaVals)
        colNames{j} = sprintf('sigma_%d', round(sigmaVals(j)*10));
    end

    T_psnr = array2table(sensPSNR_local, 'VariableNames', colNames, 'RowNames', rowNames);
    T_ssim = array2table(sensSSIM_local, 'VariableNames', colNames, 'RowNames', rowNames);

    label = imageLabels{n};

    writetable(T_psnr, fullfile(outputFolder, sprintf('Sensitivity_PSNR_%s.xlsx', label)), ...
        'WriteRowNames', true);
    writetable(T_ssim, fullfile(outputFolder, sprintf('Sensitivity_SSIM_%s.xlsx', label)), ...
        'WriteRowNames', true);
end

fprintf('\nTum islemler tamamlandi.\n');
fprintf('Figure 2, Figure 3 ve Appendix figürleri son haliyle olusturuldu.\n');
fprintf('Cikti klasoru: %s\n', outputFolder);

end