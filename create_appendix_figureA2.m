function create_appendix_figureA2( ...
    GT_all, Noisy_all, Proposed_all, TV_all, FISTA_all, LR_all, Wiener_all, ...
    PSNR_all, SSIM_all, selectedIdx, rowTags, zoomRects, savePath)
% MATLAB 2016a compatible
% Appendix Figure A2
% Full images + zoom patches

colTitles = {'GT', 'Blurred', 'Proposed', 'TV', 'FISTA', 'LR', 'Wiener'};
numRows = length(selectedIdx);
numCols = 7;

hFig = figure('Color','w','Position',[50 50 2200 1100]);

for r = 1:numRows
    n = selectedIdx(r);
    rect = zoomRects{r};

    GT   = GT_all{n};
    BLR  = Noisy_all{n};
    PROP = Proposed_all{n};
    TVI  = TV_all{n};
    FIS  = FISTA_all{n};
    LRI  = LR_all{n};
    WIE  = Wiener_all{n};

    imgs = {GT, BLR, PROP, TVI, FIS, LRI, WIE};

    % ---------------- FULL IMAGES ----------------
    for c = 1:numCols
        ax = subplot(numRows*2, numCols, (r-1)*numCols + c);

        imshow(imgs{c}, [], 'Parent', ax);
        axis(ax, 'off');

        hold(ax, 'on');
        rectangle('Position', rect, 'EdgeColor', 'y', 'LineWidth', 1.8, 'Parent', ax);
        hold(ax, 'off');

        if r == 1
            title(colTitles{c}, ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman');
        end

        if c >= 3
            [ps, ss] = get_metric_pair_A2(c, n, PSNR_all, SSIM_all);

            text(6, size(imgs{c},1)-10, ...
                sprintf('PSNR: %.2f dB | SSIM: %.3f', ps, ss), ...
                'Color', 'k', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman', ...
                'BackgroundColor', 'w', ...
                'Margin', 2, ...
                'VerticalAlignment', 'bottom', ...
                'Parent', ax);
        end

        if c == 1
            pos = get(ax, 'Position');
            annotation('textbox', ...
                [pos(1)-0.035, pos(2)+pos(4)/2-0.015, 0.03, 0.03], ...
                'String', rowTags{r}, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 11, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman');
        end
    end

    % ---------------- ZOOM PATCHES ----------------
    for c = 1:numCols
        ax = subplot(numRows*2, numCols, numRows*numCols + (r-1)*numCols + c);

        patch = crop_patch_A2(imgs{c}, rect);
        imshow(patch, [], 'Parent', ax);
        axis(ax, 'off');

        if c == 1
            pos = get(ax, 'Position');
            annotation('textbox', ...
                [pos(1)-0.03, pos(2)+pos(4)/2-0.012, 0.025, 0.025], ...
                'String', 'Zoom', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman');
        end
    end
end

annotation('textbox', [0 0.955 1 0.035], ...
    'String', 'Fig. S3. Additional qualitative comparisons for the remaining test images.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 13, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

set(hFig, 'PaperPositionMode', 'auto');
print(hFig, savePath, '-dpng', '-r400');
close(hFig);

end

function patch = crop_patch_A2(I, rect)
x = round(rect(1));
y = round(rect(2));
w = round(rect(3));
h = round(rect(4));

[H, W] = size(I);

x1 = max(1, x);
y1 = max(1, y);
x2 = min(W, x + w - 1);
y2 = min(H, y + h - 1);

patch = I(y1:y2, x1:x2);
end

function [ps, ss] = get_metric_pair_A2(c, n, PSNR_all, SSIM_all)
switch c
    case 3
        ps = PSNR_all(n,1); ss = SSIM_all(n,1); % Proposed
    case 4
        ps = PSNR_all(n,4); ss = SSIM_all(n,4); % TV
    case 5
        ps = PSNR_all(n,5); ss = SSIM_all(n,5); % FISTA
    case 6
        ps = PSNR_all(n,3); ss = SSIM_all(n,3); % LR
    case 7
        ps = PSNR_all(n,2); ss = SSIM_all(n,2); % Wiener
    otherwise
        ps = NaN;
        ss = NaN;
end
end