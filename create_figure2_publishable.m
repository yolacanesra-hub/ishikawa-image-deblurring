function create_figure2_publishable( ...
    GT_all, Noisy_all, Proposed_all, TV_all, FISTA_all, LR_all, Wiener_all, ...
    PSNR_all, SSIM_all, selectedIdx, rowTags, zoomRects, savePath)
% MATLAB 2016a compatible publishable Figure 2
% Layout:
% Row block 1: Full images
% Row block 2: Zoom patches
% Row block 3: Error maps
%
% Columns:
% GT | Blurred | Proposed | TV | FISTA | LR | Wiener

colTitles = {'GT', 'Blurred', 'Proposed', 'TV', 'FISTA', 'LR', 'Wiener'};
numRows = length(selectedIdx);
numCols = 7;

hFig = figure('Color','w','Position',[30 30 2200 1400]);

lastErrAx = [];

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

    diffMaps = { ...
        zeros(size(GT)), ...
        abs(BLR  - GT), ...
        abs(PROP - GT), ...
        abs(TVI  - GT), ...
        abs(FIS  - GT), ...
        abs(LRI  - GT), ...
        abs(WIE  - GT)};

    % ---------------- FULL IMAGES ----------------
    for c = 1:numCols
        ax = subplot(numRows*3, numCols, (r-1)*numCols + c);

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
            [ps, ss] = get_metric_pair(c, n, PSNR_all, SSIM_all);
            text(6, size(imgs{c},1)-10, ...
                sprintf('PSNR: %.2f dB\nSSIM: %.3f', ps, ss), ...
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
                [pos(1)-0.040, pos(2)+pos(4)/2-0.015, 0.035, 0.03], ...
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
        ax = subplot(numRows*3, numCols, numRows*numCols + (r-1)*numCols + c);

        patch = crop_patch(imgs{c}, rect);
        imshow(patch, [], 'Parent', ax);
        axis(ax, 'off');

        if c == 1
            pos = get(ax, 'Position');
            annotation('textbox', ...
                [pos(1)-0.030, pos(2)+pos(4)/2-0.012, 0.025, 0.025], ...
                'String', 'Zoom', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman');
        end
    end

    % ---------------- ERROR MAPS ----------------
    for c = 1:numCols
        ax = subplot(numRows*3, numCols, 2*numRows*numCols + (r-1)*numCols + c);

        if c == 1
            imshow(zeros(size(GT)), [], 'Parent', ax);
            axis(ax, 'off');
            text(size(GT,2)/2, size(GT,1)/2, 'N/A', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman', ...
                'Parent', ax);
        else
            imagesc(diffMaps{c}, 'Parent', ax);
            axis(ax, 'image');
            axis(ax, 'off');
            colormap(ax, gray);
            caxis(ax, [0 0.4]);
            lastErrAx = ax;
        end

        if c == 1
            pos = get(ax, 'Position');
            annotation('textbox', ...
                [pos(1)-0.035, pos(2)+pos(4)/2-0.012, 0.03, 0.025], ...
                'String', 'Error Map', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'FontName', 'Times New Roman');
        end
    end
end

% ---------------- COLORBAR ----------------
if ~isempty(lastErrAx)
    cb = colorbar('peer', lastErrAx);
    set(cb, 'Units', 'normalized');
    set(cb, 'Position', [0.92 0.12 0.012 0.25]);
    set(cb, 'FontSize', 10, 'FontName', 'Times New Roman');
    ylabel(cb, '|Error|', 'FontName', 'Times New Roman');
end

annotation('textbox', [0 0.965 1 0.03], ...
    'String', 'Fig. 2. Comparative restoration results, zoomed details, and error maps.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 13, ...
    'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

set(hFig, 'PaperPositionMode', 'auto');
print(hFig, savePath, '-dpng', '-r400');
close(hFig);

end

function patch = crop_patch(I, rect)
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

function [ps, ss] = get_metric_pair(c, n, PSNR_all, SSIM_all)
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
        ps = NaN; ss = NaN;
end
end