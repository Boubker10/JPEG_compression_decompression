function main(varargin)
    % Parse input arguments
    p = inputParser;
    addOptional(p, 'imagePath', 'chadli.jpg', @ischar);
    addOptional(p, 'quality', 'high', @ischar);
    parse(p, varargin{:});

    imagePath = p.Results.imagePath;
    quality = p.Results.quality;

    % Load image
    img = imread(imagePath);

    % Process image
    fullJPEGProcess(img, quality);
end

function fullJPEGProcess(img, quality)
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    img = double(img);

    % Calcul de la table de quantification et compression/décompression
    Q = getQuantizationTable(quality);
    [compressedImg, decompressedImg] = compressDecompress(img, Q);

    decompressedImageToShow = uint8(decompressedImg + 128);

    % Calcul de PSNR et SSIM
    psnrValue = calculatePSNR(uint8(img), decompressedImageToShow);
    ssimValue = calculateSSIM(uint8(img), decompressedImageToShow);

    % Estimation de la taille de l'image compressée
    nonZeroCoeffs = nnz(compressedImg);
    estimatedCompressedSize = nonZeroCoeffs * 8;
    originalSize = numel(img) * 8;
    reductionPercentage = 100 * (1 - estimatedCompressedSize / originalSize);

    % Préparation des textes à afficher
    textStringPSNR = sprintf('PSNR: %f dB', psnrValue);
    textStringSSIM = sprintf('SSIM: %f', ssimValue);
    textStringSize = sprintf('Original: %d bits, Compressed: %d bits', originalSize, estimatedCompressedSize);
    textStringRe = sprintf('Reduction: %f%%', reductionPercentage);

    % Affichage des images et informations dans une seule fenêtre
    figure;
    subplot(1, 2, 1), imshow(uint8(img)), title('Image Originale');
    subplot(1, 2, 2), imshow(decompressedImageToShow), title(['Image Après (', quality, ' qualité)']);

    % Utilisation d'annotation pour afficher les informations textuelles
    dim = [.2 .1 .3 .3];
    str = {textStringPSNR, textStringSSIM, textStringSize, textStringRe};
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
end

function [compressedImg, decompressedImg] = compressDecompress(img, Q)
    imgCentered = img - 128;
    [rows, cols] = size(img);
    compressedImg = zeros(size(img));
    decompressedImg = zeros(size(img));

    for i = 1:8:rows-7
        for j = 1:8:cols-7
            block = imgCentered(i:i+7, j:j+7);
            dctBlock = dct2_custom(block);
            
            quantizedBlock = round(dctBlock ./ Q);
            compressedImg(i:i+7, j:j+7) = quantizedBlock;

            dequantizedBlock = quantizedBlock .* Q;
            decompressedBlock = idct2_custom(dequantizedBlock);
            decompressedImg(i:i+7, j:j+7) = decompressedBlock;
        end
    end
end

function Q = getQuantizationTable(quality)
    Q_high = [16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 14 13 16 24 40 57 69 56; 14 17 22 29 51 87 80 62; 18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 49 64 78 87 103 121 120 101; 72 92 95 98 112 100 103 99];

    % Ajuste les coefficients pour un contraste plus marqué en PSNR
    switch quality
        case 'high'
            scaleFactor = 0.4; % Pour une compression minimale et une qualité maximale
        case 'medium'
            scaleFactor = 7.5; % Un équilibre entre qualité et compression
        case 'low'
            scaleFactor = 20; % Maximise la compression pour minimiser la taille du fichier
        otherwise
            scaleFactor = 1.5; % Par défaut à la qualité moyenne si l'option n'est pas reconnue
    end
    Q = Q_high * scaleFactor;
end

function dctBlock = dct2_custom(block)
    N = 8;
    dctBlock = zeros(N, N);
    for u = 0:N-1
        for v = 0:N-1
            sum = 0;
            for x = 0:N-1
                for y = 0:N-1
                    sum = sum + block(x+1, y+1) * cos((2*x+1)*u*pi/(2*N)) * cos((2*y+1)*v*pi/(2*N));
                end
            end
            cu = 1; cv = 1;
            if u == 0, cu = 1/sqrt(2); end
            if v == 0, cv = 1/sqrt(2); end
            dctBlock(u+1, v+1) = 0.25 * cu * cv * sum;
        end
    end
end

function block = idct2_custom(dctBlock)
    N = 8;
    block = zeros(N, N);
    for x = 0:N-1
        for y = 0:N-1
            sum = 0;
            for u = 0:N-1
                for v = 0:N-1
                    cu = 1; cv = 1;
                    if u == 0, cu = 1/sqrt(2); end
                    if v == 0, cv = 1/sqrt(2); end
                    sum = sum + cu * cv * dctBlock(u+1, v+1) * cos((2*x+1)*u*pi/(2*N)) * cos((2*y+1)*v*pi/(2*N));
                end
            end
            block(x+1, y+1) = sum * 0.25;
        end
    end
end

function psnrValue = calculatePSNR(original, compressed)
    mse = mean((original(:) - compressed(:)) .^ 2);
    if mse == 0
        psnrValue = 100;
    else
        maxPixel = 255.0;
        psnrValue = 10 * log10((maxPixel^2) / mse);
    end
end

function window = createGaussianWindow(size, sigma)
    [x, y] = meshgrid(-floor(size/2):floor(size/2), -floor(size/2):floor(size/2));
    window = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    window = window / sum(window(:));
end

function ssimValue = calculateSSIM(img1, img2)
    K = [0.01, 0.03];
    L = 255;
    size = 11;
    sigma = 1.5; 
    window = createGaussianWindow(size, sigma); 
    C1 = (K(1) * L)^2;
    C2 = (K(2) * L)^2;
    img1 = double(img1);
    img2 = double(img2);
    
    mu1 = filter2(window, img1, 'valid');
    mu2 = filter2(window, img2, 'valid');
    mu1_sq = mu1 .* mu1;
    mu2_sq = mu2 .* mu2;
    mu1_mu2 = mu1 .* mu2;
    
    sigma1_sq = filter2(window, img1 .* img1, 'valid') - mu1_sq;
    sigma2_sq = filter2(window, img2 .* img2, 'valid') - mu2_sq;
    sigma12 = filter2(window, img1 .* img2, 'valid') - mu1_mu2;
    
    ssim_map = ((2 * mu1_mu2 + C1) .* (2 * sigma12 + C2)) ./ ((mu1_sq + mu2_sq + C1) .* (sigma1_sq + sigma2_sq + C2));
    ssimValue = mean(ssim_map(:));
end

% Pour exécuter cette fonction à partir de la ligne de commande :
% matlab -batch "main('path/to/image.jpg', 'quality')"

