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

% Include other functions here as previously defined...

% To run this function from the command line:
% matlab -batch "main('path/to/image.jpg', 'quality')"
