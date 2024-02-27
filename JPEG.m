
image = imread('incredible0.png');


fullJPEGProcess(image, 'high');   % Pour une haute qualité
fullJPEGProcess(image, 'medium'); % Pour une qualité moyenne
fullJPEGProcess(image, 'low');    % Pour une faible qualité






function fullJPEGProcess(img, quality)
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    img = double(img);

    figure, imshow(uint8(img)), title('Image Originale');

    Q = getQuantizationTable(quality); % Obtenir la table de quantification avec qualité
    [compressedImg, decompressedImg] = compressDecompress(img, Q);

    figure, imshow(uint8(decompressedImg + 128)), title(['Image Après Compression et Décompression (', quality, ' qualité)']);
    

    psnrValue = calculatePSNR(uint8(img), uint8(decompressedImg + 128));
    fprintf('PSNR de l''image compressée (%s qualité): %f dB\n', quality, psnrValue);

    nonZeroCoeffs = nnz(compressedImg); % Nombre de coefficients non nuls dans l'image compressée
    estimatedCompressedSize = nonZeroCoeffs * 8; % Estimation de la taille en bits (8 bits par coefficient non nul)
    originalSize = numel(img) * 8; % Taille originale en bits (8 bits par pixel)

    fprintf('Taille originale de l''image: %d bits\n', originalSize);
    fprintf('Taille estimée de l''image compressée: %d bits\n', estimatedCompressedSize);
    fprintf('Réduction de taille: %f%%\n', 100 * (1-estimatedCompressedSize / originalSize));

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
            scaleFactor = 0.9; % Pour une compression minimale et une qualité maximale
        case 'medium'
            scaleFactor = 6.5; % Un équilibre entre qualité et compression
        case 'low'
            scaleFactor = 15; % Maximise la compression pour minimiser la taille du fichier
        otherwise
            scaleFactor = 3.5; % Par défaut à la qualité moyenne si l'option n'est pas reconnue
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
        psnrValue = 20 * log10(maxPixel / sqrt(mse));
    end
end

