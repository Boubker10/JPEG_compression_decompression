# JPEG Compression and Decompression in MATLAB

This project demonstrates a JPEG-like image compression and decompression process using MATLAB. The main function loads an image, processes it using Discrete Cosine Transform (DCT), applies quantization based on the specified quality, and then reconstructs the image. The decompressed image is saved to disk, and metrics such as PSNR and SSIM are calculated to evaluate the compression quality.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
  - [main](#main)
  - [fullJPEGProcess](#fulljpegprocess)
  - [saveImage](#saveimage)
  - [compressDecompress](#compressdecompress)
  - [getQuantizationTable](#getquantizationtable)
  - [dct2_custom](#dct2_custom)
  - [idct2_custom](#idct2_custom)
  - [calculatePSNR](#calculatepsnr)
  - [calculateSSIM](#calculatessim)
  - [createGaussianWindow](#creategaussianwindow)
- [Example](#example)
- [License](#license)

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/jpeg-compression-matlab.git
    ```
2. Open MATLAB and navigate to the cloned repository folder.

## Usage

To compress and decompress an image, run the following command in MATLAB:

```matlab
main('path/to/image.jpg', 'quality');
