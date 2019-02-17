%% shot noise for Nikon D5100
%% parameters: brightness_level, QE_curve
clc
clear 
close all

%% load data
load('camss_nikonD5100_IEmeas.mat');
wavelength = camss_NikonD5100(5:65, 1); % 400:5:700
delta_wavelength = (wavelength(end)- wavelength(1))/(length(wavelength)-1); % 5nm
rgb = camss_NikonD5100(5:65, 2:4); % sensor sensitivity

load('MCCCreflectance_380to780_5.mat');
reflectance = MCCCreflectance(5:65, 2:25);
load_cie_illuminants; % load illuminants
light_spd = IllA(21:5:321, 2);
% light_spd = D65(21:5:321, 2);
light_spd = light_spd/max(light_spd); % normalize

% CIE luminosity curve
load('1931_CMF_380to780_5.mat')
luminousFunction =  CMF_2(5:65, 3);
% constants
h_planck = 6.6260693*10^(-34); % unit:  J.s
velocity = 2.99792458*10^(8); % unit: m/s

% camera parameters
ADU = 0.60 ; % e-/ADU
gamma = 2.2;
bit_depth = 14;
pixleArea = 23.6*15.6*10^(-6) / (16.9*10^6)*0.50; % assume 55% of effective pixel area

%% photon and electron number calculation
% assumptions
exp_t = 0.001 ;% exposure time: s
QE_curve = rgb;% assume spectral sensitivity approximates QE
illuminance = 500; % lux, baseline illuminance level
brightness_level = [0.01 0.03 0.1 0.3 1 3 10 30 100];
well_capacity  = 12000; % well capacity

irradiance_560 = 1/683 * illuminance/(delta_wavelength*luminousFunction'*light_spd);
irradiance_spectral = irradiance_560* light_spd;
irradiance = sum(irradiance_spectral*delta_wavelength);
illuminance_spectral = 1/683*irradiance_spectral.*luminousFunction;

color_signal_spectral = diag(irradiance_spectral) * reflectance; 

% calculate average photon and electron number: integers
mu_p =round(pixleArea * exp_t / (h_planck*velocity) * (color_signal_spectral .* delta_wavelength).* wavelength*10^-9); 
mu_e = round(1/4 * mu_p' * QE_curve);  % electron number, 1/4 is a result of bayes pattern
    
%% add shot noise to Macbeth chart
% to form a 800*1200 figure presenting 24 patches where each patch has 200*200 pixels
k = length(brightness_level);
psnr_MCCC = zeros(k, 1); % measure of noise level

maxDigitValue = 2^bit_depth-1;
RGB_resultsOrigin = zeros(800, 1200, 3, k);
RGB_noisedResultsOrigin = zeros(800, 1200, 3, k);
RGB_resultsDisplay = zeros(800, 1200, 3, k); % rgb values with gain and gamma correction 
RGB_noisedResultsDisplay =zeros(800, 1200, 3, k);
for j =1:k
    electronMap_noised = zeros(800, 1200, 3);
    for i = 1:24
        patch = brightness_level(j)* reshape(repmat(mu_e(i, :), 200*200, 1), 200, 200, 3);
        % reference image without shot noise
        electronMap(1+200*floor((i-1)/6): 200*(floor((i-1)/6)+1),1+200*(mod(i-1, 6)):200*(mod(i-1, 6)+1), :) = patch;
        % image with shot noise shot noise
        electronMap_noised(1+200*floor((i-1)/6): 200*(floor((i-1)/6)+1),1+200*(mod(i-1, 6)):200*(mod(i-1, 6)+1), :) = round(random('Poisson', patch, size(patch)));
    end
    % convert electron map into a digital RGB image
    electronMap(electronMap>well_capacity) = well_capacity;
    digitMap = round(electronMap/ADU);  % quantization
    digitMap(digitMap>maxDigitValue) = maxDigitValue;
    RGB = digitMap./maxDigitValue;
    RGB_stretch = RGB./max(RGB(:)); % this can be improved by ISO setting
    RGB_gammaC = RGB_stretch.^(1/gamma);
    
    % process the image with shot noise
    electronMap_noised(electronMap_noised > well_capacity) = well_capacity;
    electronMap_noised(electronMap_noised < 0) = 0;
    digitMap_noised = round(electronMap_noised/ADU);  
    digitMap_noised(digitMap_noised > maxDigitValue) = maxDigitValue;
    imRGB_noise = digitMap_noised./maxDigitValue;
    imRGB_noise_stretch = imRGB_noise./max(imRGB_noise(:));
    imRGB_noise_gammaC = imRGB_noise_stretch.^(1/gamma); % gamma correction
    
    % calculate the noise level: we look into the images before any stretch or gamma correction 
    psnr_MCCC(j) = psnr(RGB, imRGB_noise);
    
    % restore the whole results
    RGB_resultsOrigin(:,:,:, j) =  RGB;
    RGB_noisedResultsOrigin(:,:,:, j) =  imRGB_noise;
    RGB_resultsDisplay(:,:,:, j) =  RGB_gammaC;
    RGB_noisedResultsDisplay(:,:,:, j) =  imRGB_noise_gammaC;
    
    figure
%     imshow(RGB_gammaC); 
    imshow(imRGB_noise_gammaC); 
    title(['brightness level of ' num2str(illuminance * brightness_level(j)) 'lux']);
%     print(['shotNoiseResults/MCCC_nikonD5100_' num2str(illuminance * brightness_level(j)) 'lux'], '-depsc');
%     print(['shotNoiseResults/MCCC_nikonD5100_illumA_' num2str(illuminance * brightness_level(j)) 'lux'], '-depsc');    
pause(3);
end
psnr_MCCC

%% adding filter
load('filter_graham.mat');
filter = r(5:65)/max(r(5:65)); % normalization
% calculate average photon number
mu_p_filt = round(pixleArea * exp_t / (h_planck*velocity) * (diag(filter)*color_signal_spectral .* delta_wavelength).* wavelength*10^-9); 
mu_e_filt = round(1/4 * mu_p_filt' * QE_curve);  % electron number, 1/4 is a result of bayes pattern

%% add shot noise to Macbeth chart
% to form a 800*1200 figure presenting 24 patches where each patch has 200*200 pixels
k = length(brightness_level);
psnr_MCCC_filtered = zeros(k, 1); % measure of noise level
maxDigitValue = 2^bit_depth-1;
RGB_resultsOrigin = zeros(800, 1200, 3, k);
RGB_noisedResultsOrigin = zeros(800, 1200, 3, k);
RGB_resultsDisplay = zeros(800, 1200, 3, k);
RGB_noisedResultsDisplay =zeros(800, 1200, 3, k);

for j =1:k 
    electronMap_noised = zeros(800, 1200, 3);
    for i = 1:24
        patch = brightness_level(j)* reshape(repmat(mu_e_filt(i, :), 200*200, 1), 200, 200, 3);
        % reference image without shot noise
        electronMap(1+200*floor((i-1)/6): 200*(floor((i-1)/6)+1),1+200*(mod(i-1, 6)):200*(mod(i-1, 6)+1), :) = patch;
        % image with shot noise shot noise
        electronMap_noised(1+200*floor((i-1)/6): 200*(floor((i-1)/6)+1),1+200*(mod(i-1, 6)):200*(mod(i-1, 6)+1), :) = round(random('Poisson', patch, size(patch)));
    end
    % convert electron map into a digital RGB image
    electronMap(electronMap>well_capacity) = well_capacity;
    digitMap = round(electronMap/ADU);  
    digitMap(digitMap>maxDigitValue) = maxDigitValue;
    RGB = digitMap./maxDigitValue;
    RGB_stretch = RGB./max(RGB(:)); % this can be improved by ISO setting
    RGB_gammaC = RGB_stretch.^(1/gamma);
    
    % process the image with shot noise
    electronMap_noised(electronMap_noised > well_capacity) = well_capacity;
    electronMap_noised(electronMap_noised < 0) = 0;
    digitMap_noised = round(electronMap_noised/ADU);  
    digitMap_noised(digitMap_noised > maxDigitValue) = maxDigitValue;
    imRGB_noise = digitMap_noised./maxDigitValue;
    imRGB_noise_stretch = imRGB_noise./max(imRGB_noise(:));
    imRGB_noise_gammaC = imRGB_noise_stretch.^(1/gamma); % gamma correction
    
    % calculate the noise level: we look into the images before any stretch or gamma correction 
    psnr_MCCC_filtered(j) = psnr(RGB, imRGB_noise);
    
    % restore the whole results
    RGB_resultsOrigin(:,:,:, j) =  RGB;
    RGB_noisedResultsOrigin(:,:,:, j) =  imRGB_noise;
    RGB_resultsDisplay(:,:,:, j) =  RGB_gammaC;
    RGB_noisedResultsDisplay(:,:,:, j) =  imRGB_noise_gammaC;
    
    figure
    imshow([RGB_gammaC imRGB_noise_gammaC]);
    imshow(imRGB_noise_gammaC);
    title(['brightness level of ' num2str(illuminance * brightness_level(j)) 'lux']);
%     print(['shotNoiseResults/noised_MCCC_nikonD5100_' num2str(illuminance * brightness_level(j)) 'lux'], '-depsc');
     print(['shotNoiseResults/noised_MCCC_nikonD5100_illumA_' num2str(illuminance * brightness_level(j)) 'lux'], '-depsc');
    pause(3);
end
psnr_MCCC_filtered