% =============================================================
% simulate_PSFs.m
% Generate 2D PSFs for a 10L4 probe using Field II
% =============================================================

clear all; close all;
field_init(0);

% -------------------------------------------------------------
% Transducer parameters (approximate for Siemens 10L4)
% -------------------------------------------------------------
f0 = 5e6;              % center frequency [Hz]
fs = 100e6;            % sampling frequency [Hz]
c  = 1540;             % speed of sound [m/s]
lambda = c/f0;

no_elem = 128;         % subset for speed
width  = 0.3e-3;       % element width [m]
height = 5e-3;         % elevation size [m]
kerf   = 0.02e-3;      % kerf [m]
focus  = [0 0 40e-3];  % geometric focus

Th = xdc_linear_array(no_elem, width, height, kerf, 1, 1, focus);

% Excitation and impulse
excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(Th, excitation);
impulse = gauspuls(-2/f0:1/fs:2/f0, f0, 0.6);
xdc_impulse(Th, impulse);

% -------------------------------------------------------------
% Imaging grid for PSF
% -------------------------------------------------------------
Nx = 65;                         % lateral samples
Nz = 200;                        % depth samples
lat_range = linspace(-5e-3, 5e-3, Nx);    % +/- 5 mm
z_range   = linspace(20e-3, 50e-3, Nz);   % 20–50 mm

% Helper function for PSF simulation
simulate_PSF_func = @(apo) simulate_psf_grid(Th, apo, lat_range, z_range, fs, c);

% -------------------------------------------------------------
% Original PSF (rectangular apodization)
% -------------------------------------------------------------
apo_rect = ones(1,no_elem);
xdc_apodization(Th,0,apo_rect);
PSF_orig = simulate_PSF_func(apo_rect);

% -------------------------------------------------------------
% Apodized PSF (Hanning example)
% -------------------------------------------------------------
apo_hann = gausswin(no_elem)';    %hamming(no_elem)';
x= 0:no_elem-1;
NW = 0.25;    % Time-bandwidth product (adjust for trade-off between main lobe width and side lobe suppression)
slepian_win = dpss(no_elem, NW,1)';

%apo_hann = 1./(((x-no_elem/2)/8).^2+1).*mean(slepian_win,1);         %1./(((x-128/2)/14).^2+1);                   %hanning(N_elements);          %1 - ((x - N_elements/2).^2) / (N_elements/2)^2;            %hanning(N_elements);;                 %ones(1, N_elements);  % <- Insert your apodization

xdc_apodization(Th,0,apo_hann);
PSF_hann = simulate_PSF_func(apo_hann);

% -------------------------------------------------------------
% Save for later
% -------------------------------------------------------------
save('PSFs_10L4.mat','PSF_orig','PSF_hann','lat_range','z_range','Nx','Nz');

% -------------------------------------------------------------
% Plot
% -------------------------------------------------------------
figure;
subplot(1,2,1);
imagesc(lat_range*1e3, z_range*1e3, 20*log10(PSF_orig./max(PSF_orig(:))));
title('Rect Apodization PSF'); colormap gray; colorbar; caxis([-60 0]);
xlabel('Lateral [mm]'); ylabel('Depth [mm]');

subplot(1,2,2);
imagesc(lat_range*1e3, z_range*1e3, 20*log10(PSF_hann./max(PSF_hann(:))));
title('Hanning Apodization PSF'); colormap gray; colorbar; caxis([-60 0]);
xlabel('Lateral [mm]'); ylabel('Depth [mm]');

disp('PSFs saved to PSFs_10L4.mat');
