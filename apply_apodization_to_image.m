function I_new = apply_apodization_to_image(I_in, PSF_orig, PSF_apo,z_range,lat_range)
% Apply approximate apodization effect using precomputed PSFs
% 
% I_in   : grayscale ultrasound phantom image (2D matrix)
% PSF_orig : 2D PSF from Field II with rectangular apodization
% PSF_apo  : 2D PSF from Field II with desired apodization
%
% Returns: reconstructed/apodized phantom image

    % Normalize input image
    I = double(I_in);
    I = I ./ max(I(:));

    [Nz,Nx] = size(I);

    % ------------------------------------------------------
    % Fourier-domain filtering
    % ------------------------------------------------------
    
    % Forward transforms
F_I        = ifftshift(fft2((I)));
F_PSF_orig = ifftshift(fft2((PSF_orig), size(I,1), size(I,2)));
F_PSF_apo  = ifftshift(fft2((PSF_apo),  size(I,1), size(I,2)));

    eps_val = 1e-6;
    F_mod = F_I.* F_PSF_apo;
    
    % Reconstruct apodized image
    I_new = circshift((abs((ifft2((F_mod))))),-190, 1);
    I_new = I_new ./ max(I_new(:));

    % ------------------------------------------------------
    % Display comparison
    % ------------------------------------------------------
    figure;
    subplot(1,2,1);
    imagesc(lat_range*1e3, z_range*1e3, I); axis image; colormap gray;colorbar
    colorbar;
    xlabel('Lateral [mm]'); ylabel('Depth [mm]');
    title('Original Phantom Image');

    [Gx, Gy] = gradient(I); % Compute gradients in x and y directions
    gradient_magnitude_orig = sqrt(Gx.^2 + Gy.^2); % Compute gradient magnitude
    
    % Compute sharpness as the mean gradient magnitude
    sharpness_values = mean(gradient_magnitude_orig(:))
    subplot(1,2,2);
    imagesc(lat_range*1e3, z_range*1e3,I_new); axis image; colormap gray;colorbar
    colorbar; 
    xlabel('Lateral [mm]'); ylabel('Depth [mm]');
    title('Reconstructed with Gaussian');
    [Gx, Gy] = gradient(I_new); % Compute gradients in x and y directions
    gradient_magnitude_apod = sqrt(Gx.^2 + Gy.^2); % Compute gradient magnitude
    
    % Compute sharpness as the mean gradient magnitude
    sharpness_values_apod = mean(gradient_magnitude_apod(:))
end

