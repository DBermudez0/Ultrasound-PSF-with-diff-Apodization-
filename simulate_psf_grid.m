function PSF_2D = simulate_psf_grid(Th, apo, lat_range, z_range, fs, c)
    Nx = length(lat_range);
    Nz = length(z_range);
    PSF_2D = zeros(Nz, Nx);

    for ix = 1:Nx
        for iz = 1:Nz
            scatter_pos = [lat_range(ix), 0, z_range(iz)];
            scatter_amp = 1.0;
            [rf_psf, tstart] = calc_scat_multi(Th, Th, scatter_pos, scatter_amp);
            rf_sum = rf_psf * apo(:); % apply apodization explicitly
            env = abs(hilbert(rf_sum));
            z_axis = (0:length(env)-1)/fs*c/2 + tstart*c/2;
            env_interp = interp1(z_axis, env, z_range, 'linear', 0);
            PSF_2D(:,ix) = env_interp(:);
        end
    end
    PSF_2D = PSF_2D ./ max(PSF_2D(:)); % normalize
end
