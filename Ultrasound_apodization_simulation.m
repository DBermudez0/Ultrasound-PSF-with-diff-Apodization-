c = 1540;            % Speed of sound [m/s].
f0 = 9e6;            % Signal center frequency [Hz].
lambda = c / f0;     % Wavelength [m]

pitch = lambda;      % Distance between transducer elements. [m]
kerf = pitch * 0.1;  % Empty space between transducer elements [m]
width = pitch * 0.9; % Width of elements [m]
height = width;      % Height of elements [m]

N_elements = 128; % Number of transducer elements.
element_x = ((-65:64)*pitch)+pitch/2; % 1xN_elements array of element lateral position [m]
% The following assertions will give hints if there was a mistake.
assert(element_x(1) < 0, ['The first element lateral position "element_x(1)" should be negative. (it is currently: ' num2str(element_x(1)) ')']);
assert(all(diff(element_x) > 0), 'Element lateral positions should increase monotonically');
assert(element_x(end/2) < 0, ['The ' num2str(N_elements/2) 'th element lateral position should be negative (it is currently: ' num2str(element_x(end/2)) ')']);
assert(element_x(end/2+1) > 0, ['The ' num2str(N_elements/2+1) 'th element lateral position should be positive (it is currently: ' num2str(element_x(end/2+1)) ')']);
assert(all(abs(element_x(1:end/2) + element_x(end:-1:end/2+1)) <= eps('single')), 'The element positions should show symmetry around x = 0.');
assert(all(abs(diff(element_x) - pitch) <= eps('single')), 'The spacing between element positions should be equal to the pitch.');
disp('Element positions OK!');  % No assertion failed :)
% Active emitting element (choose one).
active_element_idx = N_elements/2;

% Scatterer position (choose where you would like it to be)
scatterer_pos_x = 0;
scatterer_pos_z = 20/1000;

active_element_x = element_x(active_element_idx);  % 1x1 array of active element lateral position [m].

% Setup Field-II simulation
if ~exist('field_initialized', 'var')
    field_init(0);
    fs = f0 * 30;  % Sampling rate [Hz].
    set_sampling(fs);
    field_initialized = true;
end

% Aperture geometry.
Tx = xdc_linear_array (N_elements, width, height, kerf, 1, 1, [0, 0, 0]); % Transmit aperture.
Rx = xdc_linear_array (N_elements, width, height, kerf, 1, 1, [0, 0, 0]); % Receive aperture.
xdc_focus_times(Tx, 0, zeros(1, N_elements));  % Clear focusing delays, as we implement them manually.
xdc_focus_times(Rx, 0, zeros(1, N_elements));

% Set active element using Tx apodization mask.
element_mask = zeros(1, N_elements);
element_mask(active_element_idx) = 1;
xdc_apodization(Tx, 0, element_mask);

% A 4-cycle sinusoid emission (you may vary the number of cycles).
cycles = 4;
excitation_osci = sin((0:1/fs:cycles/f0)'*f0*pi*2);  % Pulse oscillation.  
excitation_envelope = hann(numel(excitation_osci)); % Pulse envelope.
excitation = excitation_osci .* excitation_envelope;
excitation = excitation - mean(excitation);
xdc_excitation(Tx, excitation');

% Set scatterer positions.
scatter_pos = [scatterer_pos_x(:), scatterer_pos_x(:)*0, scatterer_pos_z(:)];
scatter_amplitude = ones(size(scatter_pos, 1), 1);

% Use Field-II to simulate the received signals from the emission!
[rf_data, t0] = calc_scat_multi(Tx, Rx, scatter_pos, scatter_amplitude);
rf_data = rf_data ./ max(abs(rf_data(:))); % Rescale simulation output.

% Time axis for excitation.
t_axis0 = (0:1/fs:(numel(excitation)-1)/fs)';

% Move t = 0 to the middle of the excitation.
t_offset = mean(t_axis0);
t0 = t0 - t_offset;
t_axis = t0+(0:size(rf_data, 1)-1)'/fs;
t_axis0 = t_axis0 - t_offset;

figure(2);
subplot(211);
plot(t_axis0*1e6, excitation, 'DisplayName', 'Excitation');
grid on;
xlabel('time [\mus]');
ylabel('Amplitude [V]');
title('Transmitted signal');
legend;
subplot(212);
imagesc(1:N_elements, t_axis * 1e6, rf_data, [-1, 1]);
% imagesc(1:N_elements, t_axis * 1e6, abs(hilbert(rf_data)));
hold on;
tof = plot(nan, nan, '.-r', 'DisplayName', 'Predicted time of flight');
hold off;
colorbar;
xlabel('Element number');
ylabel('time [\mus]');
title('Received signal from each element');
x=0:128;
apodization = 1./(((x-128/2)/26).^2+1);          %1./(((x-128/2)/14).^2+1);             %1 - ((x - N_elements/2).^2) / (N_elements/2)^2;              %hanning(N_elements);           %ones(1, N_elements);

% Define image dimensions.
h = 40/1000;  % [m]
w = (element_x(end)-element_x(1))*1.5;
Nz = ceil(8*h/lambda);  % Sample at a proper rate.
Nx = ceil(4*w/lambda);
x = linspace(-w/2, w/2, Nx);  % Lateral positions.
z = linspace(0, h, Nz)';  % Axial positions.
% Active emitting element (choose one).
active_element_idx = N_elements/2;
% Scatterer position (choose where you would like it to be)
scatterer_pos_x = 0;
scatterer_pos_z = 20/1000;

active_element_x = element_x(active_element_idx);  % 1x1 array of active element lateral position [m].
dist_tx = sqrt((active_element_x - scatterer_pos_x)^2 + scatterer_pos_z^2);
dist_rx = sqrt((element_x - scatterer_pos_x).^2 + scatterer_pos_z^2);

% Run the following assertion first. It gives a hint if there was a mistake.
assert(abs(dist_tx-dist_rx(active_element_idx)) <= eps('single'), ...
'"dist_tx" should be equal to "dist_rx(active_element_idx)"');
rf_data_delayed = zeros(numel(t_axis0), size(rf_data, 2));
% Compute the delays.
distance = dist_tx+dist_rx;  % distance [m]
delays = distance/c;  % time of flight [s]
% Delay and sum
for i = 1:numel(delays)-2
    % Perform the delay.
    rf_data_delayed(:, i) = interp1(t_axis, rf_data(:, i), delays(i)+t_axis0, 'spline', 0);
    % rf_data_delayed(:, i) = rf_data_delayed(:, i) .* dist_rx(i) .* dist_tx;
end
% Perform the sum.
DAS_result = sum(rf_data_delayed, 2);
DAS_result = DAS_result / max(abs(DAS_result)) * max(abs(excitation));  % Scale to match first plot.

rf_data_delayed = rf_data_delayed / max(abs(rf_data_delayed(:)));

figure(3);
subplot(211);

imagesc(1:N_elements, t_axis0*1e6, rf_data_delayed, [-1, 1]);
title('Received signals after applying delays.');
colorbar;
xlabel('Element number');
ylabel('time w.r.t. delay [\mus]');

%%
figure(4);
DAS_result_analytic = hilbert(DAS_result);

plot(t_axis0*1e6, real(DAS_result_analytic), 'Display', 'real(analytic(DAS)) = DAS');
hold on;
plot(t_axis0*1e6, imag(DAS_result_analytic), 'Display', 'imag(analytic(DAS))');
plot(t_axis0*1e6, abs(DAS_result_analytic), 'Display', 'abs(analytic(DAS)) = detected envelope');
plot(t_axis0*1e6, excitation_envelope, '--', 'Display', 'True envelope');
hold off;
ylabel('Amplitude [$V/V$]', 'Interpreter', 'latex');
xlabel('time w.r.t. delay [\mus]');
legend;
grid on;
title('Envelope detection');


subplot(212);
plot(t_axis0*1e6, excitation, '--', 'DisplayName', 'Excitation');
hold on;
plot(t_axis0*1e6, DAS_result, 'DisplayName', 'Delay-and-sum result');
hold off;
title('Beamforming (Sum of signals after applying delays)');
ylabel('Amplitude [$V/V$]', 'Interpreter', 'latex');
xlabel('time w.r.t. delay [\mus]');
grid on;
legend;

% Compute the delays.
distance = dist_tx+dist_rx;  % distance [m]
delays = distance/c;  % time of flight [s]

figure(2);
subplot(212);
tof.XData = 1:N_elements;
tof.YData = delays * 1e6;  % Plot predicted time of flight in microseconds.
legend;
% Prepare animated figure.
figure(6);
sc = imagesc(x*1e3, z*1e3, 0);
hold on;
plot(element_x(active_element_idx)*1e3, element_x(active_element_idx)*0, 'ok', 'DisplayName', 'Transmit element');
plot(scatterer_pos_x(:)*1e3, scatterer_pos_z(:)*1e3, 'xb', 'DisplayName', 'Scatterer');
plot(element_x(:)*1e3, element_x(:)*0, '.r', 'DisplayName', 'Receive elements');
hold off;
legend;
axis image;
colormap gray;
xlim(scatterer_pos_x*1e3 + [-2, 2]);  % Zoom on scatterer (disable if you want).
ylim(scatterer_pos_z*1e3 + [-2, 2]);
xlabel('x [mm]');
ylabel('z [mm]');

% Delay and sum beamformer.
rf_data_analytic = hilbert(rf_data);
beamformed_image = 0;
for i = 1:size(rf_data, 2)
    tx_pos_x = element_x(active_element_idx);
    rx_pos_x = element_x(i);
    
    % Compute the distance from tx_pos to (x, z).
    dist_tx = sqrt((tx_pos_x - x).^2 + z.^2);;
    % Compute the distance from rx_pos to (x, z).
    dist_rx = sqrt((rx_pos_x - x).^2 + z.^2);
    % Hint:
    %  The two arrays should be in meters, and the size should be: [Nz, Nx].
    %  Use the same calculation as in parts 6 and 7 of the exercises, 
    %   but replace "scatterer_pos_x" with "x", and "active_element_x" with "tx_pos_x" and so on.
    %  Computations should be vectorized for efficiency (don't use loops!).
    %   This means you need to use element-wise operations, e.g. .^ instead of ^.
    delays = (dist_tx + dist_rx) * 1/c;
    
    backprojection = interp1(t_axis, rf_data_analytic(:, i), delays, 'spline', 0);
    beamformed_image = beamformed_image + apodization(i) .* backprojection;
    
    % Update plot
    sc.CData = abs(beamformed_image);
    title(['Delay and sum beamforming for element 1 to ' num2str(i)]);
    drawnow;
end
title(['Delay and sum beamforming result']);

% Compute the envelope image (magnitude of the beamformed image)
envelope_image = abs(beamformed_image);

% Log compression
B_mode = 20 * log10(envelope_image);

% Normalize by subtracting the peak value
B_mode = B_mode - max(B_mode(:));

figure(7);
imagesc(x*1e3, z*1e3, B_mode, [-60, 0]);
hold on;
plot(element_x(active_element_idx)*1e3, element_x(active_element_idx)*0, 'ok', 'DisplayName', 'Transmit element');
plot(scatterer_pos_x(:)*1e3, scatterer_pos_z(:)*1e3, 'xb', 'DisplayName', 'Scatterer');
plot(element_x(:)*1e3, element_x(:)*0, '.r', 'DisplayName', 'Receive elements');
hold off;
legend;
axis image;
colormap gray;
xlim(scatterer_pos_x*1e3 + [-2, 2]);
ylim(scatterer_pos_z*1e3 + [-2, 2]);
xlabel('x [mm]');
ylabel('z [mm]');
title('B-mode image [dB]');
colorbar;

%% Point Phantom
scatterer_pos_x = 0.5*linspace(-15, 15, 5)/1000;
scatterer_pos_z = 0.5*linspace(10, 40, 5)'/1000;

% Convert the two lines to their cartesian product grid.
scatterer_pos_x = scatterer_pos_x + 0 * scatterer_pos_z;
scatterer_pos_z = 0 * scatterer_pos_x + scatterer_pos_z;

% Set scatterer positions.
scatter_pos = [scatterer_pos_x(:), scatterer_pos_x(:)*0, scatterer_pos_z(:)];
scatter_amplitude = ones(size(scatter_pos, 1), 1);
[rf_data, t0] = calc_scat_multi(Tx, Rx, scatter_pos, scatter_amplitude);
rf_data = rf_data ./ max(abs(rf_data(:))); % Rescale simulation output.

% Time axis for excitation.
t_axis0 = (0:(numel(excitation)-1))'/fs;
% Move t = 0 to the middle of the excitation.
t0 = t0 - mean(t_axis0);
t_axis = t0+(0:size(rf_data, 1)-1)'/fs;
t_axis0 = t_axis0 - mean(t_axis0);

figure(8);
imagesc(1:N_elements, t_axis * 1e6, abs(hilbert(rf_data)));
colorbar;
xlabel('Element number');
ylabel('time [\mus]');
title('Received signal envelope from each element.');

x= 0:128-1;
NW = 3.5;    % Time-bandwidth product (adjust for trade-off between main lobe width and side lobe suppression)
slepian_win = dpss(N_elements, NW,1)';

combined_apod =  1./(((x-128/2)/30).^2+1).*mean(slepian_win,1);         %1./(((x-128/2)/14).^2+1);                   %hanning(N_elements);          %1 - ((x - N_elements/2).^2) / (N_elements/2)^2;            %hanning(N_elements);;                 %ones(1, N_elements);  % <- Insert your apodization

apodization = combined_apod./max(combined_apod);
% Delay and sum beamformer.
w = diff(scatterer_pos_x([1, end]))*1.2;
h = scatterer_pos_z(end)*1.2;

Nz = ceil(8*h/lambda);
Nx = ceil(4*w/lambda);
x = linspace(-w/2, w/2, Nx);  % Lateral positions.
z = linspace(0, h, Nz)';  % Axial positions.

rf_data_analytic = hilbert(rf_data);

figure(9);
sc = imagesc(x*1e3, z*1e3, 0);
hold on;
plot(element_x(active_element_idx)*1e3, element_x(active_element_idx)*0, 'ok', 'DisplayName', 'Transmit element');
plot(scatterer_pos_x(:)*1e3, scatterer_pos_z(:)*1e3, 'xb', 'DisplayName', 'Scatterer');
plot(element_x(:)*1e3, element_x(:)*0, '.r', 'DisplayName', 'Receive elements');
hold off;
xlabel('x [mm]');
ylabel('z [mm]');
axis image; axis xy;
colormap gray;
legend;
drawnow;

beamformed_image = 0;
for i = 1:size(rf_data, 2)
    tx_pos_x = element_x(active_element_idx);
    rx_pos_x = element_x(i);
    
     dist_tx = sqrt((tx_pos_x - x).^2 + z.^2);
    % Compute the distance from rx_pos to (x, z).
    dist_rx = sqrt((rx_pos_x - x).^2 + z.^2);

    delays = (dist_tx + dist_rx) * 1/c;
    backprojection = interp1(t_axis, rf_data_analytic(:, i), delays, 'spline', 0);
    beamformed_image = beamformed_image + apodization(i) .* backprojection .* dist_rx .* dist_tx;
    if mod(i-1, 1) == 0
        sc.CData = abs(beamformed_image);
        drawnow;
    end
end
colorbar;
title('Press any key to continue to B-mode display.');
pause
% Compute the envelope image (magnitude of the beamformed image)
envelope_image = abs(hilbert(beamformed_image));  % Or rf_data_analytic
B_mode = 20 * log10(envelope_image / max(envelope_image(:)));  % Normalize before log

% Normalize by subtracting the peak value
B_mode = B_mode - max(B_mode(:));
sc.CData = B_mode;
imagesc(x*1e3,z*1e3, sc.CData, [-60, 0])

xlabel('x [mm]');
ylabel('z [mm]');
axis image; axis xy;
colormap gray;
title('B-mode of wire phantom.');

% normalize = true;  % Normalize to show resolution, not ampltiude uniformity.
% depths = unique(scatterer_pos_z);
% laterals = unique(scatterer_pos_x);
% xline_axis = linspace(-1, 1, 401)/1000;
% figure;
% for i = 1:numel(depths)
%     x_line = B_mode(find(z>=depths(i), 1), :);
% 
%     for j = 1:numel(laterals)
%         x_line_sub = interp1(x, x_line, laterals(j) + xline_axis, 'spline', 0);
%         if normalize
%             x_line_sub = x_line_sub - max(x_line_sub(:));
%         end
%         subplot(100+numel(laterals)*10+j);
%         plot(xline_axis*1e3, x_line_sub, 'DisplayName', ['z = ' num2str(depths(i)*1e3) 'mm']);
%         xlabel('x offset [mm]');
%         title(['PSF at x = ' num2str(laterals(j)*1e3) 'mm']);
%         ylim([20*log10(1/2), 0]);
%         hold on;
%     end
% end
% subplot(100+numel(laterals)*10+1);
% ylabel('Amplitude [dB]');
% for j = 1:numel(laterals)
%     subplot(100+numel(laterals)*10+j); hold off;
% end
% legend('Location', 'South');
%% CNR
% A = B_mode( ~any( isnan( B_mode ) | isinf( B_mode ), 2 ),: );
% 
% 
% grayImage=A;
% subplot(2, 2, 1);
% imagesc(grayImage);
% caxis([-60, 0]);
% g = gcf;
% g.WindowState = 'maximized';
% fontSize = 14;
% title('Gray Scale Image', 'FontSize', fontSize);
% uiwait(helpdlg('Click and drag to draw the inner circular region.'));
% roi = drawcircle('Color', 'red');
% innerMask = roi.createMask;
% innerXY = roi.Vertices;
% uiwait(helpdlg('Click and drag to draw the outer circular region.'));
% roi = drawcircle('Color', 'green');
% outerMask = roi.createMask;
% outerXY = roi.Vertices;
% % Exclude the inner mask from the outer mask.
% outerMask = xor(outerMask, innerMask);
% % Get the mean in the circles.
% innerMeanGL = mean(grayImage(innerMask));
% outerMeanGL = mean(grayImage(outerMask));
% 
% innerVarianceGL=std(grayImage(innerMask));
% outerVarianceGL=std(grayImage(outerMask));
% stdard=std(grayImage(innerMask));
% CSR =abs(outerMeanGL-innerMeanGL)/sqrt(1/2*(outerVarianceGL^2));
% dB=20*log(CSR);
% 
% subplot(2, 2, 2);
% imshow(innerMask);
% title('Inner Mask Image', 'FontSize', fontSize);
% subplot(2, 2, 3);
% imshow(outerMask);
% title('Outer Mask Image', 'FontSize', fontSize);
% subplot(2, 2, 4);
% imshow(grayImage);
% hold on;
% % Draw the circles.
% plot(innerXY(:, 1), innerXY(:, 2), 'r-', 'LineWidth', 2);
% plot(outerXY(:, 1), outerXY(:, 2), 'g-', 'LineWidth', 2);
% caption = sprintf('Inner mean GL = %.3f, Outer mean GL = %.3f,Outer var GL = %.3f,inner var GL = %.3f, CSR = %.3f', innerMeanGL, outerMeanGL, innerVarianceGL, outerVarianceGL, CSR);
% title(caption, 'FontSize', fontSize);
% 
% %% sharpeness
% A = B_mode( ~any( isnan( B_mode ) | isinf( B_mode ), 2 ),: );
% % by Tolga Birdal
% % Sharpness estimation using gradients
% % This is the most simple and basic measure of sharpness. Although many
% % great other techniques exist, even this primitive estimation algorithm
% % finds itself decent practical usage.
% % The file reads an image and iteratively smooths it more and more to
% % present how sharpness reduces.
% % The method is tested against mean and unsharp (sharpen) filtering.
% %
% % Usage: Just execute this file. (Remember to have lena.jpg next to this
% % file)
% % All results will be printed to matlab console.
% %
% % Copyright (c) 2011, Tolga Birdal <http://www.tbirdal.me>
% function [sharpness]=estimate_sharpness_test(A)
% 
% 
% G=double(A);
% % measure the sharpness of original image
% sharpness=estimate_sharpness(G);
% disp(['Sharpness of original image: ' num2str(sharpness)]);
% % iteratively smooth and measure sharpness
% for i=3:2:11
%     F=imfilter(G, fspecial('average',i), 'replicate');
%     sharpness=estimate_sharpness(F);
%     disp(['Sharpness after mean filtering: ' num2str(sharpness) '.  Kernel Size: ' num2str(i)]);
% end
% disp(' ');
% % iteratively sharpen the image and measure sharpness
% for i=1:-0.25:0
%     F=imfilter(G, fspecial('unsharp',i), 'replicate');
%     sharpness=estimate_sharpness(F);
%     disp(['Sharpness after unsharp masking: ' num2str(sharpness) '.  Alpha: ' num2str(i)]);
% end
% end
% % Estimate sharpness using the gradient magnitude.
% % sum of all gradient norms / number of pixels give us the sharpness
% % metric.
% function [sharpness]=estimate_sharpness(G)
% [Gx, Gy]=gradient(G);
% S=sqrt(Gx.*Gx+Gy.*Gy);
% sharpness=sum(sum(S))./(numel(Gx));
% end
% 
% [sharpness]=estimate_sharpness_test(A)
% 
