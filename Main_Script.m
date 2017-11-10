% Inverse Radon Transform Script

%% Part 1: Generation of the phantom.

% Here, the domain is a cube. We are concerne with half the width of the
% domain.
d_half = 200;

% Since it's a discrete grid, there needs to be a grid spacing.
% For simplicity, the voxels are cubes...
[x, y, z] = meshgrid(-d_half:d_half, -d_half:d_half, -d_half:d_half);

% Below is the catalogue of available phantoms to use. Feel free to make 
% your own if you wish:

% Gaussian sphere
% phantom_grid = (1 - heaviside(sqrt(x.^2+y.^2+z.^2) - 70)) .* exp(-(x.^2+y.^2+z.^2) / 1000);

% Monochromatic sphere
phantom_grid = (1 - heaviside(sqrt(x.^2+y.^2+z.^2) - 20));
% phantom_grid = heaviside(sqrt(x.^2+y.^2+z.^2) - 30);

% Rectangle
% phantom_grid = (heaviside(x + 100) - heaviside(x - 100)) .* (heaviside(y + 100) - heaviside(y - 100)) .* (heaviside(z + 100) - heaviside(z - 100));
% phantom_grid = 1 - (heaviside(x + 10) - heaviside(x - 10)) .* (heaviside(y + 10) - heaviside(y - 10)) .* (heaviside(y + 10) - heaviside(y - 10));

% Dual x-pillars
% phantom_grid = (heaviside(y - 10) - heaviside(y + 10));

% NOTE: Due to the fact that the grid is discretised, the grid index
% gives the lower left location of each voxel.

% phantom_grid = repmat(phantom('Modified Shepp-Logan', 2*d_half), [1, 1, 2*d_half]);

figure(1)
imagesc(phantom_grid(:,:,d_half))
colormap gray

% Slice and dice... then plot!
% figure(1)
% slice(x, y, z, phantom_grid, 0, 0, 0)
% colormap gray
% xlabel x
% ylabel y

%% Part 2: Forward Problem

% We allocate an array that stores the number of results for the forward
% problem. 

voxel_iterations = length(phantom_grid) - 1;

% The angular resolution...
d_theta = pi / 80;
angle_iterations = length(0:d_theta:(pi/2-d_theta));
forward_array  = zeros(voxel_iterations, angle_iterations);

for i = 1:2*d_half
    
    voxel_centre = -d_half + ([i, 1, d_half] - 1) + 0.5;
    voxel_value = 0;
    for t = 2:angle_iterations-1
        
        angle = d_theta * (t - 1);
                
        % Compute the rays...
        v_0 = [tan(angle), 1, 0];
        v_0 = v_0 / norm(v_0);
        
        v_1 = [-v_0(1), v_0(2), v_0(3)];
        
        % To compute the end points, let's throw them far away from
        % the domain and let the clipping take care of matters...
        r_0 = voxel_centre + d_half * sqrt(6) * v_0;
        r_1 = voxel_centre + d_half * sqrt(6) * v_1;
        
        % We now perform the ray integrals...
        [ray_01, ray_11] = clip_ray(r_0, voxel_centre, d_half);
        [ray_12, ray_02] = clip_ray(voxel_centre, r_1, d_half);
        
        voxel_value = SLLineIntegralApproximator(phantom_grid, d_half, ray_01, ray_11);
        voxel_value = voxel_value + SLLineIntegralApproximator(phantom_grid, d_half, ray_02, ray_12);
        
        forward_array(i, t) = voxel_value;
        
    end
end

%% Part 3: Inverse Problem

inverse_array = zeros(voxel_iterations, voxel_iterations);

for i = 1:2*d_half
   for j = 1:2*d_half
       
       voxel_value = 0;
       
       for angle_index = 2:angle_iterations-1
           
           theta = d_theta * (angle_index - 1);
           t = tan(theta);
           
           if round(i + j * t) <= 2*d_half
               temp_val = forward_array(round(i + j * t), angle_index) / sqrt(1 + t^2) * tan(d_theta);
               voxel_value = voxel_value + temp_val;
           end
           
           if round(i - j * t) >= 1
               temp_val = forward_array(round(i - j * t), angle_index) / sqrt(1 + t^2) * tan(d_theta);
               voxel_value = voxel_value + temp_val;
           end
           
       end
       
       inverse_array(i, j) = voxel_value;
       
   end
end

figure(2)
imagesc(inverse_array)
colormap gray
