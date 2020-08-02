function [ Ls_dc, Rs_dc ] = get_CTMC_decorr_boundaries(sigmas, C, T, num_devs, sigma_dc, num_devs_Y)

Ls_dc_0 = zeros(1,length(sigmas));
Rs_dc_0 = zeros(1,length(sigmas));

if nargin < 6
    num_devs_Y = num_devs;
end

Ls_dc = zeros(size(sigmas));
Rs_dc = zeros(size(sigmas));
for i = 1:length(sigmas)
    Ls_dc(i) = min([Ls_dc_0(i), Rs_dc_0(i), -num_devs_Y*sigma_dc(i)*sqrt(T)]);
    Rs_dc(i) = max([Ls_dc_0(i), Rs_dc_0(i), num_devs_Y*sigma_dc(i)*sqrt(T)]);
end


end



