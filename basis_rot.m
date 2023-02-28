function [ D_rot ] = basis_rot( D, angle )

%Basis transformation of the distortion matrix by coordinate axes
%rotation
%   inputs:
%       D -- Distortion matrix 
%       angle -- angle by which you want to rotate the basis vectors of the
%       distortion matrix

%   outputs:
%       D_rot -- Distortion matrix after rotating the basis vectors by the
%       given angle

%%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University. Function added by Hari (Muller group).
%Last updated Sept 17, 2021.

[N_x1, N_x2] = size(D.m11);

D_rot = struct('m11', nan(N_x1,N_x2), 'm22', nan(N_x1,N_x2),'m12', nan(N_x1,N_x2), 'm21', nan(N_x1,N_x2));

mat = zeros(2,2);
rot = zeros(2,2);
rot(1,1) = cosd(angle);
rot(1,2) = -sind(angle);
rot(2,1) = sind(angle);
rot(2,2) = cosd(angle);

for i =1:N_x1
    for j = 1:N_x2
       
       mat(1,1) = D.m11(i,j);
       mat(1,2) = D.m12(i,j);
       mat(2,1) = D.m21(i,j);
       mat(2,2) = D.m22(i,j);
       
       temp_mat = rot * mat / rot; 
       
       D_rot.m11(i,j) = temp_mat(1,1);
       D_rot.m12(i,j) = temp_mat(1,2);
       D_rot.m21(i,j) = temp_mat(2,1);
       D_rot.m22(i,j) = temp_mat(2,2);
    end
end
end