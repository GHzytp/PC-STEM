function [ D_final ] = basis_non_orth( D_rot, angle )

%Basis transformation of the distortion matrix to transform from a pair of
%non-orthogonal basis vectors to orthogonal basis vectors. Changes basis vector #1 to be
% perpendicular to #2.

%   inputs:
%       D_rot -- Distortion matrix that has been rotated to orient basis vector #2 along the required direction 
%       angle -- angle between basis vectors #1 and #2

%   outputs:
%       D_final -- Final distortion matrix with basis vectors orthogonal to
%       each other

%%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University. Function added by Hari (Muller group).
%Last updated Sept 17, 2021.

[N_x1, N_x2] = size(D_rot.m11);

D_final = struct('m11', nan(N_x1,N_x2), 'm22', nan(N_x1,N_x2),'m12', nan(N_x1,N_x2), 'm21', nan(N_x1,N_x2));

transform = zeros(2,2);
transform(1,1) = sind(angle);
transform(1,2) = 0;
transform(2,1) = cosd(angle);
transform(2,2) = 1;

for i =1:N_x1
    for j = 1:N_x2
       
       mat(1,1) = D_rot.m11(i,j);
       mat(1,2) = D_rot.m12(i,j);
       mat(2,1) = D_rot.m21(i,j);
       mat(2,2) = D_rot.m22(i,j);
       
       temp_mat = transform * mat / transform;
       
       D_final.m11(i,j) = temp_mat(1,1);
       D_final.m12(i,j) = temp_mat(1,2);
       D_final.m21(i,j) = temp_mat(2,1);
       D_final.m22(i,j) = temp_mat(2,2);
    end
end
end