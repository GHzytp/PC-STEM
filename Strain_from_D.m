function [ StrainComponents,StrainTensors] = Strain_from_D( D, latticeCoords )
% Makes maps of strain tensor components from the distortion matrix
%   inputs:
%       D -- distortion matrix after rotation and transformation to correct
%       basis vectors
%       latticeCoords -- choose coordinate system to reference
%                        strain to.  Value 0 selects "image 
%                        coordinates", with the 1,2 directions 
%                        corresponding to image x1,x2.  Value 1 selects
%                        "lattice coordinates", where the 1,2 directions
%                        refer to the reference directions.
%   outputs:
%       StrainComponents -- struct containing 2D arrays (images) of the
%               strain tensor components, rotation, and strain ellipse 
%               parameters, including:
%                   Eps11 -- The first diagonal element of the strain tensor
%                   Eps22 -- The second diagonal element of the strain tensor
%                   Eps12 -- The 2D shear component of the strain tensor (1,2) off
%                            diagonal
%                   Theta -- The rotation relative to the reference, decomposed from
%                            the strain tensor, in degrees
%                   minAx -- semi-minor axis length of the strain ellipse
%                   majAx -- semi-major axis length of the strain ellipse
%                   strainAngle -- angle of the strain ellipse semi-major
%                                  axis
%       StrainTensors -- struct array containing the strain tensor 'E' and
%                        rotation matrix 'R' describing the distortion at each
%                        real space position, along with the eigenvalues
%                        'evals' and eigenvectors 'evecs' of the strain tensor
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated Sept 17, 2021.

%initialize data info
[Nx1,Nx2] = size(D.m11);

%prealocate results structs
StrainComponents = struct('Eps11', nan(Nx1,Nx2), 'Eps22', nan(Nx1,Nx2),...
                          'Eps12', nan(Nx1,Nx2), 'Theta', nan(Nx1,Nx2),...
                          'majAx', nan(Nx1,Nx2), 'minAx', nan(Nx1,Nx2),...
                          'strainAngle', nan(Nx1,Nx2));
StrainTensors = struct('E',cell(Nx1,Nx2),'R',cell(Nx1,Nx2),...
    'evals',cell(Nx1,Nx2),'evecs',cell(Nx1,Nx2));

%calculate strain across map
for j=1:Nx1
    for k=1:Nx2
        
        X = zeros(2,2);
        X(1,1) = D.m11(j,k);
        X(1,2) = D.m12(j,k);
        X(2,1) = D.m21(j,k);
        X(2,2) = D.m22(j,k);
                   
        [R,U,V] = poldecomp(X); % D = R*U = V*R.  R is a rotation.
        %V is the strain in "world coordinates" and is independent of
        %the reference orientation. U is the strain in "lattice
        %coordinates", corresponding to the reference directions

        if latticeCoords
            E = U-eye(2); % strain tensor
        else
            E = V-eye(2); % strain tensor
        end

        StrainTensors(j,k).R = R;
        StrainTensors(j,k).E = E;

        %get strain tensor elements
        StrainComponents.Eps11(j,k) = E(1,1); %fractional stretch
        StrainComponents.Eps22(j,k) = E(2,2); %fractional stretch
        StrainComponents.Eps12(j,k) = E(1,2); %fractional shear
        StrainComponents.Theta(j,k) = atan2d(R(2,1),R(1,1)); %rotation in degrees

        %get strain elipse parameters
        [v,d]=eig(E);
        [evals,sorting] = sort([d(1,1),d(2,2)]); %sorts in ascending order
        evecs = v(:,sorting);
        StrainComponents.minAx(j,k) = evals(1); 
        StrainComponents.majAx(j,k) = evals(2); 
        StrainComponents.strainAngle(j,k) = atand(evecs(2,2)/evecs(1,2));

        StrainTensors(j,k).evecs = evecs(:,2);
        StrainTensors(j,k).evals = evals';
            
    end
end
end



function [R, U, V] = poldecomp(F)
%POLDECOMP  Performs the polar decomposition of a regular square matrix.
%   [R U V] = POLDECOMP(F) factorizes a non-singular square matrix F such
%   that F=R*U and F=V*R, where
%   U and V are symmetric, positive definite matrices and
%   R is a rotational matrix
%
%   See also EIG, DIAG, REPMAT


% This kind of decomposition is often used in continuum mechanics so it is
% convenient to comment the code that way. From now, we use the matrix
% formalism of tensors. C is the right Cauchy-Green deformation tensor,
% F is the deformation tensor, lambda is the stretch.
%
%Copyright (c) 2014, Zoltán Csáti

% Check input
[m n] = size(F);
if m ~= n
    error('Matrix must be square.');
end

C = F'*F;
[Q0, lambdasquare] = eig(C);
lambda = sqrt(diag((lambdasquare))); % extract the components
% Uinv is the inverse of U and is constructed with the help of Q0. Uinv is
% produced in the same base as F not in the base of its eigenvectors.
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
% Using the definition, R, U and V can now be calculated
R = F*Uinv;
U = R'*F;
V = F*R';
end