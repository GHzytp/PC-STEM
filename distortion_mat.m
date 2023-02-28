function [ D ] = distortion_mat( spotMaps, spotReferences)
%From the EWPC peak positions at each pixel stored in spotMaps, calculates the
%affine transformation required to map to the reference positions.
%   inputs:
%       spotMaps -- struct array containing maps of maximum spot index Q1,Q2
%                 for each spot in spotList. Fields are:
%                 'id' -- a name identifying the spot.
%                 'Q1map' -- map of Q1 index value at peak maximum.
%                 'Q2map' -- map of Q2 index value at peak maximum.
%                 'x1map' -- map of x1 vector component for the spot.
%                 'x2map' -- map of x2 vector component for the spot.
%                 'spotlength' -- map of the spot vector length.
%                 'spotangle' -- map of the spot vector angle in degrees.
%       spotReferences -- struct array containing reference points
%                         corresponding to points in Qmaps. Fields are:
%                         'id' -- a name identifying the spot.
%                         'point' -- reference spot location [q1,q2]
%   outputs:
%       D -- Distortion matrix that represents the affine
%       transformation between EWPC peak positions at each pixel to the
%       reference positions.
%
%This function is part of the PC-STEM Package by Elliot Padgett in the 
%Muller Group at Cornell University.  Last updated June 26, 2019.

%initialize data info
[Nx1,Nx2] = size(spotMaps(1).Q1map);

%prealocate results structs
D = struct('m11', nan(Nx1,Nx2), 'm22', nan(Nx1,Nx2),...
                          'm12', nan(Nx1,Nx2), 'm21', nan(Nx1,Nx2));

%prepare reference point list
referencePoints = [];
for s = 1:length(spotReferences)
    referencePoints = [referencePoints; spotReferences(s).point];
end

%calculate strain across map
for j=1:Nx1
    for k=1:Nx2
        
        dataPoints = [];
        for s = 1:length(spotReferences)
            %center spot
            q1c = spotMaps(s).VectorX1(j,k);
            q2c = spotMaps(s).VectorX2(j,k);
            
            %include in list for tranformation calculation
            dataPoints = [dataPoints; [q1c,q2c]];
        end
        if sum(isnan(dataPoints(:)))
            %nan values mean values are unknown and this point should
            %be skipped

            D.m11(j,k) = nan;
            D.m12(j,k) = nan;
            D.m21(j,k) = nan;
            D.m22(j,k) = nan;

        else
            %Calculate distorion/deformation matrix
            
            M = dataPoints /(referencePoints);
            D.m11(j,k) = M(1,1);
            D.m12(j,k) = M(1,2);
            D.m21(j,k) = M(2,1);
            D.m22(j,k) = M(2,2);
            
        end
    end
end