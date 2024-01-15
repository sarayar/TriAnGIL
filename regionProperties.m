function P = regionProperties(I, L)

PADDING = 1.5;

props = {
    'Area'
    'Perimeter'
    'Eccentricity'
    'MajorAxisLength'
    'MinorAxisLength'
    'EquivDiameter'
    'Solidity'
    'Orientation'
    'Centroid'
    'WeightedCentroid'
    };

% NOTE: the image is inverted
P = regionprops(L, 255-I, props);

% NOTE: the image is inverted
% extract separate images and label images (masks) for each region
[R M D] = labelmatrixToROI(255-I, L, PADDING, [], P);

% calculate additional properties
for k = 1:length(P)
    
    c = P(k).Centroid;
    w = P(k).WeightedCentroid;
    d = P(k).EquivDiameter;
    a = P(k).Area;
    l = P(k).Perimeter;
    
    region = R{k};
    regionMask = M{k};
    distance = D{k};
            
    pixelValues = region(regionMask);
        
    % circularity:
    circularity = 4 * pi * a / l^2;
    
    % elliptical deviation:
    % difference between the region and the elliptical approximation of the
    % region
    ellipse = distance <= 1;
    intersection = ellipse & regionMask;
    ellipticalDeviation = 1 - 2*sum(intersection(:)) / (P(k).Area + sum(ellipse(:)));
    
    % mass displacement:
    % distance betweend centroid and weighted centroid normalized with
    % equivalent diameter
    massDisplacement = sqrt(sum((c - w).^2))/d;
    
    % integrated intensity:
    integratedIntensity = sum(pixelValues);
    
    % mean intensity:
    meanIntensity = mean(pixelValues);
    
    % intensity deviation:
    intensityDeviation = std(pixelValues);
    
    % intensity range:
    intensityRange = prctile(pixelValues, 97.5)...
        - prctile(pixelValues, 2.5);
    
    % the inside boundary is defined as the residual between the image
    % region and its erosion with an isotropic strucutring element with
    % radius 1/8 of the equivalent diameter of the region
    se = strel('disk', round(d/8), 0);
    insideBoundary = xor(regionMask, imerode(regionMask, se));
    outsideBoundary = xor(regionMask, imdilate(regionMask, se));
    
    insideBoundaryValues = region(insideBoundary);
    outsideBoundaryValues = region(outsideBoundary);
    
    % inside boundary intensity statistics:
    meanInsideBoundaryIntensity = mean(insideBoundaryValues );
    insideBoundaryIntensityDeviation = std(insideBoundaryValues );
    insideBoundaryIntensityRange = prctile(insideBoundaryValues , 97.5) - prctile(insideBoundaryValues , 2.5);
    normalizedInsideBoundaryIntensity = meanInsideBoundaryIntensity / meanIntensity;
    
    % outside boundary intensuty statistics:
    meanOutsideBoundaryIntensity = mean(outsideBoundaryValues );
    outsideBoundaryIntensityDeviation = std(outsideBoundaryValues );
    outsideBoundaryIntensityRange = prctile(outsideBoundaryValues , 97.5) - prctile(outsideBoundaryValues , 2.5);
    normalizedOutsideBoundaryIntensity = meanOutsideBoundaryIntensity / meanIntensity;
    
    % boundary saliency:
    boundarySaliency = meanInsideBoundaryIntensity - meanOutsideBoundaryIntensity;
    normalizedBoundarySaliency = normalizedInsideBoundaryIntensity - normalizedOutsideBoundaryIntensity;
        
    % add to structure
    P(k).Circularity = circularity;
    P(k).EllipticalDeviation = ellipticalDeviation;
    P(k).MassDisplacement = massDisplacement;
    P(k).IntegratedIntensity =  integratedIntensity;
    P(k).MeanIntensity = meanIntensity;
    P(k).IntensityDeviation = intensityDeviation;
    P(k).IntensityRange = intensityRange;
    P(k).MeanInsideBoundaryIntensity = meanInsideBoundaryIntensity;
    P(k).InsideBoundaryIntensityDeviation = insideBoundaryIntensityDeviation;
    P(k).InsideBoundaryIntensityRange = insideBoundaryIntensityRange;
    P(k).NormalizedInsideBoundaryIntensity = normalizedInsideBoundaryIntensity;
    P(k).MeanOutsideBoundaryIntensity = meanOutsideBoundaryIntensity;
    P(k).OutsideBoundaryIntensityDeviation = outsideBoundaryIntensityDeviation;
    P(k).OutsideBoundaryIntensityRange = outsideBoundaryIntensityRange;
    P(k).NormalizedOutsideBoundaryIntensity = normalizedOutsideBoundaryIntensity;
    P(k).BoundarySaliency = boundarySaliency;
    P(k).NormalizedBoundarySaliency = normalizedBoundarySaliency;
    
end

end

function [R M D] = labelmatrixToROI(I, L, padding, interpType, P)
% labelmatrixToROI: Extracts regions from an image given with a labelmatrix
% as separate images such that the major and minor axes of the regions are
% aligned with the image axes.

if  ~exist('padding', 'var') || isempty(padding)
    padding = 1;
end

if  ~exist('interpType', 'var') || isempty(interpType)
    interpType = 'linear';
end

% the properties matrix can be given as an input if already computed
if  ~exist('P', 'var') || isempty(P)
    props = {
        'Centroid'
        'Orientation'
        'MajorAxisLength'
        'MinorAxisLength'
        };
    
    P = regionprops(L, props);
end

R = cell(length(P), 1);
M = cell(length(P), 1);
D = cell(length(P), 1);

iI = griddedInterpolant(I, interpType);
iL = griddedInterpolant(L, 'nearest');

warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

for k = 1:length(P)
    
    % location, size and orientation of the image
    c = P(k).Centroid;
    a = P(k).MajorAxisLength;
    b = P(k).MinorAxisLength;
    t = P(k).Orientation;
    
    w = round(a * padding);
    h = round(b * padding);
    phi = degToRad(t);
    
    % form a sampling grid
    x = -(h/2-0.5):(h/2-0.5);
    y = -(w/2-0.5):(w/2-0.5);
        
    [X Y] = ndgrid(x, y);
    
    % rotate and center
    A = rotationTranslationMatrix(phi, c);
    T = maketform('affine', A);
    
    [V U] = tformfwd(T, Y, X);
    
    % sample image and label matrix
    %region = interp2(I, U, V, interpType, 0);
    %regionMask = interp2(L, U, V, '*nearest', 0) == k;
    region = iI(U, V);
    regionMask = iL(U, V)==k;
    
    % NOTE: The sampling of the region mask assumes propper region
    % ordering: k-th region has label k.
    
    % normalized distance from center
    distance = sqrt((2*X/a).^2 + (2*Y/b).^2);
    
    region(isnan(region)) = 0;
    regionMask(isnan(regionMask)) = false;
    
    R{k} = region;
    M{k} = regionMask;
    D{k} = distance;
    
end

warning('on', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

end

function radians = degToRad(degrees)
% degToRad: Converts degrees to radians.

radians = degrees * pi / 180;

end
function A = rotationTranslationMatrix(theta, c)
% rotationTranslationMatrix: Returns a rotation ans translation matrix for
% a given angle and offset.

if ~exist('theta', 'var') || isempty(theta)
    theta = 0;
end

if ~exist('c', 'var') || isempty(c)
    c = zeros(1, 2);
end

A = [ cos(theta) -sin(theta) 0  ;
      sin(theta)  cos(theta) 0  ;
      c(1)        c(2)       1 ];

end

function x = roundToOdd(x)
% roundToOdd: Rounds the number to the nearest odd integer.

x = ceil(x);
x = x + mod(x-1,2);

end

%==========================================================================
% Parameters parsing

function p = getParameters(varargin)

ip = inputParser;

% scales for multi-scale processing in micrometers:
ip.addOptional('scales', 6:13, @(x)checkNumericBounds(x, 0));
%1:5
%6:13

% fast radial symmetry transform properties (see Loy & Zelinski paper):
ip.addOptional('alpha', 1, @(x)checkNumericBounds(x, 0));

ip.addOptional('beta', 60, @(x)checkNumericBounds(x, 0));

ip.addOptional('kappa', 10, @(x)checkNumericBounds(x, 0));

% h parameter for extended regional minima (see imextendedmin.m):
ip.addOptional('h', 0.6, @(x)checkNumericBounds(x, 0));

% do not calculate background markers (faster):
ip.addOptional('noBackgroundMarkers', false, @islogical);

% properties ranges for rule based rejection:
ranges = {'Solidity', [0.875 1],...    
          'MassDisplacement', [0 0.08],...
          'BoundarySaliency', [20 255]
    };

ip.addOptional('featureRanges', ranges, @iscell);

% threshold value for contour merging:
ip.addOptional('Th', 0.2,...
    @(x)checkNumericBounds(x, 0));

% parse properties
ip.parse(varargin{:});
p = ip.Results;

end