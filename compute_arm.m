function momentArm = compute_arm(info)

% Constants
rho50 = 10;  % g/in3
rhoBallast = 20; % g/in3
rhoWater = 16;   % g/in3
nPoints = 50;
ballastHeight = 0.5; % in
ballastDepth = -1; % in

pitch = info(1);
roll = info(2);

% Define the shape (inches)
xVals = linspace(-1, 1, nPoints); % Wing to wing
yVals = linspace(0, 2.5, nPoints); % Vertical
zVals = linspace(-1.5, 1.5, nPoints); % Beak to butt
[x,y,z] = meshgrid(xVals, yVals, zVals);
points = [x(:), y(:), z(:)];
% Volumes
inBody = points(:,1) <= 0.75 & points(:,1) >= -0.75 & points(:,2) <= 1.5 & points(:,3) <= 0.5 & points(:,3);
inWings = (points(:,1) >= 0.75 | points(:,1) <= -0.75) & points(:,2) >= 0.5 & points(:,2) <= 1.5 & points(:,3) >= -1.25 & points(:,3) <= 0.25;
inHead = points(:,1) >= -0.5 & points(:,1) <= 0.5 & points(:,2) >= 1 & points(:,3) >= 0.25 & points(:,3) <= 1;
inBeak = points(:,1) >= -0.5 & points(:,1) <= 0.5 & points(:,2) >= 1.5 & points(:,2) <= 2 & points(:,3) >= 1;
inWattle = points(:,1) >= -0.25 & points(:,1) <= 0.25 & points(:,2) >= 1 & points(:,2) <= 1.5 & points(:,3) >= 1 & points(:,3) <= 1.25;
inChicken = inBody | inWings | inHead | inBeak | inWattle;
% Density (with ballast)
densityFunction = @(y, z) (rhoBallast * (y < ballastHeight) + rho50 * (y >= ballastHeight) ...
                          + rhoBallast * (z < ballastDepth) + rho50 * (z >= ballastDepth)) / 2;
densityFunctionVector = densityFunction(points(:,2),points(:,3));

% Rotate points
rotation = [cosd(roll), -sind(roll) 0; sind(roll), cosd(roll) 0; 0, 0, 1]...
    * [1, 0, 0; 0, cosd(pitch), -sind(pitch); 0, sind(pitch), cosd(pitch)];
points = points * rotation;

% Center of Mass
dX = xVals(2) - xVals(1);
dY = yVals(2) - yVals(1);
dZ = zVals(2) - zVals(1);
dV = dX * dY * dZ;
boatMasses = inChicken * dV .* densityFunctionVector;
boatMass = sum(boatMasses);
centerOfMass = sum(points .* boatMasses) / boatMass;

% Center of Buoyancy
function waterMasses = water_masses(d)
    inChickenAndBelowWater = inChicken & (points(:,2) <= d);
    waterMasses = inChickenAndBelowWater * dV * rhoWater;
end
function massDiff = compute_diff(d)
    massDiff = boatMass - sum(water_masses(d));
end
d = fzero(@compute_diff, 0);
compute_diff(d);
waterMasses = water_masses(d);
centerOfBuoyancy = sum(points .* waterMasses) / sum(waterMasses);

momentArm = [centerOfMass(1) - centerOfBuoyancy(1); centerOfMass(3) - centerOfBuoyancy(3)];

end