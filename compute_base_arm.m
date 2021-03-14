function momentArm = compute_base_arm(pitch, roll, print)

pitch = - pitch;

% Create points
maxX = 3;
maxY = 3;
maxZ = 3;

nPoints = 50;
xVals = linspace(-maxX, maxX, nPoints);
yVals = linspace(0, maxY, nPoints);
zVals = linspace(-maxZ, maxZ, nPoints);
[x,y,z] = meshgrid(xVals, yVals, zVals);
points = [x(:), y(:), z(:)];

zCurveBase = maxY * (points(:,1) / maxX) .^ 2;
xCurve = zCurveBase - points(:,2) <= 0;
zCurve = maxY * (points(:,3) / maxZ) .^ 2 + zCurveBase - points(:,2) <= 0;
inChicken = xCurve & zCurve;

% Density
rhoWater = 16;
rhoBallast = 15;
rhoPartial = 2;
ballastHeight = 2;
rho = rhoPartial * (points(:,2) > ballastHeight) + rhoBallast * (points(:,2) <= ballastHeight);

% Rotate points
rotation = [cosd(roll), -sind(roll) 0; sind(roll), cosd(roll) 0; 0, 0, 1]...
    * [1, 0, 0; 0, cosd(pitch), -sind(pitch); 0, sind(pitch), cosd(pitch)];
points = points * rotation;

% Plot
if print
    figure(1);
    clf;
    surf(reshape(inChicken .* points(:,1), size(x, 1), []), reshape(inChicken .* points(:,3), size(z, 1), []), reshape(inChicken .* points(:,2), size(y, 1), []));
    axis equal
end

% Center of Mass
dX = xVals(2) - xVals(1);
dY = yVals(2) - yVals(1);
dZ = zVals(2) - zVals(1);
dV = dX * dY * dZ;
boatMasses = inChicken .* dV .* rho;
boatMass = sum(boatMasses);
centerOfMass = sum(points .* boatMasses) / boatMass;
if print
    boatMass = boatMass
    centerOfMass = centerOfMass
    volume = sum(inChicken .* dV)
    rhoTotal = boatMass / volume
    displacementRatio = rhoTotal / rhoWater
end

% Center of Buoyancy
function waterMasses = water_masses(d)
    inChickenAndBelowWater = inChicken & (points(:,2) <= d);
    waterMasses = inChickenAndBelowWater * dV * rhoWater;
end
function massDiff = compute_diff(d)
    massDiff = boatMass - sum(water_masses(d));
end
d = fzero(@compute_diff, 0);
waterMasses = water_masses(d);
centerOfBuoyancy = sum(points .* waterMasses) / sum(waterMasses);
if print
    d = d
    centerOfBuoyancy = centerOfBuoyancy
end

momentArm = [centerOfBuoyancy(1) - centerOfMass(1); centerOfBuoyancy(3) - centerOfMass(3)];

end

