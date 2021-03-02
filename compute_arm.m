function momentArm = compute_arm(info)

% Constants
rho50 = 10;  % g/in3
rhoWater = 16;   % g/in3
nPoints = 50;

pitch = info(1);
roll = info(2);

% Define the shape (inches)
xVals = linspace(-1, 1, nPoints); % Wing to wing
yVals = linspace(0, 2.5, nPoints); % Vertical
zVals = linspace(-1.5, 1.5, nPoints); % Beak to butt
[x,y,z] = meshgrid(xVals, yVals, zVals);
points = [x(:), y(:), z(:)];

inBody = points(:,1) <= 0.75 & points(:,1) >= -0.75 & points(:,2) <= 1.5 & points(:,3) <= 0.5 & points(:,3);
inWings = (points(:,1) >= 0.75 | points(:,1) <= -0.75) & points(:,2) >= 0.5 & points(:,2) <= 1.5 & points(:,3) >= -1.25 & points(:,3) <= 0.25;
inHead = points(:,1) >= -0.5 & points(:,1) <= 0.5 & points(:,2) >= 1 & points(:,3) >= 0.25 & points(:,3) <= 1;
inBeak = points(:,1) >= -0.5 & points(:,1) <= 0.5 & points(:,2) >= 1.5 & points(:,2) <= 2 & points(:,3) >= 1;
inWattle = points(:,1) >= -0.25 & points(:,1) <= 0.25 & points(:,2) >= 1 & points(:,2) <= 1.5 & points(:,3) >= 1 & points(:,3) <= 1.25;
inChicken = inBody | inWings | inHead | inBeak | inWattle;

% Rotate points
rotation = [cosd(roll), -sind(roll) 0; sind(roll), cosd(roll) 0; 0, 0, 1]...
    * [1, 0, 0; 0, cosd(pitch), -sind(pitch); 0, sind(pitch), cosd(pitch)];
points = points * rotation;

% Center of Mass
dX = xVals(2) - xVals(1);
dY = yVals(2) - yVals(1);
dZ = zVals(2) - zVals(1);
dV = dX * dY * dZ;
boatMasses = inChicken * dV * rho50;
boatMass = sum(boatMasses);
centerOfMass = sum(points .* boatMasses) / boatMass;

% Center of Buoyancy
waterMasses = 0;
waterMass = 0;
function massDiff = compute_diff(d)
    inChickenAndBelowWater = inChicken & (points(:,2) <= d);
    waterMasses = inChickenAndBelowWater * dV * rhoWater;
    waterMass = sum(waterMasses);
    massDiff = boatMass - waterMass;
end
d = fzero(@compute_diff, 0);
compute_diff(d);
centerOfBuoyancy = sum(points .* waterMasses) / waterMass;

momentArm = [centerOfMass(1) - centerOfBuoyancy(1); centerOfMass(3) - centerOfBuoyancy(3)];

end


