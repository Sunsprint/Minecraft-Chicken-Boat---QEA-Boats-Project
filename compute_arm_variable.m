function momentArm = compute_arm_variable(pitch, roll, rhoHollow, rhoBack, ballastHeight, ballastDepth, print)

pitch = -pitch; % Positive is forward

% Create points
nPoints = 50;
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
rhoBottom = rhoHollow; % g/in3
rhoWater = 16;

inBack = points(:,3) <= ballastDepth;
inBottom = ~inBack & (points(:,2) <= ballastHeight);
hollowPoints = ~inBack & ~inBottom;
rho = rhoBack * inBack + rhoBottom * inBottom + rhoHollow * hollowPoints;

% Rotate points
rotation = [cosd(roll), -sind(roll) 0; sind(roll), cosd(roll) 0; 0, 0, 1]...
    * [1, 0, 0; 0, cosd(pitch), -sind(pitch); 0, sind(pitch), cosd(pitch)];
points = points * rotation;

% Plot
if print
    figure(1);
    clf;
    scatter3(points(inChicken,1), points(inChicken,3), points(inChicken,2), 0.1, rho(inChicken));
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

