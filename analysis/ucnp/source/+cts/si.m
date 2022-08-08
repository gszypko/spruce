% defines constants in si units
classdef si
    properties (Constant)
        kB = 1.380649e-23; % Boltzmann constant
        mE = 9.1093837e-31; % electron mass
        mI = 1.455e-25; % Sr+ mass
        c = 2.99792458e8; % speed of light
        e = 1.60217663e-19; % electric charge
        eps0 = 8.85418782e-12; % permittivity of free space
        u0 = 4e-7*pi; % permeability of free space
    end
end