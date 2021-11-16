function [a_ws] = getWignerSeitzRadius(n)
% n (mat double): plasma number density (cm^-3)
% a_ws (mat double): Wigner-Seitz radius (cm) - i.e., average interparticle spacing 

a_ws = (3./(4*pi.*n)).^(1/3);

end