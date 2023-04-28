% User input parameters
Vinc = input('Enter the amplitude of the incident wave in milli volts: ')*1e-3; % Amplitude of the incident wave in V
Vref = input('Enter the amplitude of the reflected wave in milli volts: ')*1e-3; % Amplitude of the reflected wave in V
a = input('Enter the width of the waveguide in mm: ')*1e-3; % Width of the waveguide in m
b = input('Enter the height of the waveguide in mm: ')*1e-3; % Height of the waveguide in m
epsr = 1; % Relative permittivity of the waveguide material

% Calculation of Characteristic Impedance, Reflection Coefficient, Load Impedance
Z0 = (60 / sqrt(epsr)) * log10(b/a + sqrt((b/a)^2 - 1)); % Characteristic Impedance of the waveguide
Gamma = (Vref / Vinc); % Reflection Coefficient
ZL = Z0 * (1 + Gamma) / (1 - Gamma); % Load Impedance

% Output Results
fprintf('Characteristic Impedance (Z0) = %f ohm\n', Z0);
fprintf('Reflection Coefficient (Gamma) = %f\n', Gamma);
fprintf('Load Impedance (ZL) = %f + j%f ohm\n', real(ZL), imag(ZL));