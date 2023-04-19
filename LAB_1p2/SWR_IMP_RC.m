% Prompt user to input values
a = input('Enter waveguide width (mm): ')/1000;
b = input('Enter waveguide height (mm): ')/1000;
f = input('Enter frequency (GHz): ')*1e9;
Vmax = input('Enter maximum voltage (mV): ')/1000;
Vmin = input('Enter minimum voltage (mV): ')/1000;
d = input('Enter distance from load to measurement plane (mm): ')/1000;

% Constants
eps_r = 1; % relative permittivity of air
c = 299792458; % speed of light in vacuum
mu = 4*pi*10^-7; % permeability of vacuum

% Calculate cutoff frequency
fc = c/(2*pi)*1/sqrt((a/b)^2 - 1);
disp(['Cutoff frequency (fc) = ' num2str(fc/1e9) ' GHz']);

% Calculate propagation constant and wave impedance
if a > b
    k = 2*pi*f/c*sqrt(eps_r)*(sqrt(1 - (b/a)^2));
    Z0 = 377/sqrt(eps_r)*(b/a)/sqrt(1 - (b/a)^2);
else
    k = 2*pi*f/c*sqrt(eps_r)*(sqrt((a/b)^2 - 1));
    Z0 = 377*sqrt(eps_r)/(sqrt((a/b)^2 - 1));
end

% Calculate reflection coefficient
if f < fc
    Gamma = 0;
    SWR = 1;
    ZL = Inf;
else
    beta_L = sqrt((a/b)^2 - 1)*k;
    if d > 0
        ZL = -1j*Z0*tan(beta_L*d);
    else
        ZL = 1j*Z0*cot(beta_L*d);
    end
    Gamma = (ZL - Z0)/(ZL + Z0);
    SWR = (1 + abs(Gamma))/(1 - abs(Gamma));
end

% Display results
disp(['Waveguide width (a) = ' num2str(a*1000) ' mm']);
disp(['Waveguide height (b) = ' num2str(b*1000) ' mm']);
disp(['Frequency (f) = ' num2str(f/1e9) ' GHz']);
disp(['Maximum voltage (Vmax) = ' num2str(Vmax*1000) ' mV']);
disp(['Minimum voltage (Vmin) = ' num2str(Vmin*1000) ' mV']);
disp(['Distance from load to measurement plane (d) = ' num2str(d*1000) ' mm']);
disp(['Wave impedance (Z0) = ' num2str(Z0) ' ohms']);
disp(['Reflection coefficient (Gamma) = ' num2str(Gamma)]);
disp(['SWR = ' num2str(SWR)]);
disp(['Load impedance (ZL) = ' num2str(ZL) ' ohms']);
