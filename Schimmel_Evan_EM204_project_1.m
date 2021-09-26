clear all
close all
clc

% Define power quantities
Php = 180; % hp
Pinlbs = Php * 550 * 12; % in-lb/s

% Define rotational speed quantities
omegaRPM = 2700; % rpm
omegaRAD = omegaRPM * 1/60 * (2*pi); % rad/s

% Calculate torque in shaft
T = (Pinlbs / omegaRAD) * 1/1000; % kip-in

% Define FOS and calculate max allowable shear stress
FOS = 1.5; % unitless
tau_shear = 26; % ksi
tau_max = (tau_shear / FOS); % ksi

% Define 6061-T6 aluminum density
rho = 0.098; % lb/in^3

% Define loop parameters for testing inner radius values
min_test = 1.0; % in
max_test = 1.4; % in
step = 1e-5; % in

i=0;
k=0;
for d = min_test:step:max_test
    k=k+1;
    
    % Calculate radius value for ease
    r = d/2; % in
    
    % Create diameter vector for plotting
    diam(k) = d; % in
    
    % Loop for 0.090 wall thickness
    for w = 0.090 % in
        i = i+1;
        J = (pi*(((r+w)^4)-(r^4)))/2; % in^4
        tau090 = (T*(r+w))/J; % ksi
        weight = rho * ((pi*r^2)-(pi*(r-w)^2)) * 12; % lb/ft
        
        %  Add values to vector for all shaft options
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau090;
        option(i,4) = weight;
        
        % Create shear stress vector for plotting
        tau_090(k) = tau090;
    end
    
    % Loop for 0.100 wall thickness
    for w = 0.100 % in
        i = i+1;
        J = (pi*(((r+w)^4)-(r^4)))/2; % in^4
        tau100 = (T*(r+w))/J; % ksi
        weight = rho * ((pi*r^2)-(pi*(r-w)^2)) * 12; % lb/ft
        
        % Add values to vector for all shaft options
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau100;
        option(i,4) = weight;
        
        % Create shear stress vector for plotting
        tau_100(k) = tau100;
    end
    
    % Loop for 0.125 wall thickness
    for w = 0.125 % in
        i = i+1;
        J = (pi*(((r+w)^4)-(r^4)))/2; % in^4
        tau125 = (T*(r+w))/J; % ksi
        weight = rho * ((pi*r^2)-(pi*(r-w)^2)) * 12; % lb/ft
        
        % Add values to vector for all shaft options
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau125;
        option(i,4) = weight;
        
        % Create shear stress vector for plotting
        tau_125(k) = tau125;
    end
    
end

% Plot diameter and shear stress data for given wall thickness options
figure
plot(diam,tau_090,diam,tau_100,diam,tau_125)
yline(tau_max,'k','Max Allowable Shear Stress')
xlabel('Radius [in]')
ylabel('Shear Stress [ksi]')
legend('0.090" Wall','0.100" Wall','0.125" Wall')
set(gcf, 'color', 'w')

% Identify all shaft options that have shear stress under allowable value
n=0;
for m = 1:length(option)
    if option(m,3) < tau_max
        n=n+1;
        valid(n,1) = option(m,1);
        valid(n,2) = option(m,2);
        valid(n,3) = option(m,3);
        valid(n,4) = option(m,4);
    end
end

% Identify shaft combo with the lowest weight from valid options
[best_weight,index] = min(valid(:,4));

% Write best measurements to new variables
new_diam = valid(index,1);
new_wall = valid(index,2);
new_tau = valid(index,3);
new_weight = valid(index,4);
new_FOS = tau_shear / new_tau;

% Print ideal shaft information to command window
fprintf('--- Ideal shaft information --- \n\n');
fprintf('Inner Diameter: %7.5f in \n',new_diam);
fprintf('Wall Thickness: %7.5f in \n',new_wall);
fprintf('Shear Stress: %7.5f ksi \n',new_tau);
fprintf('Weight: %7.5f lb/ft \n',new_weight);
fprintf('Factor of Safety: %7.5f \n',new_FOS);