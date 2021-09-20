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
T = (Pinlbs / omegaRAD) * 1/1000; % ksi

% Define FOS and calculate max allowable shear stress
FOS = 1.5; % unitless
tau_shear = 26; % ksi
tau_max = (tau_shear / FOS); % ksi

% Define 6061-T6 aluminum density
rho = 0.098; % lb/in^3

% Define loop parameters for testing inner diameter values
min_test = 0.6;
max_test = 1.8;
step = 1e-5;

i=0;
k=0;
for d = min_test:step:max_test
    k=k+1;
    C = d/2;
    diam(k) = d;
    
    for w = 0.090
        i = i+1;
        J = (pi*(((C+w)^4)-(C^4)))/2; % in^4
        tau090 = (T*C)/J; % ksi
        weight = rho * ((pi*C^2)-(pi*(C-w)^2)) * 12; % lb/ft
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau090;
        option(i,4) = weight;
        tau_090(k) = tau090;
    end
    
    for w = 0.100
        i = i+1;
        J = (pi*(((C+w)^4)-(C^4)))/2; % in^4
        tau100 = (T*C)/J; % ksi
        weight = rho * ((pi*C^2)-(pi*(C-w)^2)) * 12; % lb/ft
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau100;
        option(i,4) = weight;
        tau_100(k) = tau100;
    end
    
    for w = 0.125
        i = i+1;
        J = (pi*(((C+w)^4)-(C^4)))/2; % in^4
        tau125 = (T*C)/J; % ksi
        weight = rho * ((pi*C^2)-(pi*(C-w)^2)) * 12; % lb/ft
        option(i,1) = d;
        option(i,2) = w;
        option(i,3) = tau125;
        option(i,4) = weight;
        tau_125(k) = tau125;
    end
    
end

% Plot diameter and shear stress data for various wall thicknesses
figure
plot(diam,tau_090,diam,tau_100,diam,tau_125)
yline(tau_max,'k','Max Allowable Shear Stress')
xlabel('Diameter [in]')
ylabel('Shear Stress [ksi]')
legend('0.090 Wall','0.100 Wall','0.125 Wall')
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
fprintf('Factor of Safety: %7.5f ksi \n',new_FOS);
fprintf('Weight: %7.5f lb/ft \n',new_weight);