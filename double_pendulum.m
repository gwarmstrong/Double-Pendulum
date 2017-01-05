function dsdt = double_pendulum( t, s, m, L)
% REDO DOCUMENTATION !!!!!!!!!!!!!!!!!!!!!!!!
% pendulum(): calculates the derivative for a single time and position of 
%              a damped driven pendulum
% inputs:   t:  the time at which the derivative is being calculated
%           s:  [omega theta] where omega is the current angular velocity
%               and omega is the current angle
%           beta: damping coefficient
%           omega_0: initial angular velocity
%           F:  @(t) a function handle for the driving force at t
% outputs:  [omega_dot; theta_dot] = a vector with the numerical values of
%           the temporal derivative of the position and angular velocity
%   

% z1 = theta1, z2 = theta2, z3 = theta_dot1, z4 = theta_dot2
m1 = m(1);
m2 = m(2);
L1 = L(1);
L2 = L(2);
g = 1;%9.81;
theta1 = s(1);
theta2 = s(2);
theta_dot1 = s(3);
theta_dot2 = s(4);

z1 = theta1;
z2 = theta2;
z3 = theta_dot1;
z4 = theta_dot2;

dsdt = zeros(4,1);
dsdt(1) = theta_dot1;
dsdt(2) = theta_dot2;
dsdt(3) = (-m2*L1*z4^2 * sin(z1-z2) * cos(z1-z2) + ...
    g * m2 * sin(z2)*cos(z1-z2) - ...
    m2 * L2 * z4^2 * sin(z1 - z2) - ...
    (m1+m2)*g*sin(z1))/(L1*(m1+m2) - ...
    m2*L1*cos(z1-z2)^2);

dsdt(4) = (m2*L1*z4^2 * sin(z1-z2) * cos(z1-z2) +...
    g * sin(z1)*cos(z1-z2)*(m1+m2)+ ...
    L1*z4^2 * sin(z1-z2)*(m1 + m2)- ...
    g * sin(z2) * (m1 + m2)) / (L2*(m1+m2) - ...
    m2*L2*cos(z1-z2)^2);




end


