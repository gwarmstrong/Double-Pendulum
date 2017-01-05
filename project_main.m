clear all;
options = odeset('RelTol', 1e-8,'AbsTol',1e-12);

%%
% Author(s):    George Armstrong <garmstrong@colgate.edu>
%               Erik Pohl <epohl@colgate.edu>
% Institution:  Colgate University 
% Description:  Scripts we used for modeling and analyzing a double pendulum
%               system for our final project for Computational
%               Mechanics (PHYS 451). We modeled both a simple double
%               pendulum system and a damped double pendulum. Our analysis
%               includes phase-space plots, Poincare maps, and first-flip 
%               diagrams for the simple double pendulum.
% Original:     December 2016
% Modified:     January 2017
% References:
%{
Christian, W., & Hampton, D. (2009). Double Pendulum With Poincare Map 
    Model (Version 1.0) [Computer software]. Retrieved January 5, 2017, 
    from http://www.compadre.org/Repository/document/ServeFile.cfm?ID=9384&DocID=1289
Szuminski, Wojciech. "Dynamics of Multiple Pendula." (n.d.): n. pag. 
    University of Wisconsin,	Milwaukee. Web.
%}


%% simple undamped double pendulum
theta0 = [pi/2; pi/2+0.01; 1; 1]; % change IC's here
m = [1 1]; % [mass1 mass2]
L = [1 1]; % [Length1 Length2]
t_span = 0:0.001:700;
[t_sol, s_sol] = ode45(@(t,s) double_pendulum( t, s, m, L), t_span, theta0, options);
theta1 = s_sol(:,1);
theta2 = s_sol(:,2);
theta1_dot=s_sol(:,3);
theta2_dot=s_sol(:,4);

x1 = L(1)*sin(theta1);
x2 = L(1)*sin(theta1) + L(2)*sin(theta2);
y1 = -L(1)*cos(theta1);
y2 = -L(1)*cos(theta1) - L(2)*cos(theta2);

x1_pos = zeros(1,length(x1));
x1_pos(1) = x1(1);
for i = 2:length(x1_pos)
    x1_pos(i) = abs(x1(i)-x1(i-1)) + x1_pos(i-1);
end

mov = VideoWriter('double_pendlulum_m_1_1_l_1_1.avi', 'Uncompressed AVI');%'Motion JPEG 2000');
open(mov);
% pendulum movie
for i = 1:200:length(t_span)/2 % did not use all 700 s b/c the movie was longer than we needed
    clf
    hold on
    plot([x1(i), x2(i)],[y1(i), y2(i)], 'b')
    plot([x1(i), 0],[y1(i), 0], 'b')
    plot(x1(max(1,i-36000):i),y1(max(1,i-36000):i),'r')
    plot(x2(max(1,i-36000):i),y2(max(1,i-36000):i),'b')
    axis([-3,3,-3,3])
    axis square
    xlabel('x position (m)')
    ylabel('y position (m)')
    title('Undamped Simple Double Pendulum')
    hold off
    M = getframe(gcf);
    writeVideo(mov, M);

end
close(mov);

% phase space plot
plot (mod(theta1+pi,2*pi)-pi,theta1_dot,'.')

%% phase space plots at low energies
theta0 = [theta1s(1),theta2s(1),0,0];
[t_sol, s_sol] = ode45(@(t,s) double_pendulum( t, s, m, L), [0 t_end], theta0, options);
theta1 = s_sol(:,1);
theta2 = s_sol(:,2);
theta1_dot=s_sol(:,3);
theta2_dot=s_sol(:,4);
plot (mod(theta1+pi,2*pi)-pi,theta1_dot,'.')
ylabel('theta1 dot')
xlabel('theta1')
title('Normalized Phase Space Plot of Undamped Pendulum with Low Energy')
plot (mod(theta2+pi,2*pi)-pi,theta2_dot,'.')
ylabel('theta2 dot')
xlabel('theta2')
title('Normalized Phase Space Plot of Undamped Pendulum with Low Energy')


%% poincare maps
% constants
g = 1; % have to change g in double_pendulum.m to reflect this
m = [3 1]; L = [2 1]; m1 = m(1); m2 = m(2); l1 = L(1); l2 = L(2);

% when looking at different energy levels, we have to find acceptable
% combinations of initial thetas. For a single poincare map, only a single
% of the following blocks of codes for energy levels should be used

% first energy level
E = -8.95;
theta1s = -0.1:0.008:0.1;

% second energy level
E = -7.3;
numsteps = 20;
min = 0.5;
max = 0.66;
step_size = (max-min)/numsteps;
first_half = min:step_size:max;
second_half = [];
theta1s = [first_half, second_half];

% third energy level
E = -6.9;
numsteps = 10;
min = 0.5;
max = 0.74;
step_size = (max-min)/numsteps;
theta1s = min:step_size:max;

% for energy levels with no initial velocity
% e.g. -8.95, 0, 1
E = -5.3;
numsteps = 20;
minimum = min([acos((-1-E)/8),real(acos((1-E)/8))]) + 0.01;%1.45;
maximum = max([acos((-1-E)/8),real(acos((1-E)/8))]) - 0.01;%1.68
step_size = (maximum-minimum)/numsteps;
first_half = minimum:step_size:maximum;
second_half = -maximum:step_size:-minimum;
theta1s = [first_half, second_half];

theta2_fh =@(theta1) acos((E+m(1)*L(1)*cos(theta1)+m(2)*L(1)*cos(theta1))/(-L(2)*m(2)));
theta2s = theta2_fh(theta1s);

% for higher energy levels with theta1=theta2=pi/2
E = 150;
thetadot1s = -1.2:0.2:1.2; %start with low initial velocities
a = 2*E - m1*l1^2.*thetadot1s + m2*l1^2.*thetadot1s.^2;
b = (2*m2*l2*l1.*thetadot1s)./(2*E-m1*l1^2.*thetadot1s+m2*l1^2.*thetadot1s.^2);
first_half = 2*sqrt(a)./(sqrt(a.*b.^2 + 4*m2*l2) + sqrt(a).*b);
second_half = 2*sqrt(a)./(-sqrt(a.*b.^2 + 4*m2*l2) + sqrt(a).*b);
thetadot1s = [thetadot1s, thetadot1s];
thetadot2s = [first_half, second_half];


% solving Poincare maps
t_span = 0:0.001:700;
t_end = 700;
p1 = @(theta1,theta2,thetadot1,thetadot2) ...
    (m(1)+m(2))*L(1)^2*thetadot1 + m(2) * L(1) * L(2) * thetadot2 * cos(theta1-theta2);
p2 = @(theta1,theta2,thetadot1,thetadot2) ...
    m(2)*L(2)^2*thetadot2 + m(2)*L(1)*L(2)*thetadot1*cos(theta1-theta2);

hold on
for k = 1:length(theta1s)
    theta0 = [theta1s(k),theta2s(k),0,0]; % for energy levels with no initial velocity
    %theta0 = [pi/2, pi/2, thetadot1s(k),thetadot2s(k)]; % use this for 
    [t_sol, s_sol] = ode45(@(t,s) double_pendulum( t, s, m, L), [0 t_end], theta0, options);
    preserved_s = zeros(1,2);
    epsilon = 1e-3;
    for i = 1:length(t_sol)
        % Plotting conditions from Christian & Hampton (2009)
%         if abs(s_sol(i,2)) < epsilon && s_sol(i,3) * L(1)/L(2) * cos(s_sol(i,1)) > 0
%             preserved_s(i,1) = s_sol(i,1);
%             preserved_s(i,2) = s_sol(i,3);
%         end
        % Plotting conditions from from Szuminski n.d.
        momentum1 = p1(s_sol(i,1),s_sol(i,2),s_sol(i,3),s_sol(i,4));
        momentum2 = p2(s_sol(i,1),s_sol(i,2),s_sol(i,3),s_sol(i,4));
        if abs(s_sol(i,1)) < epsilon && momentum1 > 0
            preserved_s(i,1) = s_sol(i,2);
            preserved_s(i,2) = s_sol(i,4);%momentum2;
        end
    end
    plot(mod(preserved_s(:,1)+pi,2*pi)-pi,preserved_s(:,2),'.')
    %axis([-pi pi -5 10]) %for larger energy values it will help to adjust
    % the axes
    title('Poincare map for E = -5.3 J')
    xlabel('\theta_2')
    ylabel('\omega_2')

    sprintf('k = %d',k)
end
hold off



%% first flip map
% These diagrams map how long it takes for a particular set of initial
% conditions to result in the first pendulum flipping and the second
% pendulum flipping. We condsider a flip to be a rotation past pi in either
% direction. Different conditons could be used to quantitatively determine 
% if a pendulum has flipped.

%constants/initialization
m = [1 1];
L = [1 1];
t_end = 1000;
theta_step = 0.15;
initial_thetas1 = -3:theta_step:3;
initial_thetas2 = 0:theta_step:3;

n = length(initial_thetas1);
p = length(initial_thetas2);
first_flip1 = zeros(n,p);
first_flip2 = zeros(n,p);
counter = 0;

% first flip solutions
for i = 1:n
    for j = 1:p
        theta0 = [initial_thetas1(i); initial_thetas2(j); 0; 0];
        [t_sol, s_sol] = ode45(@(t,s) double_pendulum( t, s, m, L), [0 t_end], theta0, options);
        theta1 = s_sol(:,1);
        theta2 = s_sol(:,2);
        % flip conditions
        flip1 = find(abs(theta1) > pi,1);
        flip2 = find(abs(theta2) > pi,1);
        if ~isempty(flip1)
            first_flip1(i,j) = log10(t_sol(flip1));
        end
        if ~isempty(flip2)
            first_flip2(i,j) = log10(t_sol(flip2));
        end
        counter = counter + 1;
        sprintf('finished trial: %d/%d',counter,n*p)
    end
end

% adjust value of IC's that did not flip so the coloring looks proper
for i = 1:n
    for j = 1:p
        if first_flip1(i,j) == 0
            first_flip1(i,j) = 3;
        end
        if first_flip2(i,j) == 0
            first_flip2(i,j) = 3;
        end
    end
end

figure
subplot(1,2,1)
imagesc(initial_thetas1,initial_thetas2,first_flip1')
xlabel('Initial \theta_1')
ylabel('Initial \theta_2')
title('Time until first flip for pendulum 1')
%axis square
caxis auto
subplot(1,2,2)
imagesc(initial_thetas1,initial_thetas2,first_flip2')
xlabel('Initial \theta_1')
ylabel('Initial \theta_2')
title('Time until first flip for pendulum 2')
%axis square
caxis auto


%% damped pendulum
theta0 = [pi/2; pi/2; 0; 0];
m = [1 1];
L = [1 1];
D = [50, 50];
t_span = 0:0.001:400;
[t_sol, s_sol] = ode45(@(t,s) dd_pendulum( t, s, m, L, D), t_span, theta0, options);
theta1 = s_sol(:,1);
theta2 = s_sol(:,2);
theta1_dot=s_sol(:,3);
theta2_dot=s_sol(:,4);

x1 = L(1)*sin(theta1);
x2 = L(1)*sin(theta1) + L(2)*sin(theta2);
y1 = -L(1)*cos(theta1);
y2 = -L(1)*cos(theta1) - L(2)*cos(theta2);

hold on
plot(t_sol(1:30000),theta1(1:30000))
plot(t_sol(1:30000),theta2(1:30000))
hold off

mov = VideoWriter('damped_d1_0_d2_0.5.avi', 'Uncompressed AVI');
open(mov);
% damped driven pendulum movie
for i = 1:200:200*400
    clf
    hold on
    plot([x1(i), x2(i)],[y1(i), y2(i)], 'k')
    plot([x1(i), 0],[y1(i), 0], 'k')
    plot(x1(max(1,i-36000):i),y1(max(1,i-36000):i),'r')
    plot(x2(max(1,i-36000):i),y2(max(1,i-36000):i),'b')
    axis([-3,3,-3,3])
    axis square
    xlabel('x position (m)')
    ylabel('y position (m)')
    title('Damped Double Pendulum')
    hold off
    M = getframe(gcf);
    writeVideo(mov, M);
end
close(mov);
hold on
plot (t_span(1:round(length(t_span)/8)),x1(1:round(length(t_span)/8)))
title ('Double Pendulum: drag constants: 50, 50')
xlabel('timestep')
ylabel('x1 (m)')
plot (t_span(1:round(length(t_span)/8)),x2(1:round(length(t_span)/8)))
title ('Double Pendulum: drag constants: 50, 50')
xlabel('timestep')
ylabel('x position (m)')
legend('x1','x2')
hold off

hold on
plot (t_span(1:round(length(t_span)/30)),x1(1:round(length(t_span)/30)))
title ('Double Pendulum: drag constants: 100, 100')
xlabel('timestep')
ylabel('x1 (m)')
plot (t_span(1:round(length(t_span)/30)),x2(1:round(length(t_span)/30)))
title ('Double Pendulum: drag constants: 100, 100')
xlabel('timestep')
ylabel('x position (m)')
legend('x1','x2')
hold off
