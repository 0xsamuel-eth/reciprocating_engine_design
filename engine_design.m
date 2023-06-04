%% Problem Setup
 
%% Given Values
 
m1 = 150; %lbm
l = 100; %in
a = 20; %in
t = 0.5; %in
E = 30e6; %psi
c1 = 1; %lbf-s/ft
F0 = 50; %lbf
%wop = 50:60; %rad/s (range)
 
%% Conversions
 
m1 = m1/32.2; %lbm to slugs
l = l/12; %in to ft
a = a/12; %in to ft
t = t/12; %in to ft
E = E*144; %psi to lbf/ft^2
 
%% Part I #1
 
 
% Find k
k1 = 16*E*a*(t^3)/(l^3); %lbf/ft
 
% Frequency Response
 
wop = 1:.001:100;
xtilde1 = zeros(1,100);
x1 = 1:.001:100;
counter = 1;
 
for i = 1:.001:100
    
    xtilde1(counter) = F0./sqrt( (k1-m1*(i.^2)).^2+(c1*i).^2 );
    
    counter = counter + 1;
 
end
 
% Phase Angle
 
phi1 = zeros(1, 100);
counter = 1;
 
for i = 1:100
    
    phi1(counter) = atan( (c1*i)/(k1-m1*(i.^2)) );
    
    if (k1-m1*(i.^2)) < 0
       phi1(i) = phi1(i) + pi;
    end
    counter = counter + 1;
    
end  
 
%% Plots
 
% Phase Angle
 
x1_phi = 1:100;
 
figure(1);
plot(x1_phi, phi1, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Phase Angle {\phi} (rad)');
title('Phase Angle vs. Frequency');
grid on;
 
%% Part I #2
 
% Given Values
w = 55.6; %rad/s
time = linspace(0,10*pi/w,100); %sec
 
% Steady-State Response
xtilde1_state = F0/sqrt( (k1-m1*(w^2))^2+(c1*w)^2 );
phi1_state = atan( (c1*w)/(k1-m1*(w^2)) );
xp1 = xtilde1_state * exp(1i*(w*time-phi1_state));
 
 
%% Plot
 
% Steady State
 
figure(2);
plot(time, xp1, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Time (sec)');
ylabel('Steady State x_{p}(t) (m)');
title('Steady State Response vs. Time');
grid on;
 
%% Part II 
 
% Givens
 
w = 55.6; %rad/sec
m2 = 20/32.2; %lbm to slugs
k2 = m2*(w^2); %lbf/ft
 
% Impedance Matrix
 
%Z11 =  -m1*(w^2) + k1 + k2 + 1i*w*c1;
%Z12 = -k2;
%Z21 = Z12;
%Z22 = -m2*(w^2) + k2;
 
% Xtilde
 
%xtilde21 = ( Z22*F0 ) / ( Z11*Z22 - (Z12^2) );
%xtilde22 = ( -Z21*F0 ) / ( Z11*Z22 - (Z12^2) );
 
 
wop = 1:.001:100;
xtilde21 = zeros(1,100);
xtilde22 = zeros(1,100);
x1 = 1:.001:100;
counter = 1;
 
for i = 1:.001:100
    
    Z11 =  -m1*(i.^2) + k1 + k2 + 1i*i*c1;
    Z12 = -k2;
    Z21 = Z12;
    Z22 = -m2*(i.^2) + k2;
    
    xtilde21(counter) = ( Z22*F0 ) / ( Z11*Z22 - (Z12^2) );
    xtilde22(counter) = ( -Z21*F0 ) / ( Z11*Z22 - (Z12^2) );
 
    counter = counter +1;
 
end
 
xtilde21 = abs(xtilde21);
xtilde22 = abs(xtilde22);
 
%% Part III
 
% Givens
c2 = 5; %lbf-s/ft
 
% Impedance Matrix
 
%Z11 =  -m1*(w^2) + k1 + k2 + 1i*w*(c1+c2);
%Z12 = -k2 - 1i*c2;
%Z21 = Z12;
%Z22 = -m2*(w^2) + k2 + 1i*w*c2;
 
% Xtilde
 
%xtilde31 = ( Z22*F0 ) / ( Z11*Z22 - (Z12^2) );
%xtilde32 = ( -Z21*F0 ) / ( Z11*Z22 - (Z12^2) );
 
wop = 1:.001:100;
xtilde31 = zeros(1,100);
xtilde32 = zeros(1,100);
x1 = 1:.001:100;
counter = 1;
 
for i = 1:.001:100
    
    Z11 =  -m1*(i.^2) + k1 + k2 + 1i*i*(c1+c2);
    Z12 = -k2 - 1i*c2;
    Z21 = Z12;
    Z22 = -m2*(i.^2) + k2 + 1i*i*c2;
 
    xtilde31(counter) = ( Z22*F0 ) / ( Z11*Z22 - (Z12^2) );
    xtilde32(counter) = ( -Z21*F0 ) / ( Z11*Z22 - (Z12^2) );
    
    counter = counter + 1;
end
 
xtilde31 = abs(xtilde31);
xtilde32 = abs(xtilde32);
 
%% Frequency Spectra of Systems (Fig.7)
 
figure(3);
 
subplot(3,1,1)
plot(x1,xtilde1, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Frequency Spectra - Original System Engine');
ylim([0,1]);
grid on;
 
subplot(3,1,2)
plot(x1,xtilde21, 'r','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Frequency Spectra - Undamped System Engine');
ylim([0,1.2]);
grid on;
 
subplot(3,1,3)
plot(x1,xtilde31, 'm','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Frequency Spectra - Damped System Engine');
ylim([0,0.03]);
grid on;
 
%% System Amps Separate
 
figure(4);
plot(x1,xtilde1, 'b','LineWidth',2,'MarkerSize', 15);
hold on;
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Original System Frequency Spectra');
ylim([0,1.5]);
grid on;
 
figure(5);
plot(x1,xtilde21, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Undamped System Frequency Spectra');
ylim([0,1.5]);
grid on;
 
figure(6);
plot(x1,xtilde31, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Damped System Frequency Spectra');
ylim([0,1.5]);
grid on;
 
%% Damped & Un-damped Plots (Fig.8)
 
figure(7);
 
subplot(2,1,1)
plot(x1,xtilde22, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Undamped Absorber');
ylim([0,4]);
grid on;
 
subplot(2,1,2)
plot(x1,xtilde32, 'r','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Damped Absorber');
ylim([0,0.07]);
grid on;
 
%% Absorber Amps Separate
 
figure(8);
plot(x1,xtilde22, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Undamped Absorber Frequency Spectra');
ylim([0,1.5]);
grid on;
 
figure(9);
plot(x1,xtilde32, 'b','LineWidth',2,'MarkerSize', 15);
xlabel('Frequency (rad/sec)');
ylabel('Amplitude (m)');
title('Damped Absorber Frequency Spectra');
ylim([0,1.5]);
grid on;
