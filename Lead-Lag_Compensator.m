%Gs = 31.98/(s*(0.096*s+1));
GPole1 = -1/0.096;
numOfGs = 31.98/0.096;
Gpole2 = 0;

%% changable params
PO = 10; %% in percent
Ts = 0.3; %% in sec


%% find loc of pole
Z = -log(PO/100)/(sqrt(pi^2+log(PO/100)^2));
W = 4/(Ts*Z);

poleRe = -W*Z;
poleIM = W*(sqrt(1-Z^2));

% put in complex num
desired_pole = poleRe + 1i * poleIM;

% find angle
theta1 = 180 - atand(abs(poleIM/poleRe));

distfromPole = abs(poleRe-GPole1); % distance from desired pole RE to natural system pole

theta2 = 180 - atand(abs(poleIM/distfromPole));
thetac = 180+theta1+theta2-360;

thetap = 90-thetac;

intermediateRe = poleIM/tand(thetap);


%% NEW pole for controller that allows for desired root locus
compensatorPole = -(poleRe-intermediateRe);
compensatorZero = -poleRe;
%% find transfer function
s = tf('s');

% magnitude critera
dCompPoletodRe = compensatorPole-(-poleRe); % the distance from compensator pole to pole re
h1 = sqrt(poleIM^2+dCompPoletodRe^2);
h2 = sqrt(poleRe^2+poleIM^2);
h3 = sqrt(distfromPole^2+poleIM^2);

Kc = h1*h2*h3/(poleIM*numOfGs);


%% find transfer function

% Lead only
Gc = (s+compensatorZero)/(s+compensatorPole);
G = 31.98/(s*(0.096*s+1));
openLoop = Gc*G*Kc;
H = 1;

% Lead-lag
Gc2 = ((s+compensatorZero)*(s+0.1))/((s+compensatorPole)*(s+0.01));
G2 = 31.98/(s*(0.096*s+1));
openLoop2 = Gc2*G*Kc;

figure(1)
rlocus(openLoop);
hold on;
title('Lead Compensated Root Locus');

figure(2)
rlocus(openLoop2)
hold on
title('Lead-lag Compensated Root Locus');


sys_cl = feedback(openLoop,H);
sys_cl_min = minreal(sys_cl); % remove zero/poles that cancel

sys_cl2 = feedback(openLoop2,H);
sys_cl_min2 = minreal(sys_cl2);


stepinfo(sys_cl_min);
stepinfo(sys_cl_min2);

% plot step response of CL system
figure(3)
step(sys_cl_min, sys_cl_min2)
title('Step Response of Lead vs Lead-lag')
legend('Lead','Lead-lag')
grid on;

% Plot Ramp Response

% Find slowest pole
poles_lead = pole(sys_cl_min);
poles_leadlag = pole(sys_cl_min2);
all_poles = [poles_lead; poles_leadlag];

% Time constant
tau_slowest = 1 / min(abs(real(all_poles)));
t_end = 6 * tau_slowest; % Six time constants

% r(6t)
t = 0:0.01:t_end; 
r = 6 * t;

% Plot
y_lead_ramp = lsim(sys_cl_min, r, t);
y_leadlag_ramp = lsim(sys_cl_min2, r, t);
figure(4)
plot(t, r, 'k--', 'LineWidth', 1.5) 
hold on
plot(t, y_lead_ramp, 'b', 'LineWidth', 1.5)
plot(t, y_leadlag_ramp, 'r', 'LineWidth', 1.5)
hold off
title('Ramp Response of Lead vs Lead-lag')
xlabel('Time (seconds)')
ylabel('Amplitude')
legend('Reference', 'Lead', 'Lead-lag')
grid on;
