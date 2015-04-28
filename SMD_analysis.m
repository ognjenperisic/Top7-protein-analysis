cd c:\projects\colaboration\QYS\QYS_pulling_072406
close all
clear all
clc

fig=0;

SMD=load('smd_072606.txt');

hb_all=load('nmb_hb_all.dat');

hb_first  = load('nmb_hb_first.dat');
hb_second = load('nmb_hb_second.dat');
hb_third  = load('nmb_hb_third.dat');
hb_fourth = load('nmb_hb_fourth.dat');

wm_first = load('nmb_wm_r.dat');
wm_first=wm_first-wm_first(1);


%SMD=smd_062206;

nx=0.99;
ny=-0.06;
nz=0.02;
vel=0.0001;
n_tmstp=length(SMD);  % number of time steps
position=zeros(length(SMD),3); % positions

step=zeros(1,n_tmstp);
for i=1:n_tmstp
   step(i)=SMD(i,1)*vel;
   %position(i,1)=SMD(i,2);
   %position(i,2)=SMD(i,3);
   %position(i,3)=SMD(i,4);
end


force=nx*SMD(:,5)+ny*SMD(:,6)+nz*SMD(:,7);

fig=fig+1;
figure(fig)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
grid on
axis([-inf inf 1.5*min(force) 1.3*max(force)])

%figure(2)
%plot3(position(:,1),position(:,2),position(:,3))
%grid on

stepk=step(n_tmstp)/length(hb_all);

i=1:length(hb_all);
step2=i*stepk;

fig=fig+1;
figure(fig)
plot(step2,hb_first)
hold on
plot(step2,hb_second,'g')
hold on
plot(step2,hb_third,'r','Linewidth',2)
hold on
plot(step2,hb_fourth,'c')
hold on
plot(step2,hb_all,'k')
legend('first','second','third','fourth','all',1)
xlabel('extension (A)')
ylabel('hudrogen bonds - Beta sheet')
grid on
axis([-inf inf 0 1.3*max(hb_all)])

fig=fig+1;
figure(fig)
plot(step2,wm_first)
axis([-inf inf 0 1.3*max(wm_first)])
title('number of water molecules around residues 10 29 33 36 40 50 52 63 67 71')
xlabel('extension (A)')
grid on

fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,hb_first)
hold on
plot(step2,hb_second,'g')
hold on
plot(step2,hb_third,'r')
hold on
plot(step2,hb_fourth,'c')
hold on
plot(step2,hb_all,'k')
legend('first','second','third','fourth','all',1)
xlabel('extension (A)')
ylabel('hudrogen bonds - Beta sheet')
grid on
axis([-inf inf 0 1.3*max(hb_all)])


fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,wm_first)
axis([-inf inf 0 1.3*max(wm_first)])
xlabel('extension (A)')
title('number of water molecules around residues 10 29 33 36 40 50 52 63 67 71')
grid on