% 1QYS pulling  0.0001   <- 

cd c:\projects\colaboration\QYS\QYS_pulling_080406
close all
clear all
clc

fig=0;

SMD=load('smd_092106.txt');

hb_all    = load('nmb_hb_all.dat');
hb_first  = load('nmb_hb_first.dat');
hb_second = load('nmb_hb_second.dat');
hb_third  = load('nmb_hb_third.dat');
hb_fourth = load('nmb_hb_fourth.dat');
wm_first  = load('nmb_wm_r.dat');

sb1       = load('nmb_sb_1.dat');
sb2       = load('nmb_sb_2.dat');

helix1    = load('helix1.dat');
helix2    = load('helix2.dat');


if (2>3)   
   hb_all_r1    = load('nmb_hb_allr1.dat');
   hb_first_r1  = load('nmb_hb_firstr1.dat');
   hb_second_r1 = load('nmb_hb_secondr1.dat');
   hb_third_r1  = load('nmb_hb_thirdr1.dat');
   hb_fourth_r1 = load('nmb_hb_fourthr1.dat');
   wm_first_r1  = load('nmb_wm_rr1.dat');
   sb1_r1       = load('nmb_sb_1_r1.dat');
   sb2_r1       = load('nmb_sb_2_r1.dat');
   
   hb_all_r2    = load('nmb_hb_allr2.dat');
   hb_first_r2  = load('nmb_hb_firstr2.dat');
   hb_second_r2 = load('nmb_hb_secondr2.dat');
   hb_third_r2  = load('nmb_hb_thirdr2.dat');
   hb_fourth_r2 = load('nmb_hb_fourthr2.dat');
   wm_first_r2  = load('nmb_wm_rr2.dat');
   sb1_r2       = load('nmb_sb_1_r2.dat');
   sb2_r2       = load('nmb_sb_2_r2.dat');
   
   hb_all    = [hb_all' hb_all_r1' hb_all_r2'];
   hb_first  = [hb_first' hb_first_r1' hb_first_r2'];
   hb_second = [hb_second' hb_second_r1' hb_second_r2'];
   hb_third  = [hb_third' hb_third_r1' hb_third_r2'];
   hb_fourth = [hb_fourth' hb_fourth_r1' hb_fourth_r2'];
   wm_first  = [wm_first' wm_first_r1' wm_first_r2'];
   sb1       = [sb1' sb1_r1' sb1_r2'];
   sb2       = [sb2' sb2_r1' sb2_r2'];
   
end


%sb1       = load('nmb_sb_1.dat');
%sb2       = load('nmb_sb_2.dat');
kratak=0;

if (kratak)
   SMD1=zeros(1000,7);
   hb_all    =  hb_all(1:1000);
   hb_first  =  hb_first(1:1000);
   hb_second =  hb_second(1:1000);
   hb_third  =  hb_third(1:1000);
   hb_fourth =  hb_fourth(1:1000);
   wm_first  =  wm_first(1:1000);
   sb1       =  sb1(1:1000);
   sb2       =  sb2(1:1000);
   SMD1(:,1)  =  SMD(1:1000,1);
   SMD1(:,2)  =  SMD(1:1000,2);
   SMD1(:,3)  =  SMD(1:1000,3);
   SMD1(:,4)  =  SMD(1:1000,4);
   SMD1(:,5)  =  SMD(1:1000,5);
   SMD1(:,6)  =  SMD(1:1000,6);  
   SMD1(:,7)  =  SMD(1:1000,7);
   SMD=SMD1;
end



wm_first=wm_first-wm_first(1);


%SMDVel 0.0001
%SMDDir -0.991 0.095 -0.083

nx=1.00;
ny=0.00;
nz=0.00;
vel=0.00005;
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

granica=10;
lnf=length(force);
smooth_force=zeros(lnf,1);

for i=1:lnf
   if (i<=granica)
      smooth_force(i)=mean(force(1:i+granica));
   else
      if (i>=lnf-granica)
         smooth_force(i)=mean(force(i-granica:lnf));
      else
         smooth_force(i)=mean(force(i-granica:i+granica));
      end
   end
end


granica_wm=10;
lnf=length(wm_first);
smooth_wm=zeros(lnf,1);

for i=1:lnf
   if (i<=granica_wm)
      smooth_wm(i)=mean(wm_first(1:i+granica_wm));
   else
      if (i>=lnf-granica_wm)
         smooth_wm(i)=mean(wm_first(i-granica_wm:lnf));
      else
         smooth_wm(i)=mean(wm_first(i-granica_wm:i+granica_wm));
      end
   end
end

granica_wm=10;
lnf=length(wm_first);
smooth_hb_all    = zeros(lnf,1);
smooth_hb_first  = zeros(lnf,1);
smooth_hb_second = zeros(lnf,1);
smooth_hb_third  = zeros(lnf,1);
smooth_hb_fourth = zeros(lnf,1);
smooth_wm_first  = zeros(lnf,1);
smooth_helix1    = zeros(lnf,1);
smooth_helix2    = zeros(lnf,1);

for i=1:lnf
   if (i<=granica_wm)
      %smooth_wm(i)=mean(wm_first(1:i+granica_wm));
      smooth_hb_all(i)   = mean(hb_all(1:i+granica_wm));
      smooth_hb_first(i) = mean(hb_first(1:i+granica_wm));
      smooth_hb_second(i)= mean(hb_second(1:i+granica_wm));
      smooth_hb_third(i) = mean(hb_third(1:i+granica_wm));
      smooth_hb_fourth(i)= mean(hb_fourth(1:i+granica_wm));
      smooth_wm_first(i) = mean(wm_first(1:i+granica_wm));      
      smooth_helix1(i)   = mean(helix1(1:i+granica_wm));
      smooth_helix2(i)   = mean(helix2(1:i+granica_wm));
   else
      if (i>=lnf-granica_wm)
         %smooth_wm(i)=mean(wm_first(i-granica_wm:lnf));
         smooth_hb_all(i)    = mean(hb_all(i-granica_wm:lnf));
         smooth_hb_first(i)  = mean(hb_first(i-granica_wm:lnf));
         smooth_hb_second(i) = mean(hb_second(i-granica_wm:lnf));
         smooth_hb_third(i)  = mean(hb_third(i-granica_wm:lnf));
         smooth_hb_fourth(i) = mean(hb_fourth(i-granica_wm:lnf));
         smooth_wm_first(i)  = mean(wm_first(i-granica_wm:lnf));
         smooth_helix1(i)    = mean(helix1(i-granica_wm:lnf)); 
         smooth_helix2(i)    = mean(helix2(i-granica_wm:lnf));
      else
         %smooth_wm(i)=mean(wm_first(i-granica_wm:i+granica_wm));
         smooth_hb_all(i)    = mean(hb_all(i-granica_wm:i+granica_wm));
         smooth_hb_first(i)  = mean(hb_first(i-granica_wm:i+granica_wm));
         smooth_hb_second(i) = mean(hb_second(i-granica_wm:i+granica_wm));
         smooth_hb_third(i)  = mean(hb_third(i-granica_wm:i+granica_wm));
         smooth_hb_fourth(i) = mean(hb_fourth(i-granica_wm:i+granica_wm));
         smooth_wm_first(i)  = mean(wm_first(i-granica_wm:i+granica_wm));
         smooth_helix1(i)    = mean(helix1(i-granica_wm:i+granica_wm));
         smooth_helix2(i)    = mean(helix2(i-granica_wm:i+granica_wm));
      end
   end
end





fig=fig+1;
figure(fig)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06-new')
grid on
axis([-inf inf 1.5*min(force) 1.3*max(force)])

%figure(2)
%plot3(position(:,1),position(:,2),position(:,3))
%grid on

stepk=step(n_tmstp)/length(hb_all);

i=1:length(hb_all);
step2=i*stepk;

%step2=distance;

fig=fig+1;
figure(fig)
plot(step2,hb_first)
hold on
plot(step2,hb_second,'g','Linewidth',2)
hold on
plot(step2,hb_third,'r','Linewidth',2,'Linewidth',2)
hold on
plot(step2,hb_fourth,'c')
hold on
plot(step2,hb_all,'k')
legend('first','second','third','fourth','all',1)
xlabel('extension (A)')
ylabel('hudrogen bonds - Beta sheet')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
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
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,hb_first)
hold on
plot(step2,hb_second,'g','Linewidth',2)
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
subplot(2,1,1)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,sb1)
axis([-inf inf -inf inf])
title('number of salt bridges')
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
plot(step2,wm_first)
axis([-inf inf 0 1.3*max(wm_first)])
xlabel('extension (A)')
title('number of water molecules around residues 10 29 33 36 40 50 52 63 67 71')
grid on


fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,helix1);
xlabel('extension (A)')
ylabel('number of hydrogen bonds - helix1');
axis([-inf inf 0.95*min(helix1) 1.05*max(helix1)])
grid on



fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,force,'k');
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,helix2);
xlabel('extension (A)')
ylabel('number of hydrogen bonds - helix2');
axis([-inf inf 0.95*min(helix2) 1.05*max(helix2)])
grid on




fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,smooth_force,'k')
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])

subplot(2,1,2)
plot(step2,hb_first)
hold on
plot(step2,hb_second,'g','Linewidth',2)
hold on
plot(step2,hb_third,'r','Linewidth',2,'Linewidth',2)
hold on
plot(step2,hb_fourth,'c')
hold on
plot(step2,hb_all,'k')
legend('first','second','third','fourth','all',1)
xlabel('extension (A)')
ylabel('hydrogen bonds - Beta sheet')
grid on
axis([-inf inf 0 1.3*max(hb_all)])


fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,smooth_force,'k')
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])
title('Pulling from C to N terminal, pulling velocity 0.0001A/ps');

subplot(2,1,2)
plot(step2,smooth_hb_first)
hold on
plot(step2,smooth_hb_second,'g','Linewidth',2)
hold on
plot(step2,smooth_hb_third,'r','Linewidth',2,'Linewidth',2)
hold on
plot(step2,smooth_hb_fourth,'c')
hold on
plot(step2,smooth_hb_all,'k')
legend('first','second','third','fourth','all',1)
xlabel('extension (A)')
ylabel('hydrogen bonds - Beta sheet')
grid on
axis([-inf inf 0 1.3*max(hb_all)])


fig=fig+1;
figure(fig)
subplot(2,1,1)
plot(step,smooth_force,'k')
xlabel('extension (A)')
ylabel('Force (pN)')
title('Pulling velocity 0.00005 A/fs     date 08.04.06')
grid on
axis([-inf inf 1.5*min(force) 1.1*max(force)])
title('Pulling from C to N terminal, pulling velocity 0.0001A/ps');
subplot(2,1,2)
plot(step2,smooth_helix1,'b')
xlabel('extension (A)')
ylabel('hydrogen bonds - helix1')
grid on
axis([-inf inf 0 1.3*max(smooth_helix1)])

