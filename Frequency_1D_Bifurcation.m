 clear all
close all
clc
dt=1/256; % Intergation step
time=60;    % Simulation time
time_length=time/dt;
trails=1;
cpy_py=1.89;
 cpy_i1=4;
 ci1_py=1.8;
 cre_re=0.01;
 ctc_re=10;
 cre_tc=1.4;
 cpy_tc=3;
 cpy_re=1.4;
 ctc_py=1;
 cpy_i2=1.5;
 ci2_i1=0.1;
 ci1_i2=0.5;
 ci2_py=0.05;
 ctc_i1=0.05;
 ctc_i2=0.05;
 ctc_ei=4.5;
 cei_py=0.442;
 cei_i1=0.05;
 cpy_ei=0.76;
 ci1_ei=0.3;
 % Time scale parameters
   tau1=21.5;
   tau2=31.5;
   tau3=0.1;
   tau4=3.8;
   tau5=3.9;
   tau6=4.5;
% Input parameters 
   h_py=-0.4;
   h_i1=-3.4;
   h_i2=-2;
   h_tc=-2.5;
   h_re=-3.2;
   h_ei=-1;
   sigma=250000;
   miangin_py=0.7;
   miangin_tc=0.1;
%initial conditions
   
   PY=0.2775;
   I1=  0.5345;
   I2= -1.0365;
   TC=-0.1216;
   RE=-0.0401;
   EI=0.2888;
%4th-order Runge-kutta
f_tc=linspace(0,30,800);
cortical_main=nan((time/dt)+1,length(f_tc));
f_py=0;
for j = 1 : length(f_tc)
   j
t=linspace(0,time,time_length);
for i=1:time/dt
    a_tc=0.02;
    a_py=0;
    dw_py=a_py*sin(2*pi*t(i)*f_py(j))+miangin_py;
    dw_tc=a_tc*sin(2*pi*t(i)*f_tc(j))+miangin_tc;
        %%First order 
        k11=(h_py -(PY(i)) +cpy_py*(1./(1+sigma.^(-(PY(i))))) -ci1_py*(1./(1+sigma^(-I1(i)))) -ci2_py*(1./(1+sigma^(-I2(i)))) + ctc_py*(1./(1+sigma^(-TC(i))))+cei_py*(1./(1+sigma.^(-(EI(i))))))*tau1+dw_py;
        k21=(h_i1 -I1(i) +cpy_i1*(1./(1+sigma.^(-(PY(i))))) -ci2_i1*(1./(1+sigma^(-I2(i)))) +ctc_i1*(1./(1+sigma^(-TC(i))))+cei_i1*(1./(1+sigma.^(-(EI(i))))))*tau2;
        k31=(h_i2 -I2(i) +cpy_i2*(1./(1+sigma.^(-(PY(i))))) -ci1_i2*(1./(1+sigma^(-I1(i)))) +ctc_i2*(1./(1+sigma^(-TC(i)))))*tau3; 
        k41=(h_tc -TC(i) +cpy_tc*(1./(1+sigma.^(-(PY(i))))) -cre_tc*(1./(1+sigma^(-RE(i)))))*tau4+dw_tc;
        k51=(h_re -RE(i) +cpy_re*(1./(1+sigma.^(-(PY(i))))) + ctc_re*(1./(1+sigma^(-TC(i)))) -cre_re*(1./(1+sigma^(-RE(i)))))*tau5;
        k61=(h_ei-EI(i)+cpy_ei*(1./(1+sigma.^(-(PY(i)))))+ctc_ei*(1./(1+sigma.^(-(TC(i)))))-ci1_ei*(1./(1+sigma^(-I1(i)))))*tau6;

        
        %%Second order 
        k12=(h_py -(PY(i)+dt*k11/2) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k21/2)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k31/2)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k41/2))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k61/2)))))*tau1+dw_py;
        k22=(h_i1 -(I1(i)+dt*k21/2) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k31/2))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k41/2))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k61/2)))))*tau2;
        k32=(h_i2 -(I2(i)+dt*k31/2) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k21/2))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k41/2)))))*tau3; 
        k42=(h_tc -(TC(i)+dt*k41/2) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -cre_tc*(1./(1+sigma^(-(RE(i)+dt*k51/2)))))*tau4+dw_tc;
        k52=(h_re -(RE(i)+dt*k51/2) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k41/2)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k51/2)))))*tau5;
        k62=(h_ei-(EI(i)+dt*k61/2)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k11/2))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k41/2))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k21/2)))))*tau6;
        
        %%Third order
        k13=(h_py -(PY(i)+dt*k12/2) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k22/2)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k32/2)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k42/2))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k62/2)))))*tau1+dw_py;
        k23=(h_i1 -(I1(i)+dt*k22/2) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k32/2))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k42/2))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k62/2)))))*tau2;
        k33=(h_i2 -(I2(i)+dt*k32/2) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k22/2))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k42/2)))))*tau3; 
        k43=(h_tc -(TC(i)+dt*k42/2) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -cre_tc*(1./(1+sigma^(-(RE(i)+dt*k52/2)))))*tau4+dw_tc;
        k53=(h_re -(RE(i)+dt*k52/2) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k42/2)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k52/2)))))*tau5;
        k63=(h_ei-(EI(i)+dt*k62/2)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k12/2))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k42/2))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k22/2)))))*tau6;
       
        %%Fourth order
        k14=(h_py -(PY(i)+dt*k13) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k13)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k23)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k33)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k43))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k63)))))*tau1+dw_py;
        k24=(h_i1 -(I1(i)+dt*k23) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k13)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k33))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k43))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k63)))))*tau2;
        k34=(h_i2 -(I2(i)+dt*k33) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k13))))-ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k23))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k43)))))*tau3; 
        k44=(h_tc -(TC(i)+dt*k43) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k13))))-cre_tc*(1./(1+sigma^(-(RE(i)+dt*k53)))))*tau4+dw_tc;
        k54=(h_re -(RE(i)+dt*k53) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k13)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k43)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k53)))))*tau5;     
        k64=(h_ei-(EI(i)+dt*k63)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k13))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k43))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k23)))))*tau6;
        
        PY(i+1)=PY(i)+dt*(k11+2*k12+2*k13+k14)/6;
        I1(i+1)=I1(i)+dt*(k21+2*k22+2*k23+k24)/6;
        I2(i+1)=I2(i)+dt*(k31+2*k32+2*k33+k34)/6;
        TC(i+1)=TC(i)+dt*(k41+2*k42+2*k43+k44)/6;
        RE(i+1)=RE(i)+dt*(k51+2*k52+2*k53+k54)/6;
        EI(i+1)=EI(i)+dt*(k61+2*k62+2*k63+k64)/6;
 
   
     
end 
 %Output is the mean of cortical populations
cortical_main2=(PY+I1+I2+EI)/4;
%removing DC 
cortical_main(:,j)=cortical_main2-mean(cortical_main2);
end

TD=58.58; % starting time used for analysis (The last two seconds)
peaksmax = nan(size(cortical_main));
peaksmin = nan(size(cortical_main));
DF  = nan(1,length(f_tc));
pmax1=nan(1,length(f_tc));
pmax2=nan(1,length(f_tc));
pmin1=nan(1,length(f_tc));
pmin2=nan(1,length(f_tc));
state=nan(1,length(f_tc));
eq  = nan(1,length(f_tc));
for j = 1 : length(f_tc)
 j
uxcut=cortical_main(TD/dt+1:end,j);
x=findpeaks(uxcut);
peaksmax(1 : length(x) , j) = x;
peaksmax(:,j)=sort(peaksmax(:,j));
y=findpeaks(-uxcut);
peaksmin(1 : length(y) , j) = y;
peaksmin(:,j)=sort(-peaksmin(:,j));

if (length(peaksmax(:,j))>=1)
    pmax1(1,j)=peaksmax(end,j);
    pmax2(1,j)=peaksmax(1,j);
    index=find(pmax1(1,j)-peaksmax(:,j)>=0.001);
    if (length(index)>=1)
        pmax2(1,j)=peaksmax(index(end));
    end
else
    pmax1(1,j)=uxcut(1);
    pmax2(1,j)=uxcut(2);
end


if (length(peaksmin(:,j))>=1)
    pmin1(1,j)=peaksmin(end,j);
    pmin2(1,j)=peaksmin(1,j);
else
    pmin1(1,j)=uxcut(1);
    pmin2(1,j)=uxcut(2);
end  
  
end 
M_size = 5;
peaksmax=100*peaksmax;
figure(1)
plot(f_tc, peaksmax , '.b','linewidth',1.5)
title('(b)','FontSize',16)
set(gca,'FontSize',14)
xlabel('ftc (Hz)','FontSize',16)
ylabel('Extrema mean PY,I1,I2 and EI (mv)','FontSize',16)

