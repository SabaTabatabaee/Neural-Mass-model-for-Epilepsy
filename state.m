function [pmax1,pmax2,pmin1,pmin2,state,fd]=state( cpy_ei,ci1_ei)

dt=1/256;   % Intergation step
time=80;    % Simulation time

%connectivity parameters
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
   %%% 4th-order Runge-kutta
for i=1:time/dt
    
        %%First order
        k11=(h_py -(PY(i)) +cpy_py*(1./(1+sigma.^(-(PY(i))))) -ci1_py*(1./(1+sigma^(-I1(i)))) -ci2_py*(1./(1+sigma^(-I2(i)))) + ctc_py*(1./(1+sigma^(-TC(i))))+cei_py*(1./(1+sigma.^(-(EI(i))))))*tau1+miangin_py;
        k21=(h_i1 -I1(i) +cpy_i1*(1./(1+sigma.^(-(PY(i))))) -ci2_i1*(1./(1+sigma^(-I2(i)))) +ctc_i1*(1./(1+sigma^(-TC(i))))+cei_i1*(1./(1+sigma.^(-(EI(i))))))*tau2;
        k31=(h_i2 -I2(i) +cpy_i2*(1./(1+sigma.^(-(PY(i))))) -ci1_i2*(1./(1+sigma^(-I1(i)))) +ctc_i2*(1./(1+sigma^(-TC(i)))))*tau3; 
        k41=(h_tc -TC(i) +cpy_tc*(1./(1+sigma.^(-(PY(i))))) -cre_tc*(1./(1+sigma^(-RE(i)))))*tau4+miangin_tc;
        k51=(h_re -RE(i) +cpy_re*(1./(1+sigma.^(-(PY(i))))) + ctc_re*(1./(1+sigma^(-TC(i)))) -cre_re*(1./(1+sigma^(-RE(i)))))*tau5;
        k61=(h_ei-EI(i)+cpy_ei*(1./(1+sigma.^(-(PY(i)))))+ctc_ei*(1./(1+sigma.^(-(TC(i)))))-ci1_ei*(1./(1+sigma^(-I1(i)))))*tau6;

        
        %%Second order
        k12=(h_py -(PY(i)+dt*k11/2) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k21/2)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k31/2)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k41/2))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k61/2)))))*tau1+miangin_py;
        k22=(h_i1 -(I1(i)+dt*k21/2) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k31/2))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k41/2))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k61/2)))))*tau2;
        k32=(h_i2 -(I2(i)+dt*k31/2) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k21/2))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k41/2)))))*tau3; 
        k42=(h_tc -(TC(i)+dt*k41/2) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) -cre_tc*(1./(1+sigma^(-(RE(i)+dt*k51/2)))))*tau4+miangin_tc;
        k52=(h_re -(RE(i)+dt*k51/2) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k11/2)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k41/2)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k51/2)))))*tau5;
        k62=(h_ei-(EI(i)+dt*k61/2)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k11/2))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k41/2))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k21/2)))))*tau6;
        
        %%Third order
        k13=(h_py -(PY(i)+dt*k12/2) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k22/2)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k32/2)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k42/2))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k62/2)))))*tau1+miangin_py;
        k23=(h_i1 -(I1(i)+dt*k22/2) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k32/2))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k42/2))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k62/2)))))*tau2;
        k33=(h_i2 -(I2(i)+dt*k32/2) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k22/2))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k42/2)))))*tau3; 
        k43=(h_tc -(TC(i)+dt*k42/2) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) -cre_tc*(1./(1+sigma^(-(RE(i)+dt*k52/2)))))*tau4+miangin_tc;
        k53=(h_re -(RE(i)+dt*k52/2) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k12/2)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k42/2)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k52/2)))))*tau5;
        k63=(h_ei-(EI(i)+dt*k62/2)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k12/2))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k42/2))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k22/2)))))*tau6;
       
        %%Fourth order
        k14=(h_py -(PY(i)+dt*k13) +cpy_py*(1./(1+sigma^(-(PY(i)+dt*k13)))) -ci1_py*(1./(1+sigma^(-(I1(i)+dt*k23)))) -ci2_py*(1./(1+sigma^(-(I2(i)+dt*k33)))) + ctc_py*(1./(1+sigma^(-(TC(i)+dt*k43))))+cei_py*(1./(1+sigma.^(-(EI(i)+dt*k63)))))*tau1+miangin_py;
        k24=(h_i1 -(I1(i)+dt*k23) +cpy_i1*(1./(1+sigma^(-(PY(i)+dt*k13)))) -ci2_i1*(1./(1+sigma^(-(I2(i)+dt*k33))))+ctc_i1*(1./(1+sigma^(-(TC(i)+dt*k43))))+cei_i1*(1./(1+sigma.^(-(EI(i)+dt*k63)))))*tau2;
        k34=(h_i2 -(I2(i)+dt*k33) +cpy_i2*(1./(1+sigma^(-(PY(i)+dt*k13))))-ci1_i2*(1./(1+sigma^(-(I1(i)+dt*k23))))+ctc_i2*(1./(1+sigma^(-(TC(i)+dt*k43)))))*tau3; 
        k44=(h_tc -(TC(i)+dt*k43) +cpy_tc*(1./(1+sigma^(-(PY(i)+dt*k13))))-cre_tc*(1./(1+sigma^(-(RE(i)+dt*k53)))))*tau4+miangin_tc;
        k54=(h_re -(RE(i)+dt*k53) +cpy_re*(1./(1+sigma^(-(PY(i)+dt*k13)))) + ctc_re*(1./(1+sigma^(-(TC(i)+dt*k43)))) -cre_re*(1./(1+sigma^(-(RE(i)+dt*k53)))))*tau5;     
        k64=(h_ei-(EI(i)+dt*k63)+cpy_ei*(1./(1+sigma.^(-(PY(i)+dt*k13))))+ctc_ei*(1./(1+sigma.^(-(TC(i)+dt*k43))))-ci1_ei*(1./(1+sigma^(-(I1(i)+dt*k23)))))*tau6;
        
        PY(i+1)=PY(i)+dt*(k11+2*k12+2*k13+k14)/6;
        I1(i+1)=I1(i)+dt*(k21+2*k22+2*k23+k24)/6;
        I2(i+1)=I2(i)+dt*(k31+2*k32+2*k33+k34)/6;
        TC(i+1)=TC(i)+dt*(k41+2*k42+2*k43+k44)/6;
        RE(i+1)=RE(i)+dt*(k51+2*k52+2*k53+k54)/6;
        EI(i+1)=EI(i)+dt*(k61+2*k62+2*k63+k64)/6;
 
   
    ux(i)=(PY(i)+I1(i)+I2(i)+EI(i))./4;
     
end 

%The following codes are used for bifurcation analysis
TD=78; % starting time used for analysis 
uxcut=ux(TD/dt+1:end);
peaksmax=findpeaks(uxcut);
peaksmax=sort(peaksmax);
peaksmin=findpeaks(-uxcut);
peaksmin=sort(-peaksmin);

if (length(peaksmax)>=1)
    pmax1=peaksmax(end);
    pmax2=peaksmax(1);
    index=find(pmax1-peaksmax>=0.001);
    if (length(index)>=1)
        pmax2=peaksmax(index(end));
    end
else
    pmax1=uxcut(1);
    pmax2=uxcut(2);
end


if (length(peaksmin)>=1)
    pmin1=peaksmin(end);
    pmin2=peaksmin(1);
else
    pmin1=uxcut(1);
    pmin2=uxcut(2);
end  
   
%The following codes are used for frequency analysis
if (length(peaksmax)==0||length(peaksmin)==0||abs(pmax1-pmin2)<0.01)
   fd=0;
else 
Fs=1/dt;
N_n=length(ux(50/dt+1:time/dt-1));
mean_uu=ux(50/dt+1:time/dt-1)-mean(ux(50/dt+1:time/dt-1));
X1=fft(mean_uu);
f1=(0:length(X1)/2-1)*Fs./(length(X1));
Power=X1.*conj(X1)./N_n;
index=find(Power(1:N_n/2)==max(Power(1:N_n/2)));
fd=f1(index); %computing the donamite frequency
end 
%The following codes are used for state analysis
    state=1;  % low firing 
    if (abs(pmin1-pmin2)<0.2 && abs(pmax1-pmin1)<0.12 && abs(pmax1-pmin1)>=0.01 && fd<=3.5)
    state=2;  %preictal
    elseif (abs(pmin1-pmin2)>=0.004 && fd>=2 && fd<=4)%0.01
    state=4;   % SWD
    elseif (abs(pmin1-pmin2)>=0.01 && fd>7 ) %% 
     state=5;   % atypical swd   
    elseif (abs(pmin1-pmin2)<0.15 && abs(pmax1-pmin2)>=0.01 && fd<=7)
    state=6;   % clonic dis
    elseif (abs(pmin1-pmin2)<0.01 && abs(pmax1-pmin2)>=0.01 && fd>7)
    state=7;   % tonic dis
    elseif (  pmax1<-0.1 && pmax1>=-0.8 && fd<0.1)
    state=3;   % slow rythmic activity
    end
 end

