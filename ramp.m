  function input_ramp=ramp(p0,p_end)
dt=1/256; % Intergation step
% Simulation time

time1=10;
time2=10;
time3=20;
time4=10;
time5=10;
time=time1+time2+time3+time4+time5;
time_length=(time)/dt;
ds=linspace(0,(time),time_length);
m=(p_end-p0)/time2;
f_ramp=@(t) p0+m.*(t-time1).*(t>time1)-m.*(t-(time1+time2)).*(t>time1+time2)-m.*(t-(time1+time2+time3)).*(t>time1+time2+time3)+m.*(t-(time1+time2+time3+time4)).*(t>time1+time2+time3+time4);
input_ramp= f_ramp(ds);

  end
% plot(ds,miangin)