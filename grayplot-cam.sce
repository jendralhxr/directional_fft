temp=csvRead('cam-tap', ascii(9), 'double');

sample= 600; // FFT sample window
sample_time= 10; // s
freq_sampling= 240; // Hz
freq_max=30; // Hz
freq_elem= freq_max / freq_sampling * sample; 
freqS=0:freq_max/freq_elem:freq_max-0.0001;
step_max= 15; 
time=0:sample_time/step_max:sample_time;

// parsing lateral and vertical displacement
for marker= 1:9 do
  displacement_x(:,marker)=temp(:,2*marker-1);
  displacement_y(:,marker)=temp(:,2*marker);
end


// plot frequency responses
for marker= 1:9 do
  for step=0:step_max do
    dfft_x=abs(fft(displacement_x(step*120+1:step*120+sample,marker)));
    dfft_y=abs(fft(displacement_y(step*120+1:step*120+sample,marker)));
    dfft_x(1)=0;
    dfft_y(1)=0;
    grayft_x(step+1,1:freq_elem)=dfft_x(1:freq_elem)';
    grayft_y(step+1,1:freq_elem)=dfft_y(1:freq_elem)';
  end
  
// lateral axis
  xset("colormap",jetcolormap(64));
  Sgrayplot(time,freqS,grayft_x);
  xlabel("time (s)");
  ylabel("freq (Hz)");
  label=sprintf("marker%d lateral",marker-1);
  title(label);
  filename=sprintf("cam-x%d.png",marker-1);
  xs2png(0,filename);

// vertical axis
  xset("colormap",jetcolormap(64));
  Sgrayplot(time,freqS,grayft_y);
  xlabel("time (s)");
  ylabel("freq (Hz)");
  label=sprintf("marker%d vertical",marker-1);
  title(label);
  filename=sprintf("cam-y%d.png",marker-1);
  xs2png(0,filename);

end



