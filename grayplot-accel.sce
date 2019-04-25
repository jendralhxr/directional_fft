// run it this way:
// $ scilab-cli -f grayplot-cam.sce -args accelerometerdisplacementfile

clear;
args= sciargs();
displacement=csvRead(args(5), ascii(9), 'double');

sample= 1200; // FFT sample window
sample_time= 20; // s
freq_sampling= 480; // Hz
step_size=60;
step_max= (size(displacement,1)-sample)/step_size; 

freq_max=30; // Hz
freq_elem= freq_max / freq_sampling * sample; 

freqS=0:freq_max/freq_elem:freq_max-0.0001;
time=0:sample_time/step_max:sample_time;

// plot frequency responses
for marker= 1:9 do

  for step=0:step_max do
    dfft=abs(fft(displacement(step_size*step+1:step_size*step+sample,marker)));
    dfft(1)=0;
    grayft(step+1,1:freq_elem)=dfft(1:freq_elem)';
  end

// vertical axis
  xset("colormap",jetcolormap(64));
  Sgrayplot(time,freqS,grayft);
  xlabel("time (s)");
  ylabel("freq (Hz)");
  label=sprintf("marker%d vertical",marker-1);
  title(label);
  filename=sprintf("accel-y%d.png",marker-1);
  xs2png(0,filename);
 
end

exit
