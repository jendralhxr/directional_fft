// run it this way:
// $ scilab-adv-cli -f grayplot-accel.sce -args accelerometerdisplacementfile

clear;
exec('peak_detect.sci');
funcprot(0); // peaks() is redefined
args= sciargs();
displacement=csvRead(args(5), ascii(9), 'double');

threshold_ratio= 32;
sample_window= 240; // FFT sample_window window
sample_time= 20; // s
time_max= 20// sec
time_part= time_max / sample_time; 

freq_sampling= 480; // Hz
step_size=60;
step_max= (size(displacement,1)-sample_window)/step_size; 
freq_max=30; // Hz

freq_elem= freq_max / freq_sampling * sample_window; 
freqF=0:freq_sampling/sample_window:freq_sampling-0.0001;
freqS=0:freq_max/freq_elem:freq_max-0.0001;
time=0:sample_time/step_max:sample_time;
    
// plot frequency responses
for marker= 1:9 do

  for step=0:step_max do
    dfft=abs(fft(displacement(step_size*step+1:step_size*step+sample_window,marker)));
    dfft(1)=0;
    
    // extract largest peak
    [ fft_max_value(step+1) fft_max_index(step+1) ] = max(dfft);    
	if freqF(fft_max_index(step+1)) < freq_max then
		fft_max_freq(marker,step+1)= freqF(fft_max_index(step+1));
        else fft_max_freq(marker,step+1)= 0;
	end
	
	// extract peaks
	temp_peaks= peak_detect(dfft', max(dfft)/threshold_ratio);
	for order=1:size(temp_peaks,2) do
        fft_peaks_freq(step+1,order)= freqF(temp_peaks(order));
	end		
    grayft(step+1,1:freq_elem)=dfft(1:freq_elem)';
  end

// vertical axis
  xset("colormap",jetcolormap(64));
  Sgrayplot(time(1:step_max*time_part),freqS,grayft(1:step_max*time_part,:));
  //Sgrayplot(time,freqS,grayft);
  xlabel("time (s)");
  ylabel("freq (Hz)");
  label=sprintf("marker%d vertical",marker-1);
  title(label);
  filename=sprintf("accel-y%d.png",marker-1);
  xs2png(0,filename);
  
  filename=sprintf("accel-freq%d.csv",marker-1);
  csvWrite(fft_peaks_freq, filename, ascii(9));
 
end

filename=sprintf("accel-peak.csv");
csvWrite(fft_max_freq', filename, ascii(9));

exit
