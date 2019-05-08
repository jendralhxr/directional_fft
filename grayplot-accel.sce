// run it this way:
// $ scilab-adv-cli -f grayplot-cam.sce -args accelerometerdisplacementfile

clear;
exec('peak_detect.sci');
funcprot(0); // peaks() is redifined
args= sciargs();
displacement=csvRead("tapa", ascii(9), 'double');

sample_window= 1200; // FFT sample_window window
sample_time= 20; // s
freq_sampling= 480; // Hz
step_size=60;
step_max= (size(displacement,1)-sample_window)/step_size; 
freq_max=30; // Hz
peaks_order= 10; // number of frequency peaks to extract

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
        else fft_max_freq(step+1)= 0;
	end
	
	// extract peaks
	temp_peaks= peak_detect(dfft');
	for order=1:peaks_order do
        fft_peaks_freq(step+1,order)= freqF(temp_peaks(order));
	end		
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
  
  filename=sprintf("accel-freq%d.csv",marker-1);
  csvWrite(fft_peaks_freq, filename, ascii(9));
 
end

filename=sprintf("accel-peak.csv");
csvWrite(fft_max_freq', filename, ascii(9));

exit
