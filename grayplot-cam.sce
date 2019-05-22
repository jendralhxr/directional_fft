// run it this way:
// $ scilab-adv-cli -f grayplot-cam.sce -args cameradisplacementfile

clear;
exec('peak_detect.sci');
funcprot(0); // peaks() is redifined
args= sciargs();
temp=csvRead(args(5), ascii(9), 'double');
threshold_ratio= 32;

sample_window= 120; // FFT sample_window window
sample_time= 10; // s
freq_sampling= 240; // Hz
freq_max=30; // Hz
peaks_order= 8; // number of frequency peaks to extract

freq_elem= freq_max / freq_sampling * sample_window; 
freqF=0:freq_sampling/sample_window:freq_sampling-0.0001;
freqS=0:freq_max/freq_elem:freq_max-0.0001;

// parsing lateral and vertical displacement
for marker= 1:9 do
  displacement_x(:,marker)=temp(:,2*marker-1);
  displacement_y(:,marker)=temp(:,2*marker);
end

step_size=30;
step_max=(size(displacement_y,1)-sample_window)/step_size; 
time=0:sample_time/step_max:sample_time;

// plot frequency responses
for marker= 1:9 do
  for step=0:step_max do
    dfft_x=abs(fft(displacement_x(step*step_size+1:step*step_size+sample_window,marker)));
    dfft_y=abs(fft(displacement_y(step*step_size+1:step*step_size+sample_window,marker)));
    dfft_x(1)=0;
    dfft_y(1)=0;

    // extract largest lateral peak
    [ fft_max_value_x(step+1) fft_max_index_x(step+1) ] = max(dfft_x);    
	if freqF(fft_max_index_x(step+1)) < freq_max then
		fft_max_freq_x(marker,step+1)= freqF(fft_max_index_x(step+1));
        else fft_max_freq_x(marker,step+1)= 0;
        end
    // extract largest vertical peak
    [ fft_max_value_x(step+1) fft_max_index_y(step+1) ] = max(dfft_x);    
    if freqF(fft_max_index_y(step+1)) < freq_max then
		fft_max_freq_y(marker,step+1)= freqF(fft_max_index_y(step+1));
        else fft_max_freq_y(marker,step+1)= 0;
        end
	// extract frequency peaks
	temp_peaks_x= peak_detect(dfft_x', max(dfft_x)/threshold_ratio);
	temp_peaks_y= peak_detect(dfft_y', max(dfft_y)/threshold_ratio);
	for order=1:size(temp_peaks_x,2) do
        fft_peaks_freq_x(step+1,order)= freqF(temp_peaks_x(order));
	end		// order
	for order=1:size(temp_peaks_y,2) do
        fft_peaks_freq_y(step+1,order)= freqF(temp_peaks_y(order));
	end		// order
	
    grayft_x(step+1,1:freq_elem)=dfft_x(1:freq_elem)';
    grayft_y(step+1,1:freq_elem)=dfft_y(1:freq_elem)';


  end // step
  
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

  filename=sprintf("cam-freqx%d.csv",marker-1);
  csvWrite(fft_peaks_freq_x, filename, ascii(9));
  filename=sprintf("cam-freqy%d.csv",marker-1);
  csvWrite(fft_peaks_freq_y, filename, ascii(9));
  
end // marker

filename=sprintf("cam-peakx.csv");
csvWrite(fft_max_freq_x', filename, ascii(9));
filename=sprintf("cam-peaky.csv");
csvWrite(fft_max_freq_y', filename, ascii(9));

exit
