clear;
funcprot(0);
exec('peak_detect.sci');

ANGLES= 72;
SAMPLE_POINTS= 9;
SAMPLING_FREQ= 240; // Hz

x_orig=csvRead("cn2l.txt", ascii(9), 'double');
y_orig=csvRead("cn2v.txt", ascii(9), 'double');
size_m = size(x_orig);
freq=0:SAMPLING_FREQ/(size_m(1)-1):SAMPLING_FREQ;

teta=0;
for step= 0:ANGLES do
    teta= step * %pi/ANGLES;
    mprintf("teta %f: ", teta);
    for point= 1:SAMPLE_POINTS do
        mprintf("%d ", point);
        //check about positive/negative sign
		//x_shifted(:,point) = x_orig(:,point)*cos(teta) + y_orig(:,point)*sin(teta);
        temp= x_orig(:,point)*cos(teta) + y_orig(:,point)*sin(teta);
        temp_f= abs(fft(temp))';
        temp_f(1)= 0;
        ts= max(temp_f)/20;
        temp_peaks= peak_detect(temp_f,ts);
        // plot2d(freq,temp_f,2);
        // plot2d(freq(peaks),temp_f(peaks),-3);
        
        x_freq((step+1),point,:) = temp_f;
        //peak_freq((step+1),point,:) = temp_peaks;
		// may not need the y axis
        //y_shifted(:,point)= y_orig(:,point)*cos(teta) - x_orig(:,point)*sin(teta);			
        end
    mprintf("\n");
end
