clear;
funcprot(0);
exec('peak_detect.sci');

ANGLES= 72;
SAMPLE_POINTS= 9;
SAMPLING_FREQ= 240; // Hz
PEAK_THRESHOLD= 50;
SAMPLING_LENGTH= 600; // samples
//ORDER= 3;

x_orig=csvRead("lat.txt", ascii(9), 'double');
y_orig=csvRead("ver.txt", ascii(9), 'double');
size_m = size(x_orig);
freq=0:SAMPLING_FREQ/(SAMPLING_LENGTH/2-1):SAMPLING_FREQ;

teta=0;
for step= 0:ANGLES do
    teta= step * %pi/ANGLES;
    mprintf("teta %f: ", teta);
    for point= 1:SAMPLE_POINTS do
        mprintf("%d ", point);
        //check about positive/negative sign
		//x_shifted(:,point) = x_orig(:,point)*cos(teta) + y_orig(:,point)*sin(teta);
        temp= x_orig(:,point)*cos(teta) + y_orig(:,point)*sin(teta);
        temp_freq= abs(fft(temp))';
        temp_freq(1)= 0;
        temp_freq= temp_freq(1:SAMPLING_LENGTH/2);
        
        ts= max(temp_freq)/PEAK_THRESHOLD;
        temp_peaks= peak_detect(temp_freq,ts);
        
        mul_freq(point,(step+1),1:SAMPLING_LENGTH/2) = temp_freq;
		//mul_peak((step+1),point,1:ORDER) = temp_peaks(1:ORDER);
	    mul_peak1((step+1),point) = freq(temp_peaks(1));
		mul_peak2((step+1),point) = freq(temp_peaks(2));
		mul_peak3((step+1),point) = freq(temp_peaks(3));
		
		//cylindrical plot routine here
		
		// may not need the y axis
        //y_shifted(:,point)= y_orig(:,point)*cos(teta) - x_orig(:,point)*sin(teta);			
        end
    mprintf("\n");
end
