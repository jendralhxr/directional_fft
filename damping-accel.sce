// run it this way:
// $ scilab-adv-cli -f damping-accel.sce -args accelerometerdisplacementfile
clear;

exec('peak_detect.sci');
funcprot(0); // peaks() is redifined
args= sciargs();
displacement=csvRead(args(5), ascii(9), 'double');

threshold_ratio= 8;
smoothing_ratio= 20;
freq_sampling= 480; // Hz

freqZ=0:freq_sampling/size(displacement,1):freq_sampling-0.0001;

for marker= 1:9 do
    dfft=abs(fft(displacement(:,marker)));
	dfft(1)=0;

	// find peaks from smoothed signal
	sdfft=smooth( [freqZ;dfft'], (freqZ(2)-freqZ(1))*smoothing_ratio  ); 
	temp_peaks= peak_detect(sdfft(2,:), max(dfft)/threshold_ratio);
	
	for order=1:size(temp_peaks,2)/2 do
		mid= (temp_peaks(order)-1)*smoothing_ratio;
		
		// lower band limit
		i= mid-1;
		while i>1,
			if dfft(i) > dfft(mid) then mid=i; end;
			if dfft(i) < dfft(mid)/sqrt(2) then break; end;
			i= i-1;
		end
		lo = i;
			
		// upper band limit
		i= mid+1;
		while i<size(dfft,1),
			if dfft(i) > dfft(mid) then mid=i; end;
			if dfft(i) < dfft(mid)/sqrt(2) then break; end;
			i= i+1;
		end
		hi = i;
		
		freqN(marker,order)= freqZ(mid); 
		zeta(marker,order)= ( freqZ(hi) - freqZ(lo)) / (2 * freqZ(mid));
		
	end		
end

filename=sprintf("freq.csv");
csvWrite(freqN, filename, ascii(9));
 
filename=sprintf("zeta.csv");
csvWrite(zeta, filename, ascii(9));
 
    
exit
