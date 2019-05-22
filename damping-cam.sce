// run it this way:
// $ scilab-adv-cli -f damping-cam.sce -args accelerometerdisplacementfile
clear;

exec('peak_detect.sci');
funcprot(0); // peaks() is redifined
args= sciargs();
displacement=csvRead(args(5), ascii(9), 'double');

threshold_ratio= 8;
freq_sampling= 240; // Hz

freqZ=0:freq_sampling/size(displacement,1):freq_sampling-0.0001;

for marker= 1:9 do
	displacement_x(:,marker)=displacement(:,2*marker-1);
	displacement_y(:,marker)=displacement(:,2*marker);

    dfft_x=abs(fft(displacement_x(:,marker)));
	dfft_y=abs(fft(displacement_y(:,marker)));
	dfft_x(1)=0;
	dfft_y(1)=0;
	temp_peaks_x= peak_detect(dfft_x', max(dfft_x)/threshold_ratio);
	temp_peaks_y= peak_detect(dfft_y', max(dfft_y)/threshold_ratio);
	
	// lateral	
	for order=1:size(temp_peaks_x,2) do
		mid = temp_peaks_x(order);
		// lower band limit
		i= mid-1;
		while i>1,
			if dfft_x(i) < dfft_x(mid)/sqrt(2) then break; end;
			i= i-1;
		end
		lo = i;
		// upper band limit
		i= mid+1;
		while i<size(dfft_x,1),
			if dfft_x(i) < dfft_x(mid)/sqrt(2) then break; end;
			i= i+1;
		end
		hi = i;
		freqN_x(marker,order)= freqZ(mid); 
		zeta_x(marker,order)= ( freqZ(hi) - freqZ(lo)) / (2 * freqZ(mid));
	end
	
	// vertical
	for order=1:size(temp_peaks_y,2) do
		mid = temp_peaks_y(order);
		// lower band limit
		i= mid-1;
		while i>1,
			if dfft_y(i) < dfft_y(mid)/sqrt(2) then break; end;
			i= i-1;
		end
		lo = i;
		// upper band limit
		i= mid+1;
		while i<size(dfft_y,1),
			if dfft_y(i) < dfft_y(mid)/sqrt(2) then break; end;
			i= i+1;
		end
		hi = i;
		freqN_y(marker,order)= freqZ(mid); 
		zeta_y(marker,order)= ( freqZ(hi) - freqZ(lo)) / (2 * freqZ(mid));
	end
			
end

filename=sprintf("freq-x.csv");
csvWrite(freqN_x, filename, ascii(9));
filename=sprintf("freq-y.csv");
csvWrite(freqN_y, filename, ascii(9));
 
filename=sprintf("zeta-x.csv");
csvWrite(zeta_x, filename, ascii(9));
filename=sprintf("zeta-y.csv");
csvWrite(zeta_y, filename, ascii(9));
 
    
exit
