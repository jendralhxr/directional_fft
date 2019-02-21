for marker= 1:9 do
	frlat(:,marker)=real(fft(lat(:,marker)))/1200;
	filat(:,marker)=imag(fft(lat(:,marker)))/1200;
	frver(:,marker)=real(fft(vert(:,marker)))/1200;
	fiver(:,marker)=imag(fft(vert(:,marker)))/1200;
end
