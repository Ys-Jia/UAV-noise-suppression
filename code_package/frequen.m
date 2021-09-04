function [realf, realy] = frequen(x,Fs)
  y = fft(x); 
  n = length(x);                         
  realy=2*abs(y(1:n/2+1))/n;
  realf=(0:n/2)*(Fs/n);
end
