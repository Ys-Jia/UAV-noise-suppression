function [realf, realy] = frequen(x,Fs)
y = fft(x); 
%��������ת������ʾΪƵ��f= n*(fs/N)
% f = (0:length(y)-1)*Fs/length(y);
% length(y)
n = length(x);                         
realy=2*abs(y(1:n/2+1))/n;
realf=(0:n/2)*(Fs/n);

% figure;
% plot(realf,realy);

end
