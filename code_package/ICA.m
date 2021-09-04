function [Xe,average,eigvalue,B,Q] = ICA(X)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[M,T] = size(X); 
average= mean(X,2);                                 
for i=1:M
    X(i,:)=X(i,:)-average(i)*ones(1,T);
end                                                  
Cx = cov(X',1); 
[eigvector,eigvalue] = eig(Cx); 

Q=eigvector*eigvalue^(-1/2)*eigvector'; 
Z=Q*X; 

 Maxcount=100; 
Critical=0.000001; 
m=M; 
B=rand(m);
for n=1:m
    WP=B(:,n);
    % Y=WP'*Z;
    % G=Y.^3;
    % GG=3*Y.^2;
    count=0;
    LastWP=zeros(m,1);
    B(:,n)=B(:,n)/norm(B(:,n));
    while abs(WP-LastWP) & abs(WP+LastWP)>Critical
        count=count+1
        LastWP=WP
        % WP=1/T*Z*((LastWP'*Z).^3)'-3*LastWP;
        for i=1:m
            WP(i)=mean(Z(i,:).*(tanh((LastWP)'*Z)))-(mean(1-(tanh((LastWP))'*Z).^2)).*LastWP(i);
        end
        WPP=zeros(m,1);
        for j=1:n-1
            WPP=WPP+(WP'*B(:,j))*B(:,j);
        end
        WP=WP-WPP;
        WP=WP/(norm(WP))
        if count==Maxcount
            break
%             fprintf('cannot find signals');
%             return;
        end

    end
    B(:,n)=WP;
end
C=B'
Xe=B'*Z;
end

