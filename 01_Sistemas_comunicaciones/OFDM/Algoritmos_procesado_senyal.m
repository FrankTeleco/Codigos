clc; clear; close all;

%y[n]=0.6y[n−1]+x[n]

x= [ 2 3 1 3 5 2];

L= length(x);

y= zeros(L,1);
y(1)=x(1);
for k = 2: L

    y(k)=0.6*y(k-1)  + x(k);

end


clc;clear; close all;
%y[n]=0.5y[n−1]+0.25x[n]+0.25x[n−1]

x=[ 4 2 6 1 3 5];

L= length(x);

y= zeros(L,1);

y(1)=0.25*x(1); 

for k = 2: L

y(k)= 0.5*y(k-1)+0.25*x(k)+0.25*x(k-1);
end

clc
