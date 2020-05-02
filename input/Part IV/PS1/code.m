%% Data import
load 'data.mat';
x = data(:,1); %(nX1) vector (age)
y = data(:,6); %(nX1) vector (happiness)
w = data(:,7); %(nX1) vector (income)

w_max = quantile(w,0.95); %111.5
w_min = quantile(w,0.05); %4.9
dummy = ((w<=w_max) & (w>=w_min)); 

y = y(dummy); %we are excluding some outliers 
x = x(dummy); %we are excluding some outliers 
w = w(dummy); %we are excluding some outliers 
w = log(w); %we will use log(w) instead of just w

%% 1.[a] histogram for x and y
%% 
histogram(x) % counts number of observations falling into each bin
ylabel('frequency')
xlabel('x')
title('Histogram(x)');
histogram(x,'Normalization','probability') % relative frequency
%% 


histogram(y) 
ylabel('frequency')
xlabel('y')
title('Histogram(y)');
histogram(y,'Normalization','probability') 

%% 1. [b] naive density estimator for x and y
N = 100; %number of estimation points in the domain of x and y
x_0 = linspace(min(x), max(x), N);
y_0 = linspace(min(y), max(y), N);

% naive density estimation(=Uniform kernel density estimation)
uniform=@(u)1*(abs(u.^1)<=1/2);
f_x_uniform = zeros(1,N); %(NX1)vector
f_y_uniform = zeros(1,N); %(NX1)vector
h_x = std(x)*(length(x))^(-1/5); %2.7028
h_y = std(y)*(length(y))^(-1/5); %0.2031
for k = 1:N
    z_x = uniform((x - x_0(k))/h_x);
    f_x_uniform(k)= (1/(length(x)*h_x))*sum(z_x);
    z_y = uniform((y - y_0(k))/h_y);
    f_y_uniform(k)= (1/(length(y)*h_y))*sum(z_y);
end

plot(x_0,f_x_uniform,'r')
ylabel('f(x)')
xlabel('x')
title('Naive density estimation for x');

plot(y_0,f_y_uniform,'r')
ylabel('f(y)')
xlabel('y')
title('Naive density estimation for y');

%% 1. [c] kernel density estimator for x and y
gauss=@(u)exp(-u.*u/2)/sqrt(2*pi);
epanechnikov=@(u)(3/4)*(1-u.*u).*(abs(u.^1)<= 1);
parzen=@(u)(1-6*u.*u+6*abs(u.^1).^3).*(abs(u.^1)<= 1/2)+2*((1-abs(u.^1)).^3).*(abs(u.^1)<= 1).*(abs(u.^1)> 1/2);
bartlett=@(u)(1-abs(u.^1)).*(abs(u.^1)<=1);

% work with above kernel functions to get the density estimators 
