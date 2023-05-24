clc
clear all
close all

%% b
n = 1:100;
s = cos(0.04*pi*n);
y = randn(1,100);
y = (y-mean(y))./(std(y));
MEAN_y = mean(y)
VAR_y = var(y)

%% c
b_v = [1];
a_v = [1 -0.7];

b_vp = [1];
a_vp = [1 0.5 -0.36];

v = filter(b_v,a_v,y);
vprime = filter(b_vp,a_vp,y);

z = s + v;

figure
plot(s,'r','linewidth',2)
hold on
plot(z,'b','linewidth',2)
hold on
grid on
xlabel('Samples')
ylabel('Amplitude')
title('Signal s(n) & Observations z(n)')
legend('s(n)','z(n)','location','best')

figure
plot(vprime,'linewidth',2)
hold on
grid on
xlabel('Samples')
ylabel('Amplitude')
title('Noise v_{prime}(n)')
legend('v_{prime}(n)','location','best')

%% d
for k = 1:100
    
    SUM1 = 0;
    SUM2 = 0;
    
    for j = k+1:100
        
        SUM1 = (v(j)*vprime(j-k+1))+SUM1;
        SUM2 = (vprime(j)*vprime(j-k+1))+SUM2;
        
    end
    
    Rvvprime(k) = 0.01*SUM1;
    Rvprime(k) = 0.01*SUM2;
    
end

figure
stem(Rvvprime)
hold on
autocorr(z,99)
hold on
grid on
ylabel('R_{vv_{prime}}(k) = R_z(k)')
legend('Estimated Autocorrelation','Real Autocorrelation')
title('Real & Estimated Autocorrelation for R_{vv_{prime}}')

figure
stem(Rvprime)
hold on
autocorr(vprime,99)
hold on
grid on
ylabel('R_{v_{prime}}(k)')
legend('Estimated Autocorrelation','Real Autocorrelation')
title('Real & Estimated Autocorrelation for R_{v_{prime}}')

%% e
N = 4;

for i = 1:N
RVP(1,i) = Rvprime(i);
rvvp(i,1) = Rvvprime(i);
end

% Rvp = toeplitz([Rvprime(1) Rvprime(2) Rvprime(3) Rvprime(4)])
% rvvp = [Rvvprime(1) Rvvprime(2) Rvvprime(3) Rvvprime(4)]'
Rvp = toeplitz(RVP);
h = inv(Rvp)*rvvp

%% f
% vhat = filter([h(1) h(2) h(3) h(4)],[1],vprime);
vhat = filter(h,1,vprime);

shat = z - vhat;

figure
plot(s,'r','linewidth',2)
hold on
plot(shat,'b','linewidth',2)
hold on
grid on
xlabel('Samples')
ylabel('Amplitude')
ylim([-5 4])
title(['Real Signal s(n) & Estimated Signal s_{hat}(n) for N = ',num2str(N)])
legend('s(n)','s_{hat}(n)','location','best')

figure
plot(s,'r','linewidth',2)
hold on
plot(shat,'b','linewidth',2)
hold on
plot(z,'g','linewidth',2)
hold on
grid on
xlabel('Samples')
ylabel('Amplitude')
ylim([-5 4])
title(['s(n) , s_{hat}(n) , z(n) for N = ',num2str(N)])
legend('s(n)','s_{hat}(n)','z(n)','location','best')

%% g
Eta_hat = (1/100)*sum(v);
varv = (1/100)*sum((v-Eta_hat).^2)

MSE = varv-h'*rvvp
experimental_MSE = (v-vhat)*(v-vhat)'/100

%% h
% Change N in line 85