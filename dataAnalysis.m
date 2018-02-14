clear all
close all
clc
%%%%%%%%%% Turbulence Project Wind Tunnel Data Analysis, Ka Ling WU %%%%%%%%%%
%%
Ut = load('veldata.txt');
N = length(Ut);
t = ([1:1:N]-1)/(20*1000);
%%%%% 1.1 Taylor's Hypothesis: x = Umean*t %%%%%
%%
Umean = mean(Ut); %mean velocity
x = Umean.*t;

figure('Position',[0 0 1500 280])
plot(x,Ut)
xlabel('x [m]')
ylabel('U [m/s]')
xlim([0 6700])
set(gca,'FontSize',20)

v_prime = Ut-Umean;
%Check for Taylor's Hypothesis validity
TI = sqrt(mean(v_prime.*v_prime))/Umean;

%%%%% 1.2 Correlation %%%%%
%%
[ACF,lags]=autocorr(Ut,N-1);

bound = 1/exp(1);
for i = 1:N-1
    if ACF(i) < bound
        lc = x(i-1); %Correlation length
        ii = i-1;
        break
    end
end

%%

txt1 = ['L_c=0.3568 m'];
txt2 = ['1/e'];
a(1:2)=[lc lc]; b(1:2)=[-0.2 bound];
c(1:2)=[0 lc]; d(1:2)=[bound bound];

figure('Position',[0 0 1000 500])
plot(x,ACF,'LineWidth',3)
hold on
plot(a,b,'-k','LineWidth',2)
hold on
plot(c,d,'-k','LineWidth',2)
hold on
text(0.45,-0.08, txt1,'FontSize',20);
hold on
text(0.05,0.3, txt2,'FontSize',20);
hold on
xlabel('l [m]')
ylabel('C(l)')
xlim([0 5])
set(gca,'FontSize',20)

%%%%% 1.3 Energy Spectrum of the Flow %%%%%
%%
L = x(N);
N = length(x);
Fs = 1/20000;
w = (0:N-1).*2*pi/L;

xdft = fft(v_prime)/N;

for i=1:fix(N/2)-1
     E_k(i)=1/2*(abs(1/sqrt(2*pi/L)*xdft(i+1)))^2;
end

% Since E(k) = E_tilde(k)+E_tilde(-k).
E = E_k*2;
E(N/2)=(1/2)*(abs(1*xdft(i+1)))^2;

% Check normalization, E_integrate should equal to E_parseval
E_integrate = trapz(w(1:fix(N/2)),E);
E_parseval = 0.5*mean(v_prime.*v_prime);

%Denoie data using moving averaging on the spectrum
E_smooth = smooth(E,40);

a1(1:2)=[17.61 17.61]; b1(1:2)=[10^(-20) 10^(5)];
a2(1:2)=[356 356];
a3(1:2)=[3490 3490];
index = find(E~=0);
E_index = E(index);

%%
figure
loglog(w(1:fix(N/2)),E)
hold on
loglog(w(1:fix(N/2)),10*w(1:fix(N/2)).^(-5/3),'LineWidth',2)
hold on
plot(a1,b1,'k','LineWidth',2)
hold on
plot(a2,b1,'--k','LineWidth',2)
hold on
plot(a3,b1,'--k','LineWidth',2)
hold on
ylim([10^(-18) 10^(3)])
xlim([10^(-3.5) 10^(4)])
xlabel('Wavenumber k [rad/m]')
ylabel('E(k)')
legend('Energy Spectrum','Kolmogorov Slope (-5/3)','k_c=17.61','k_{10\eta}=356')
set(gca,'FontSize',20)

%%
figure
loglog(w(1:fix(N/2)),E_smooth)
hold on
loglog(w(1:fix(N/2)),10^(-0.02)*w(1:fix(N/2)).^(-5/3),'LineWidth',2)
hold on
plot(a1,b1,'k','LineWidth',2)
hold on
plot(a2,b1,'--k','LineWidth',2)
hold on
ylim([10^(-18) 10^(3)])
xlim([10^(-3.5) 10^(4)])
xlabel('Wavenumber k [rad/m]')
ylabel('E(k)')
legend('Energy Spectrum','Kolmogorov Slope (-5/3)','k_c=17.61','k_{10\eta}=356')
set(gca,'FontSize',20)


%% Estimate Reynolds Number
Re_flow = Umean*lc/1.568e-5;
Re_integralScale = sqrt(mean(v_prime.*v_prime))*lc/1.568e-5;
Re_taylor = sqrt(Re_integralScale);
%%%%% 1.4 Velocity Increments %%%%%
%%
l = [0.5*10^-3 0.01 0.1 10]; %different length scales
index_v = round(l./(x(N)-x(N-1))); %indexing

for j = 1:4
    for i =1:N-index_v(j)
        delta_v(j,i)=v_prime(i+index_v(j))-v_prime(i);
    end
end

str = ['a','b','c','d'];
figure('Position',[0 0 1500 1000])
for j = 1:4
    subplot(4,1,j)
    plot(x(1:N-index_v(j)),delta_v(j,1:N-index_v(j)))
    title(['(' str(j) ') l=' num2str(l(j)) 'm'])
    xlim([0 6700])
    ylabel('\delta v_{||}(x)')
    xlabel('x [m]')
    set(gca,'FontSize',20)
end
%%%%% 1.5 Statistics of Velocity Increments %%%%%
%%
% Histogram for PDF
for j = 1:4
    [prob,value]=hist(delta_v(j,1:N-index_v(j)),1000);
    PDF_delta_v(j,:)=prob(:);
    value_v(j,:)=value(:);
end
% Generating Gaussian Distribution for different length scales
for j = 1:4
[mu,s,muci,sci] = normfit(delta_v(j,1:N-index_v(j)));
k = [-10:0.001:10];
norm = normpdf(k,mu,s);
gaussian(j,:) = norm(1,:);
end
%%
figure('Position',[0 0 1500 450])
for j = 1
    subplot(1,4,j)
    bar(value_v(j,:),PDF_delta_v(j,:)./trapz(value_v(j,:),PDF_delta_v(j,:)))
    hold on
    plot(k,gaussian(j,:),'r','LineWidth',2)
    title(['(' str(j) ') l=' num2str(l(j)) 'm'])
    ylabel('PDF(\delta v_{||})')
    xlabel('\delta v_{||} [m/s]')
    set(gca,'FontSize',20)
    xlim([-0.2 0.2])
end
for j = 2
    subplot(1,4,j)
    bar(value_v(j,:),PDF_delta_v(j,:)./trapz(value_v(j,:),PDF_delta_v(j,:)))
    hold on
    plot(k,gaussian(j,:),'r','LineWidth',2)
    title(['(' str(j) ') l=' num2str(l(j)) 'm'])
    ylabel('PDF(\delta v_{||})')
    xlabel('\delta v_{||} [m/s]')
    set(gca,'FontSize',20)
    xlim([-2 2])
end

for j = 3
    subplot(1,4,j)
    bar(value_v(j,:),PDF_delta_v(j,:)./trapz(value_v(j,:),PDF_delta_v(j,:)))
    hold on
    plot(k,gaussian(j,:),'r','LineWidth',2)
    title(['(' str(j) ') l=' num2str(l(j)) 'm'])
    ylabel('PDF(\delta v_{||})')
    xlabel('\delta v_{||} [m/s]')
    set(gca,'FontSize',20)
    xlim([-5 5])
end

for j = 4
    subplot(1,4,j)
    bar(value_v(j,:),PDF_delta_v(j,:)./trapz(value_v(j,:),PDF_delta_v(j,:)))
    hold on
    plot(k,gaussian(j,:),'r','LineWidth',2)
    title(['(' str(j) ') l=' num2str(l(j)) 'm'])
    ylabel('PDF(\delta v_{||})')
    xlabel('\delta v_{||} [m/s]')
    set(gca,'FontSize',20)
    xlim([-10 10])
end

%%%%% 1.6 Flatness of Velocity Increments %%%%%
%%
dx = L/N;
l = logspace(0,6,75).*dx;
index_v = round(l./(x(N)-x(N-1)));
%%
delta_v=zeros(length(index_v),N);
delta_v_4thM=zeros(length(index_v),N);
delta_v_2ndMSquare=zeros(length(index_v),N);
flatness=zeros(length(index_v),1);

for j = 1:1:length(index_v)
    for i =1:N-index_v(j)
        delta_v(j,i)=v_prime(i+index_v(j))-v_prime(i);
    end
end

for j = 1:1:length(index_v)
     delta_v_4thM(j) = mean(delta_v(j,1:N-index_v(j)).^4,2);
     delta_v_2ndMSquare(j) = mean(delta_v(j,1:N-index_v(j)).^2,2)^2;
     flatness(j)=delta_v_4thM(j)./delta_v_2ndMSquare(j);
end

%% 
a1(1:2)=[0.363 0.363]; b1(1:2)=[0 10];
a2(1:2)=[0.026 0.026];
a3(1:2)=[10^(-4) 10^3]; b3(1:2)=[3 3];
txt1 = ['L_c=0.3568 m'];
txt2 = ['10\eta=0.018 m'];

figure('Position',[0 0 1500 450])
semilogx(l,flatness,'LineWidth',2)
hold on
plot(a1,b1,'k','LineWidth',2)
hold on
plot(a2,b1,'--k','LineWidth',2)
hold on
plot(a3,b3,'-k','LineWidth',2)
hold on
text(0.45,7, txt1,'FontSize',20);
hold on
text(0.03,7, txt2,'FontSize',20);
xlabel('l [m]')
ylabel('f(l)')
ylim([2 8])
set(gca,'FontSize',20)      
ylabel('f(l)')
xlabel('l [m]')
