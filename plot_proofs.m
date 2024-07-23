% Code file for the Figures 4.1 - 4.3 and 5.1 - 5.4

clear, clc, close all
warning('off','MATLAB:integral:NonFiniteValue')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Init switches
save_results = 1; % Switch for saving data for tikz
if (save_results == 1)
    fileID = fopen('plot_proofs.txt','w');
    format = '%d %1.4e \n';
end%if
switch_case1 = 1; % Switch for computations for D1(m)
switch_case2 = 1; % Switch for computations for D2(m)

%% Computations for D1(m) with sinh and cKB

if switch_case1 ==1

%% Init

% Set testing parameters
d = [1/4, 1/2, 3/4];
a = [1, 1.1, 1.2, 1.5, 1.8, 2,3];
m_max = 20;

% Initialize vectors
int_sinh = zeros(length(a),m_max,length(d));
int_cKB = zeros(length(a),m_max,length(d));
term_sinh = zeros(length(a),m_max,length(d));
term_cKB = zeros(length(a),m_max,length(d));

%% Computation of the integrals

for j = 1:length(d)
    delta = d(j)*pi; % Set bandwidth parameter
    for m = 1:m_max % Set truncation parameter
        beta = m.*(pi-delta); % Set shape parameter
            for k = 1:length(a)
                alpha = a(k); % Set scaling factor

                % Computations for the sinh-type window function
                func_sinh = @(t) besseli(1,beta*sqrt(alpha.^2-t.^2))./sqrt(alpha.^2-t.^2);
                int_sinh(k,m,j) = integral(func_sinh,0,1);
                term_sinh(k,m,j) = 1-(alpha*beta)/sinh(alpha*beta)*int_sinh(k,m,j);
                
                % Computations for the continuous Kaiser-Bessel window function
                func_cKB = @(t) fillmissing((t<1).*(sinh(beta*alpha*sqrt(1-t.^2))./sqrt(1-t.^2)-sin(beta*alpha*t)./t),'constant',beta*alpha-sin(beta*alpha)) + (t>1).*(beta*alpha*my_sinc(beta*alpha,sqrt(t.^2-1))-sin(beta*alpha*t)./t);
                int_cKB(k,m,j) = integral(func_cKB,0,1/alpha);
                term_cKB(k,m,j) = 1-2/(pi*(besseli(0,alpha*beta)-1))*int_cKB(k,m,j);
            end%for
    end%for
end%for

%% Visualization of the integrals

% Term D_1(m) for the sinh-type window function
figure(1);
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,abs(term_sinh(1,:,j)),'-o',1:m_max,term_sinh(2,:,j),'-o',1:m_max,term_sinh(3,:,j),'-o',1:m_max,term_sinh(4,:,j),'-o',1:m_max,term_sinh(5,:,j),'-o',1:m_max,term_sinh(6,:,j),'-o',1:m_max,term_sinh(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\alpha=1$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$',num2str(d(j)),'$\pi$']);
    colororder(["#E69F00";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 4.1: Semilogarithmic plots of the term $1 - \frac{\beta}{\sinh \beta}\,\int_0^{1/\alpha} \frac{I_1(\beta \sqrt{1-\nu^2})}{\sqrt{1-\nu^2}}\,\mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.');

% Term D_1(m) for the continuous Kaiser--Bessel window function
figure(2);
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,abs(term_cKB(1,:,j)),'-o',1:m_max,term_cKB(2,:,j),'-o',1:m_max,term_cKB(3,:,j),'-o',1:m_max,term_cKB(4,:,j),'-o',1:m_max,term_cKB(5,:,j),'-o',1:m_max,term_cKB(6,:,j),'-o',1:m_max,term_cKB(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\alpha=1$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$',num2str(d(j)),'$\pi$']);
    colororder(["#E69F00";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 5.1: Semilogarithmic plots of the term $1 - \frac{2\beta}{\pi\,\big(I_0(\beta) - 1\big)}\, \int_0^{1/\alpha} \Big(\frac{\sinh\big(\beta\,\sqrt{1 - \nu^2}\big)}{\beta\,\sqrt{1 - \nu^2}} - \frac{\sin(\beta \nu)}{\beta \nu}\Big)\,\mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.')

%% Generate tables for tikz

if (save_results == 1)
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\mathcal A(m) = 1 - \\frac{\\beta}{\\sinh \\beta} \\, \\int_0^{1/\\alpha}\\frac{I_1(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;term_sinh(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for

fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\Delta_{\\mathrm{cKB},1} = 1 - \\frac{2}{\\pi(I_0(\\beta)-1)} \\, \\int_0^{1/\\alpha}\\frac{\\sinh(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}-\\frac{\\sin(\\beta t)}{t}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;term_cKB(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for
end%if
end%if

%% Computations for D2(m) with sinh and cKB

if switch_case2 ==1

% Set parameters
d = [1/4, 1/2, 3/4];
a = [1, 1.1, 1.2, 1.5, 1.8, 2,3];
m_max = 20;

% Initialize vectors
decay = zeros(length(a),m_max);
int_sinh2 = zeros(length(a),m_max);
int_sinh3 = zeros(length(a),m_max);
int_cKB2 = zeros(length(a),m_max);
int_cKB3 = zeros(length(a),m_max);
    
%% Computation of the integrals

for j = 1:length(d)
    delta = d(j);
    for k = 1:length(a)
        alpha = a(k);
        for m = 1:m_max
            beta = alpha*m*(pi-delta*pi); % Set shape parameter
            v1 = @(w) m/beta*(w+pi); % Init of auxiliary function
            decay(k,m,j) = exp(-beta/alpha); % Exponential decay to compare with
            
            % Computations for the sinh-type window function
            funci_subst = @(t) besseli(1,beta*sqrt(1-t.^2))./sqrt(1-t.^2); % Init of Bessel function
            int_sinh2(k,m,j) = beta/sinh(beta)*integral(funci_subst,1/alpha,1); % Numerical computation of integral
            int_sinh3(k,m,j) = real(beta/sinh(beta)*integral(funci_subst,1/alpha,v1(0))); % Numerical computation of integral
            
            % Computations for the continuous Kaiser--Bessel window function
            func_cKB = @(t) fillmissing((t<1).*(sinh(beta*sqrt(1-t.^2))./sqrt(1-t.^2)-sin(beta*t)./t),'constant',0) + (t==0).*(-beta+sinh(beta)) + (t>=1).*(beta*my_sinc(beta,sqrt(t.^2-1))-beta*my_sinc(beta,t));
            int_cKB2(k,m,j) = 2/(pi*(besseli(0,beta)-1))*integral(func_cKB,1/alpha,1); % Numerical computation of integral
            int_cKB3(k,m,j) = 2/(pi*(besseli(0,beta)-1))*integral(func_cKB,1/alpha,v1(0)); % Numerical computation of integral
            
            % Visualization of the integrand
            if (m==5 && alpha==1 && delta==1/4)
                x = linspace(-1,1);
                y1 = func_cKB(x)/beta;
                figure(3); plot(x,y1); xlabel('$\nu$'); set(gca,'Fontsize',16);
                xticks([-1,0,1]), yticks([0,(sinh(beta)-beta)/beta]); yticklabels({'$0$','$\frac{\sinh(\beta)-\beta}{\beta}$'})
                sgtitle('Figure 5.2: The integrand (5.17) for $m=5$, $\alpha=1$ and $\delta=\frac{\pi}{4}$');
            end%if
       end%for
    end%for
end%for

%% Visualization of the integrals

% Terms for the sinh-type window function
figure(4); 
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,decay(1,:,j),'--k',1:m_max,int_sinh3(1,:,j),'-o',1:m_max,int_sinh3(2,:,j),'-o',1:m_max,int_sinh3(3,:,j),'-o',1:m_max,int_sinh3(4,:,j),'-o',1:m_max,int_sinh3(5,:,j),'-o',1:m_max,int_sinh3(6,:,j),'-o',1:m_max,int_sinh3(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\mathrm e^{-m(\pi-\delta)}$','$\alpha=1$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$ ',num2str(d(j)),'$\pi$']);
    colororder(["#808080";"#E69F00";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 4.2: Semilogarithmic plots of the term $\frac{\beta}{\sinh \beta}\,\int_{1/\alpha}^{v_1(0)} \frac{I_1(\beta \sqrt{1-\nu^2})}{\sqrt{1-\nu^2}}\,\mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1.1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.');
figure(5); 
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,decay(1,:,j),'--k',1:m_max,int_sinh2(2,:,j),'-o',1:m_max,int_sinh2(3,:,j),'-o',1:m_max,int_sinh2(4,:,j),'-o',1:m_max,int_sinh2(5,:,j),'-o',1:m_max,int_sinh2(6,:,j),'-o',1:m_max,int_sinh2(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\mathrm e^{-m(\pi-\delta)}$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$ ',num2str(d(j)),'$\pi$']);
    colororder(["#808080";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 4.3: Semilogarithmic plots of the term $\frac{\beta}{\sinh \beta}\,\int_{1/\alpha}^1 \frac{I_1(\beta \sqrt{1-\nu^2})}{\sqrt{1-\nu^2}}\,\mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1.1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.');


% Terms for the continuous Kaiser--Bessel window function
figure(6); 
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,decay(1,:,j),'--k',1:m_max,abs(int_cKB3(1,:,j)),'-o',1:m_max,int_cKB3(2,:,j),'-o',1:m_max,int_cKB3(3,:,j),'-o',1:m_max,int_cKB3(4,:,j),'-o',1:m_max,int_cKB3(5,:,j),'-o',1:m_max,int_cKB3(6,:,j),'-o',1:m_max,int_cKB3(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\mathrm e^{-m(\pi-\delta)}$','$\alpha=1$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$ ',num2str(d(j)),'$\pi$']);
    colororder(["#808080";"#E69F00";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 5.3: Semilogarithmic plots of the term $\frac{2\beta}{\pi\,\big(I_0(\beta) - 1\big)}\,\int_{1/\alpha}^{\nu_1(0)}\Big(\frac{\sinh\big(\beta \,\sqrt{1- \nu^2}\big)}{\beta\,\sqrt{1 - \nu^2}} - \frac{\sin(\beta\,\nu)}{\beta\,\nu}\Big)\, \mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1.1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.')
figure(7); 
for j = 1:length(d)
    subplot(1,length(d),j); semilogy(1:m_max,decay(1,:,j),'--k',1:m_max,int_cKB2(2,:,j),'-o',1:m_max,int_cKB2(3,:,j),'-o',1:m_max,int_cKB2(4,:,j),'-o',1:m_max,int_cKB2(5,:,j),'-o',1:m_max,int_cKB2(6,:,j),'-o',1:m_max,int_cKB2(7,:,j),'-o');
    xlim([1;m_max]); xlabel('$m$'); set(gca,'Fontsize',16);
    legend('$\mathrm e^{-m(\pi-\delta)}$','$\alpha=1.1$','$\alpha=1.2$','$\alpha=1.5$','$\alpha=1.8$','$\alpha=2$','$\alpha=3$','Location','southwest');
    title(['$\delta=$ ',num2str(d(j)),'$\pi$']);
    colororder(["#808080";"#CC79A7";"#56B4E9";"#009E73";"#0072B2";"#F0E442";"#D55E00"])
end%for
sgtitle('Figure 5.4: Semilogarithmic plots of the term $\frac{2\beta}{\pi \,\big(I_0(\beta) - 1 \big)}\,\int_{1/\alpha}^1 \Big(\frac{\sinh\big(\beta \,\sqrt{1- \nu^2}\big)}{\beta\,\sqrt{1 - \nu^2}} - \frac{\sin(\beta \nu)}{\beta \nu}\Big)\, \mathrm{d}\nu$ for $m = 1,\ldots,20$, $\alpha \in \{1.1,\,1.2,\,1.5,\,2,\,3\}$, and $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$.')

%% Generate tables for tikz

if (save_results == 1)
fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\frac{\\beta}{\\sinh\\beta}\\int_{1/\\alpha}^1\\frac{I_1(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    fprintf(fileID,['\n\n exp(-m(pi-delta))\n\n']);
    matrix = [1:m_max;decay(1,:,j)];
    fprintf(fileID,format,matrix);
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;int_sinh2(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for

fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\frac{\\beta}{\\sinh\\beta}\\int_{1/\\alpha}^{v_1(0)}\\frac{I_1(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    fprintf(fileID,['\n\n exp(-m(pi-delta))\n\n']);
    matrix = [1:m_max;decay(1,:,j)];
    fprintf(fileID,format,matrix);
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;int_sinh3(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for

fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\frac{2\\beta}{\\pi (I_0(\\beta) - 1)}\\int_{1/\\alpha}^1\\frac{\\sinh(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}-\\frac{\\sin(\\beta t)}{t}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    fprintf(fileID,['\n\n exp(-m(pi-delta))\n\n']);
    matrix = [1:m_max;decay(1,:,j)];
    fprintf(fileID,format,matrix);
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;int_cKB2(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for

fprintf(fileID,'---------------------------------------------------------------\n');
fprintf(fileID,'---------------------------------------------------------------');
fprintf(fileID,'\n\n\\frac{2\\beta}{\\pi (I_0(\\beta) - 1)}\\int_{1/\\alpha}^{v_1(0)}\\frac{\\sinh(\\beta\\sqrt{1-t^2})}{\\sqrt{1-t^2}}-\\frac{\\sin(\\beta t)}{t}\\,\\mathrm dt\n');
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

for j = 1:length(d)
% % Test for different delta
    fprintf(fileID,['delta=',num2str(d(j)),'pi']);
    fprintf(fileID,'\n-----------------------');
    fprintf(fileID,['\n\n exp(-m(pi-delta))\n\n']);
    matrix = [1:m_max;decay(1,:,j)];
    fprintf(fileID,format,matrix);
    for k = 1:length(a)
        fprintf(fileID,['\n\n alpha=',num2str(a(k)),'\n\n']);
        matrix = [1:m_max;int_cKB3(k,:,j)];
        fprintf(fileID,format,matrix);
    end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
end%for
end%if
end%if

if (save_results == 1), fclose(fileID); end%if

%% Nested function

% Sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function

% Struve function
function f=struvem(v,x,n)
% Drew (2024). struvem.m (https://www.mathworks.com/matlabcentral/fileexchange/27218-struvem-m), MATLAB Central File Exchange. Retrieved February 19, 2024.
% 
% Calculates the Modified Struve Function L_v(x) and n is the length of
% the series calculation (n=100 if unspecified)
%
% from: Abramowitz and Stegun: Handbook of Mathematical Functions
% 		http://www.math.sfu.ca/~cbm/aands/page_498.htm
% 
if nargin<3
n=100;
end
k=0:n;
x=x(:)';
k=k(:);
xx=repmat(x,length(k),1);
kk=repmat(k,1,length(x));
TOP=1;
BOT=gamma(kk+1.5).*gamma(kk+v+1.5);
RIGHT=(xx./2).^(2.*kk+v+1);
FULL=TOP./BOT.*RIGHT;
f=sum(FULL);
end%function