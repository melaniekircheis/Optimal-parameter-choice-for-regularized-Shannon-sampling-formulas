% Code file for the Figures 3.1, 4.4, 5.5 and 6.1

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 1;

% Set parameters to study
n = 2:10; % truncation parameter
delta = [pi/4;pi/2;3/4*pi]; % bandwidth

% Initialization of error vectors
err_Gauss_opt = zeros(length(n),length(delta)); err_Gauss_smaller = zeros(length(n),length(delta)); err_Gauss_larger = zeros(length(n),length(delta)); 
err_mGauss_opt = zeros(length(n),length(delta)); err_mGauss_smaller = zeros(length(n),length(delta)); err_mGauss_larger = zeros(length(n),length(delta)); 
err_sinh_opt = zeros(length(n),length(delta)); err_sinh_smaller = zeros(length(n),length(delta)); err_sinh_larger = zeros(length(n),length(delta)); 
err_cKB_opt = zeros(length(n),length(delta)); err_cKB_smaller = zeros(length(n),length(delta)); err_cKB_larger = zeros(length(n),length(delta)); 

% Initialization of vectors for error constants
const_Gauss = zeros(length(n),length(delta)); 
const_mGauss = zeros(length(n),length(delta)); 
const_sinh = zeros(length(n),length(delta)); 
const_cKB = zeros(length(n),length(delta)); 

% Initialization of vectors for computation time
time_Gauss = zeros(length(n),length(delta)); 
time_mGauss = zeros(length(n),length(delta)); 
time_sinh = zeros(length(n),length(delta)); 
time_cKB = zeros(length(n),length(delta)); 

%% Error constants

for k = 1:length(delta)     
    % Gaussian window function
    const_Gauss(:,k) = 2*sqrt(2)./sqrt(pi*n*(pi-delta(k))).*exp(-n*(pi-delta(k))/2);
    % Modified Gaussian window function
    lambda = (pi-delta(k))./2;
    const_mGauss(:,k) = 2*sqrt(2)./sqrt(pi*n*(pi-lambda-delta(k))).*exp(-n*(pi-lambda-delta(k))/2);
    % sinh-type window function
    beta_opt = n.*(pi-delta(k)); % set shape parameter
    const_sinh(:,k) = exp(-beta_opt);
    % Continuous Kaiser-Bessel window function
    const_cKB(:,k) = (7/8*beta_opt+7/pi*beta_opt.^2).*exp(-beta_opt);
end%for

%% Reconstruction error

% Initialization of a fine grid for evaluation of reconstruction error
S = 1e5;
s = (-S:S)';
t = s/S;

% Loop for computation of the error
for k = 1:length(delta)
    % Set test function
%     f = @(x) sqrt(delta(k)/pi)*my_sinc(delta(k),x);
%     f = @(x) sqrt(3*delta(k)/(4*pi))*(my_sinc(delta(k)/2,x)).^2;
    f = @(x) 2*delta(k)/sqrt(5*delta(k)*pi+4*sin(delta(k))*pi)*(my_sinc(delta(k),x)+my_sinc(delta(k),(x-1))/2);

    % Set function evaluations for comparison
    ft = f(t);

    for i2 = 1:length(n) 
        % Set truncation parameter
        m = n(i2);

        % Set oversampling
        T = m+1;
        j = (-T:T)'; % Corresponding index set

        % Set equispaced function evaluations
        fj = f(j);

        % Initialization of vectors for reconstructions
        Rm_Gauss_opt = zeros(length(t),1); Rm_Gauss_smaller = zeros(length(t),1); Rm_Gauss_larger = zeros(length(t),1);
        Rm_mGauss_opt = zeros(length(t),1); Rm_mGauss_smaller = zeros(length(t),1); Rm_mGauss_larger = zeros(length(t),1);
        Rm_sinh_opt = zeros(length(t),1); Rm_sinh_smaller = zeros(length(t),1); Rm_sinh_larger = zeros(length(t),1);
        Rm_cKB_opt = zeros(length(t),1); Rm_cKB_smaller = zeros(length(t),1); Rm_cKB_larger = zeros(length(t),1);

        % Setup
        for i3 = 1:length(t)
            x = t(i3,:)-j.'; % evaluation points
            phi = my_sinc(pi,x); % sinc function
            ind_delta = (abs(x)-m<=eps); % characteristic function
        
            % % Evaluation of Gaussian regularization
            % Using the optimal parameter
            tic; mu_opt = sqrt(m./(pi-delta(k))); % Set variance of Gaussian
            psi_Gauss_opt = phi.*exp(-x.^2/(2*mu_opt.^2)); psi_Gauss_opt(~ind_delta) = 0;
            Rm_Gauss_opt(i3) = psi_Gauss_opt*fj; 
            time_Gauss(i2,k) = time_Gauss(i2,k) + toc;
            % Using a smaller parameter
            mu_smaller = mu_opt/2; % Set variance of Gaussian
            psi_Gauss_smaller = phi.*exp(-x.^2/(2*mu_smaller.^2)); psi_Gauss_smaller(~ind_delta) = 0;
            Rm_Gauss_smaller(i3) = psi_Gauss_smaller*fj;
            % Using a larger parameter
            mu_larger = 2*mu_opt; % Set variance of Gaussian
            psi_Gauss_larger = phi.*exp(-x.^2/(2*mu_larger.^2)); psi_Gauss_larger(~ind_delta) = 0;
            Rm_Gauss_larger(i3) = psi_Gauss_larger*fj;

            % % Evaluation of modified Gaussian regularization
            % Using the optimal parameter
            tic; lambda = (pi-delta(k))./2; mu_opt = sqrt(m./(pi-lambda-delta(k))); % Set variance of Gaussian
            psi_mGauss_opt = phi.*exp(-x.^2/(2*mu_opt.^2)).*cos(lambda*x); psi_mGauss_opt(~ind_delta) = 0;
            Rm_mGauss_opt(i3) = psi_mGauss_opt*fj; 
            time_mGauss(i2,k) = time_mGauss(i2,k) + toc;
            % Using a smaller parameter
            mu_smaller = mu_opt/2; % Set variance of Gaussian
            psi_mGauss_smaller = phi.*exp(-x.^2/(2*mu_smaller.^2)).*cos(lambda*x); psi_mGauss_smaller(~ind_delta) = 0;
            Rm_mGauss_smaller(i3) = psi_mGauss_smaller*fj;
            % Using a larger parameter
            mu_larger = 2*mu_opt; % Set variance of Gaussian
            psi_mGauss_larger = phi.*exp(-x.^2/(2*mu_larger.^2)).*cos(lambda*x); psi_mGauss_larger(~ind_delta) = 0;
            Rm_mGauss_larger(i3) = psi_mGauss_larger*fj;

            % % Evaluation of sinh-type regularization
            % Using the optimal parameter
            tic; beta_opt = m*(pi-delta(k));
            if beta_opt==0
                psi_sinh_opt = phi.*sqrt(1-(x/m).^2);
            else
                psi_sinh_opt = phi.*sinh(beta_opt*sqrt(1-(x/m).^2))/sinh(beta_opt); 
            end
            psi_sinh_opt(~ind_delta) = 0;       
            Rm_sinh_opt(i3) = psi_sinh_opt*fj; 
            time_sinh(i2,k) = time_sinh(i2,k) + toc;
            % Using a smaller parameter
            beta_smaller = 1/2*m*(pi-delta(k));
            if beta_smaller==0
                psi_sinh_smaller = phi.*sqrt(1-(x/m).^2);
            else
                psi_sinh_smaller = phi.*sinh(beta_smaller*sqrt(1-(x/m).^2))/sinh(beta_smaller); 
            end
            psi_sinh_smaller(~ind_delta) = 0;       
            Rm_sinh_smaller(i3) = psi_sinh_smaller*fj;
            % Using a larger parameter
            beta_larger = 2*m*(pi-delta(k));
            if beta_larger==0
                psi_sinh_larger = phi.*sqrt(1-(x/m).^2);
            else
                psi_sinh_larger = phi.*sinh(beta_larger*sqrt(1-(x/m).^2))/sinh(beta_larger); 
            end
            psi_sinh_larger(~ind_delta) = 0;       
            Rm_sinh_larger(i3) = psi_sinh_larger*fj;

            % % Evaluation of continuous Kaiser-Bessel regularization
            % Using the optimal parameter
            tic; beta_opt = m*(pi-delta(k));
            if beta_opt==0
                psi_cKB_opt = phi.*(1-(x/m).^2);
            else
                psi_cKB_opt = phi.*(besseli(0,beta_opt*sqrt(1-(x/m).^2))-1)/(besseli(0,beta_opt)-1);
            end
            psi_cKB_opt(~ind_delta) = 0;
            Rm_cKB_opt(i3) = psi_cKB_opt*fj; 
            time_cKB(i2,k) = time_cKB(i2,k) + toc;
            % Using a smaller parameter
            beta_smaller = 1/2*m*(pi-delta(k));
            if beta_smaller==0
                psi_cKB_smaller = phi.*(1-(x/m).^2);
            else
                psi_cKB_smaller = phi.*(besseli(0,beta_smaller*sqrt(1-(x/m).^2))-1)/(besseli(0,beta_smaller)-1);
            end
            psi_cKB_smaller(~ind_delta) = 0;
            Rm_cKB_smaller(i3) = psi_cKB_smaller*fj; 
            % Using a larger parameter
            beta_larger = 2*m*(pi-delta(k));
            if beta_larger==0
                psi_cKB_larger = phi.*(1-(x/m).^2);
            else
                psi_cKB_larger = phi.*(besseli(0,beta_larger*sqrt(1-(x/m).^2))-1)/(besseli(0,beta_larger)-1);
            end
            psi_cKB_larger(~ind_delta) = 0;
            Rm_cKB_larger(i3) = psi_cKB_larger*fj; 

        end%for

        % Computation of reconstruction errors
        err_Gauss_opt(i2,k) = norm(Rm_Gauss_opt-ft,inf); err_Gauss_smaller(i2,k) = norm(Rm_Gauss_smaller-ft,inf); err_Gauss_larger(i2,k) = norm(Rm_Gauss_larger-ft,inf);
        err_mGauss_opt(i2,k) = norm(Rm_mGauss_opt-ft,inf); err_mGauss_smaller(i2,k) = norm(Rm_mGauss_smaller-ft,inf); err_mGauss_larger(i2,k) = norm(Rm_mGauss_larger-ft,inf);
        err_sinh_opt(i2,k) = norm(Rm_sinh_opt-ft,inf); err_sinh_smaller(i2,k) = norm(Rm_sinh_smaller-ft,inf); err_sinh_larger(i2,k) = norm(Rm_sinh_larger-ft,inf);
        err_cKB_opt(i2,k) = norm(Rm_cKB_opt-ft,inf); err_cKB_smaller(i2,k) = norm(Rm_cKB_smaller-ft,inf); err_cKB_larger(i2,k) = norm(Rm_cKB_larger-ft,inf);
    end%for

    fprintf(['delta= ',num2str(delta(k)/pi),'pi done %s\n'], datestr(datetime('now')))
end%for

%% Visualization 

% Optimality of the parameter for the Gaussian window function
for k = 1:length(delta)
figure(1); subplot(1,3,k); semilogy(n,err_Gauss_smaller(:,k),':^',n,err_Gauss_opt(:,k),'-o',n,err_Gauss_larger(:,k),'--square'); 
legend('$\sigma^2=\frac{m}{4(\pi-\delta)}$','$\sigma^2=\frac{m}{\pi-\delta}$','$\sigma^2=\frac{4m}{\pi-\delta}$','Location','southwest'); 
xlabel('$m$'); title(['$\delta=$ ',num2str(delta(k)/pi),'$\pi$'])
colororder(["#0072B2";"#E69F00";"#F0E442"])
end%for
sgtitle({'Figure 3.1: Maximum approximation error (3.7) using the Gaussian','function $\varphi_{\mathrm{Gauss}}$ in (1.5) with different variances $\sigma^2 \in \big\{\frac{m}{4(\pi-\delta)}, \frac{m}{\pi-\delta}, \frac{4m}{\pi-\delta}\big\}$,','for the bandlimited function (3.8) with bandwidths $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$',' and truncation parameters $m\in\{2, 3, \ldots, 10\}$.'});

% Optimality of the parameter for the modified Gaussian window function
for k = 1:length(delta)
figure(2); subplot(1,3,k); semilogy(n,err_mGauss_smaller(:,k),':^',n,err_mGauss_opt(:,k),'-o',n,err_mGauss_larger(:,k),'--square'); 
legend('$\sigma^2=\frac{m}{4(\pi-\lambda-\delta)}$','$\sigma^2=\frac{m}{\pi-\lambda-\delta}$','$\sigma^2=\frac{4m}{\pi-\lambda-\delta}$','Location','southwest'); 
xlabel('$m$'); title(['$\delta=$ ',num2str(delta(k)/pi),'$\pi$'])
colororder(["#0072B2";"#E69F00";"#F0E442"])
end%for
sgtitle({'Maximum approximation error (3.7) using the modified Gaussian','function $\varphi_{\mathrm{mGauss}}$ in (1.6) with different parameters','$\sigma^2 \in \big\{\frac{m}{4(\pi-\lambda-\delta)}, \frac{m}{\pi-\lambda-\delta}, \frac{4m}{\pi-\lambda-\delta}\big\}$, for the bandlimited function (3.8) with','bandwidths $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$ and truncation parameters $m\in\{2, 3, \ldots, 10\}$.'});

% Optimality of the parameter for the sinh-type window function
for k = 1:length(delta)
figure(3); subplot(1,3,k); semilogy(n,err_sinh_smaller(:,k),':^',n,err_sinh_opt(:,k),'-o',n,err_sinh_larger(:,k),'--square'); 
legend('$\alpha=\frac 12$','$\alpha=1$','$\alpha=2$','Location','southwest'); 
xlabel('$m$'); title(['$\delta=$ ',num2str(delta(k)/pi),'$\pi$'])
colororder(["#0072B2";"#E69F00";"#F0E442"])
end%for
sgtitle({'Figure 4.4: Maximum approximation error (3.7) using the $\sinh$-type','window function $\varphi_{\mathrm{sinh}}$ in (1.7) with different shape parameters','$\beta=\alpha\, m (\pi-\delta)$, $\alpha\in\big\{\frac 12,1,2\big\}$, for the bandlimited function (3.8) with','bandwidths $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$ and truncation parameters $m\in\{2, 3, \ldots, 10\}$.'});

% Optimality of the parameter for the continuous Kaiser-Bessel window function
for k = 1:length(delta)
figure(4); subplot(1,3,k); semilogy(n,err_cKB_smaller(:,k),':^',n,err_cKB_opt(:,k),'-o',n,err_cKB_larger(:,k),'--square'); 
legend('$\alpha=\frac 12$','$\alpha=1$','$\alpha=2$','Location','southwest'); 
xlabel('$m$'); title(['$\delta=$ ',num2str(delta(k)/pi),'$\pi$'])
colororder(["#0072B2";"#E69F00";"#F0E442"])
end%for
sgtitle({'Figure 5.5: Maximum approximation error (3.7) using the','continuous Kaiser--Bessel window function $\varphi_{\mathrm{cKB}}$ in (1.8) with different','shape parameters $\beta=\alpha\, m (\pi-\delta)$, $\alpha\in\big\{\frac 12,1,2\big\}$, for the bandlimited','function (3.8) with bandwidths $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$','and truncation parameters $m\in\{2, 3, \ldots, 10\}$.'});

% Comparison of the different window functions
for k = 1:length(delta)
figure(5); subplot(1,3,k); semilogy(n,err_Gauss_opt(:,k),'-o',n,const_Gauss(:,k),'--',n,err_mGauss_opt(:,k),'-*',n,const_mGauss(:,k),'--',n,err_sinh_opt(:,k),'-diamond',n,const_sinh(:,k),'--',n,err_cKB_opt(:,k),'-square',n,const_cKB(:,k),'--'); 
legend('$\varphi_{\mathrm{Gauss}}$','','$\varphi_{\mathrm{mGauss}}$','','$\varphi_{\sinh}$','','$\varphi_{\mathrm{cKB}}$','','Location','southwest'); 
xlabel('$m$'); title(['$\delta=$ ',num2str(delta(k)/pi),'$\pi$'])
colororder(["#801919";"#801919";"#10520e";"#10520e";"#77AC30";"#77AC30";"#00008B";"#00008B"])
end%for
sgtitle({'Figure 6.1: Maximum approximation error (3.7) (solid) and error constants (dashed)','using $\varphi\in\{\varphi_{\mathrm{Gauss}},\varphi_{\mathrm{mGauss}},\varphi_{\sinh},\varphi_{\mathrm{cKB}}\}$, see (1.5), (1.6), (1.7), and (1.8),','for the function (3.8) with $\delta \in \big\{\frac{\pi}{4},\, \frac{\pi}{2},\,\frac{3\pi}{4}\big\}$ and $m\in\{2, 3, \ldots, 10\}$.'});

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('comparison_optimal_parameters.txt','w');
format = '%d %1.4e \n';

fprintf(fileID,'\n--------------------------\n\n');
fprintf(fileID,'Results using the optimal parameters');
fprintf(fileID,'\n--------------------------\n\n');
% Gaussian window function
fprintf(fileID,'Error for Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_Gauss_opt(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',const_Gauss(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% modified Gaussian window function
fprintf(fileID,'Error for modified Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_mGauss_opt(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for modified Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',const_mGauss(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% sinh-type window function
fprintf(fileID,'Error for sinh-type window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_sinh_opt(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for sinh-type window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',const_sinh(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Continuous Kaiser-Bessel window function
fprintf(fileID,'Error for continuous Kaiser-Bessel window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_cKB_opt(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n--------------\n\n');
fprintf(fileID,'Error constant for continuous Kaiser-Bessel window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',const_cKB(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

fprintf(fileID,'\n--------------------------\n\n');
fprintf(fileID,'Results using smaller parameters');
fprintf(fileID,'\n--------------------------\n\n');
% Gaussian window function
fprintf(fileID,'Error for Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_Gauss_smaller(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% modified Gaussian window function
fprintf(fileID,'Error for modified Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_mGauss_smaller(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% sinh-type window function
fprintf(fileID,'Error for sinh-type window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_sinh_smaller(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Continuous Kaiser-Bessel window function
fprintf(fileID,'Error for continuous Kaiser-Bessel window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_cKB_smaller(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');

fprintf(fileID,'\n--------------------------\n\n');
fprintf(fileID,'Results using larger parameters');
fprintf(fileID,'\n--------------------------\n\n');
% Gaussian window function
fprintf(fileID,'Error for Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_Gauss_larger(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% modified Gaussian window function
fprintf(fileID,'Error for modified Gaussian window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_mGauss_larger(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% sinh-type window function
fprintf(fileID,'Error for sinh-type window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_sinh_larger(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fprintf(fileID,'\n---------------------------------------------------------------\n\n');
% Continuous Kaiser-Bessel window function
fprintf(fileID,'Error for continuous Kaiser-Bessel window function');
for k = 1:length(delta)
    fprintf(fileID,['\n\n delta=',num2str(delta(k)/pi),'pi\n\n']);
    matrix = [n.',err_cKB_larger(:,k)];
    fprintf(fileID,format,matrix.');
end%for
fclose(fileID);
end%if

%% Nested functions

% Sinc function
function y = my_sinc(N,x)
    y = (sin(N*x)./(N*x));
    y(x==0) = 1; 
end%function