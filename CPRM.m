% COMPRESSIBLE PRM
% this code is for testing purposes and trying out new ideas
clc; clear; close all;
%% Inputs
St=15;
maxCompRatio=1;%4;5
deftype='shear';%'axc';%'sc';
ncases=10;%3;
spectrum='SC97';%'VK';
intType = 'quad';%'trapz';%
% low-def case inputs
% S=4.08; or 
% M0=0.4
% Lx=4pi; Ly=Lz=2pi; E(k)=16-32 tophat; sqrt(rho'rho')=sqrt(T'T')=0
% X=0; Re=864; Mrms=0.4; S*0=5.9; ReTau=200.2; Stf=24; 
nts=500; % number 0f time steps
Mt0r=0.25*ones(ncases,1);%[0.025,0.11,0.29];%[0.025,0.025,0.025,0.11,0.11,0.11,0.29,0.29,0.29];%[0.025,0.025,0.025];%
Ret0r=296*ones(ncases,1);%[358,184,500];%[358,358,358,184,184,184,500,500,500];%[358,358,358];%
beta0r=zeros(ncases,1);%[1,2,3];%
gamma=1.4;%5/3;%
R_g=286; % gas constant
T = 300; % mean temperature
nshells=60;
dilatationalRatio=0.08*ones(ncases,1);%[0.06,0.18,0.09];%[0.06,0.06,0.06,0.18,0.18,0.18,0.09,0.09,0.09];%[0,0,0];%[0.0143,]
md0_range=[2.7,4,8.3,12,16.5,24,32,42.7,53.4,66.7];%[5,87,29];%[5,5,5,87,87,87,29,29,29];%[2.7,12,32,66.7];%[0.5,0.5,0.5];% %delta Mach number range
maxk=3*ones(ncases,1);%[1.2,1.4,1.7,1.85,2,2.1,2.3,2.4,2.5,2.6];% % maximum wavevector magnitude range
npr=[1000,900,800,700,600,500,400,300,200,200];%700*ones(ncases,1);%[200 500 800 200 500 800 200 500 1000];%[800,800,800];% range of number of particles
for mm=1:ncases
Md_0=md0_range(mm);
n_particles=npr(mm);
Mt0=Mt0r(mm);
Ret0=Ret0r(mm);
beta0=beta0r(mm);
%% Pre-computations
Gij=zeros(3,3);
if strcmp(intType,'quad')
    GL_pts_wts=readmatrix('GL_pts_wts.dat');
    w=GL_pts_wts(:,2)';
    nshells=length(w);
end
a0=sqrt(gamma*R_g*T);% speed of sound
tke_init=0.5*(Mt0*sqrt(gamma*R_g*T))^2;
tke_d_init = dilatationalRatio(mm)*tke_init;
tke_s_init=tke_init-tke_d_init;
k=logspace(-3,5,500);
if strcmp(spectrum,'SC97')
    % Simone Cambon 1997 Spectrum
    kp=8;
    A = 0.00025907;
    E=A*tke_init*exp(-2*k.^2/(kp^2)).*k.^4;
else
    % Von Karman Spectrum
    k=logspace(-3,4,100);
    Cnuk=0.4843;
    E=2.0*tke_init*Cnuk*(k.^4)./((1.0 + k.^2).^(17/6));
end
eps=2*tke_init*sqrt(2*trapz(k,k.^2.*E)/Ret0);
tau=(2*tke_init)/eps;
Jinv=linspace(1,maxCompRatio,nts)';
if strcmp(deftype,'axc')
    Gamma0 = -Md_0/(tau*Mt0);
    S0=0;
    time = (Jinv.^(-1)-1)./Gamma0;
    dt = time(2:end)-time(1:end-1);
    Gij(1,1)=Gamma0; 
else
    S = Md_0/(tau*Mt0);
    time = linspace(0,St,nts);
    dt = (time(2:end)-time(1:end-1))/S;
    Gij(1,2)=S; %homogeneous shear
    % Gij(1,1)=S; Gij(2,2)=-S; %plane strain
    % Gij(1,1)=S; Gij(2,2)=-S/2; Gij(3,3)=-S/2; %axisymmetric contraction
    % Gij(1,1)=-2*S; Gij(2,2)=S; Gij(3,3)=S; %axisymmetric expansion
end

%% Vector Allocations
R=zeros(length(time),3,3); 
Rs=zeros(length(time),3,3); 
Rs2=zeros(length(time),3,3); 
Rd=zeros(length(time),3,3); 
D11=zeros(length(time),1);
Dd11=zeros(length(time),1); 
Ds11=zeros(length(time),1); 

%% Initialization

[Rhat,Rshat,Rdhat,Rsdhat,Bhat,Bshat,Bdhat,Phihat,k,K] = initVelSpecTen(n_particles,nshells,maxk(mm),tke_init,tke_d_init,a0,gamma,spectrum,intType); %initialize velocity spectrum
k2=k; Phihat2=Phihat;
n=reshape(k(:,1,:),n_particles,3);
n_init=n;
[Rschat]=initClusters(nshells,n_particles,n,tke_s_init,spectrum); %initalize clusters of energy rays

if strcmp(deftype,'sc')
    S0 = Md_0/(tau*Mt0);
    Gij(1,2)=S0;
    time = linspace(0,beta0,200);
    dt = (time(2:end)-time(1:end-1))/S0;
    for tt=2:length(time)
        [Rhat,k]=RK4Solenoidal(Rhat,k,Gij,dt(tt-1),n_particles,nshells);
        [Rschat,n]=RK4Solenoidal_cluster(Rschat,n,Gij,dt(tt-1),n_particles);
    end
    SCratio=0.1;
    Gamma0=-S0/SCratio;
    time = (Jinv.^(-1)-1)./Gamma0; % really t-t0 time
    dt = time(2:end)-time(1:end-1);
    Gij(1,2)=S0; 
    Gij(1,1)=Gamma0; 
end

E_ray=zeros(n_particles,3,3); 
Ed_ray=zeros(n_particles,3,3);
Es_ray=zeros(n_particles,3,3);
D11_ray=zeros(n_particles,1);
Ds11_ray=zeros(n_particles,1); 
Dd11_ray=zeros(n_particles,1);
if strcmp(intType,'quad')
    E_ray=reshape(sum(K.^2.*Rhat.*exp(K).*repmat(w,n_particles,1),2),n_particles,3,3);%Gauss-Laguerre
else
    for np=1:n_particles 
        E_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rhat(np,:,:,:),2); 
    end
end
R(1,:,:)=4*pi*sum(E_ray,1)/n_particles;
Rs(1,:,:)=sum(Rschat,1)/n_particles;
a=a0;
if strcmp(deftype,'axc')
if strcmp(intType,'quad')
    D11_ray=sum(k(:,:,1).^2.*K.^2.*(Rhat(:,:,1,1)+Rhat(:,:,2,2)+Rhat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
    % Ed_ray=sum(K.^2.*Rdhat.*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
    % Es_ray=sum(K.^2.*Rshat.*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
    % Ds11_ray=sum(k2(:,:,1).^2.*K.^2.*(Rshat(:,:,1,1)+Rshat(:,:,2,2)+Rshat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
    % Dd11_ray=sum(k2(:,:,1).^2.*K.^2.*(Rdhat(:,:,1,1)+Rdhat(:,:,2,2)+Rdhat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
else
    for np=1:n_particles 
        D11_ray(np) =trapz(K(np,:),k(np,:,1).^2.*K(np,:).^2.*(Rhat(np,:,1,1)+Rhat(np,:,2,2)+Rhat(np,:,3,3)),2);
        % Ed_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rdhat(np,:,:,:),2);
        % Es_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rshat(np,:,:,:),2);
        % Ds11_ray(np) =trapz(K(np,:),k2(np,:,1).^2.*K(np,:).^2.*(Rshat(np,:,1,1)+Rshat(np,:,2,2)+Rshat(np,:,3,3)),2);
        % Dd11_ray(np) =trapz(K(np,:),k2(np,:,1).^2.*K(np,:).^2.*(Rdhat(np,:,1,1)+Rdhat(np,:,2,2)+Rdhat(np,:,3,3)),2);
    end
end
Rd(1,:,:)=R(1,:,:)-Rs(1,:,:);
Ds11(1)=sum(k(:,1,1).^2.*(Rschat(:,1,1)+Rschat(:,2,2)+Rschat(:,3,3)),1)/n_particles;
D11(1)=4*pi*sum(D11_ray,1)/n_particles;
Dd11(1)=D11(1)-Ds11(1);
% Rd(1,:,:)=4*pi*sum(Ed_ray,1)/n_particles;
% Rs2(1,:,:)=4*pi*sum(Es_ray,1)/n_particles;
% Ds11(1)=4*pi*sum(Ds11_ray,1)/n_particles;
% Dd11(1)=4*pi*sum(Dd11_ray,1)/n_particles;
end

% compare initialized values
% ans1=reshape(R(1,:,:),3,3); 2*tke_init/3
%% time integration
for tt=2:length(time)
    %time-step with numerical integration
    [Rhat,Bhat,Phihat,k,K]=RK4CPRM(Rhat,Bhat,Phihat,k,K,Gij,a,dt(tt-1),n_particles,nshells);
    [Rschat,n]=RK4Solenoidal_cluster(Rschat,n,Gij,dt(tt-1),n_particles); 
    % if strcmp(deftype,'axc')
    %     [Rshat,Rdhat,Rsdhat,Bshat,Bdhat,Phihat2,k2]=RK4Decomposed(Rshat,Rdhat,Rsdhat,Bshat,Bdhat,Phihat2,K,k2,Gij,a,dt(tt-1),n_particles,nshells);
    % end
    % update scalars and states
    if strcmp(deftype,'axc')||strcmp(deftype,'sc')
    a=a0*sqrt((1+Gamma0*time(tt))/((1+Gamma0*time(tt))^gamma));
    Gij(1,1)=Gamma0/(1+Gamma0*time(tt));
    Gij(1,2)=S0/(1+Gamma0*time(tt));
    end
    % numerical integration
    if strcmp(intType,'quad')
        E_ray=reshape(sum(K.^2.*Rhat.*exp(K).*repmat(w,n_particles,1),2),n_particles,3,3);%Gauss-Laguerre
    else
        for np=1:n_particles 
            E_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rhat(np,:,:,:),2);
        end
    end
    R(tt,:,:)=4*pi*sum(E_ray,1)/n_particles;
    Rs(tt,:,:)=sum(Rschat,1)/n_particles;
    if strcmp(deftype,'axc')
        if strcmp(intType,'quad')
            D11_ray=sum(k(:,:,1).^2.*K.^2.*(Rhat(:,:,1,1)+Rhat(:,:,2,2)+Rhat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
            % Ed_ray=sum(K.^2.*Rdhat.*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
            % Es_ray=sum(K.^2.*Rshat.*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
            % Ds11_ray=sum(k2(:,:,1).^2.*K.^2.*(Rshat(:,:,1,1)+Rshat(:,:,2,2)+Rshat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
            % Dd11_ray=sum(k2(:,:,1).^2.*K.^2.*(Rdhat(:,:,1,1)+Rdhat(:,:,2,2)+Rdhat(:,:,3,3)).*exp(K).*repmat(w,n_particles,1),2);%Gauss-Laguerre
        else
            for np=1:n_particles 
                D11_ray(np) =trapz(K(np,:),k(np,:,1).^2.*K(np,:).^2.*(Rhat(np,:,1,1)+Rhat(np,:,2,2)+Rhat(np,:,3,3)),2);
                % Ed_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rdhat(np,:,:,:),2);
                % Es_ray(np,:,:) =trapz(K(np,:),K(np,:).^2.*Rshat(np,:,:,:),2);
                % Ds11_ray(np) =trapz(K(np,:),k2(np,:,1).^2.*K(np,:).^2.*(Rshat(np,:,1,1)+Rshat(np,:,2,2)+Rshat(np,:,3,3)),2);
                % Dd11_ray(np) =trapz(K(np,:),k2(np,:,1).^2.*K(np,:).^2.*(Rdhat(np,:,1,1)+Rdhat(np,:,2,2)+Rdhat(np,:,3,3)),2);
            end
        end
        Rd(tt,:,:)=R(tt,:,:)-Rs(tt,:,:);
        Ds11(tt)=sum(k(:,1,1).^2.*(Rschat(:,1,1)+Rschat(:,2,2)+Rschat(:,3,3)),1)/n_particles;
        D11(tt)=4*pi*sum(D11_ray,1)/n_particles;
        Dd11(tt)=D11(tt)-Ds11(tt);
        % Rd(tt,:,:)=4*pi*sum(Ed_ray,1)/n_particles;
        % Rs2(tt,:,:)=4*pi*sum(Es_ray,1)/n_particles;
        % Ds11(tt)=4*pi*sum(Ds11_ray,1)/n_particles;
        % Dd11(tt)=4*pi*sum(Dd11_ray,1)/n_particles;
    end

end
tke=0.5*(R(:,1,1)+R(:,2,2)+R(:,3,3));
tke_s=0.5*(Rs(:,1,1)+Rs(:,2,2)+Rs(:,3,3));
if strcmp(deftype,'axc')
    tke_d=0.5*(Rd(:,1,1)+Rd(:,2,2)+Rd(:,3,3));
end
% Energy spectrum
E=2*pi*sum(K.^2.*(Rhat(:,:,1,1)+Rhat(:,:,2,2)+Rhat(:,:,3,3)),1)./n_particles;

%% Plotting

% % plots for multiple Mach numbers
if strcmp(deftype,'axc')
figure (1)
hold on;
% M0=sqrt(2*Jinv./(2.4+Jinv-1.4*Jinv));
% plot(M0,2*tke/(2*tke(1)),'--k');
plot(Jinv,2*tke/(2*tke(1)),'--k');

figure (2)
hold on;
plot(Jinv,tke_s/tke(1),'--k');
plot(Jinv,tke_d/tke(1),'-.k');

figure (3)
hold on;
plot(Jinv,0.5*Ds11./tke_s,'--k');
plot(Jinv,0.5*Dd11./tke_d,'-.k');
% plot(Jinv,0.5*Rd(:,1,1)./tke_d,'-.k');

elseif strcmp(deftype,'sc')
figure (1)
hold on;
% plot(Jinv,2*tke/(2*tke(1)),'--');
plot(Jinv,2*tke_s/(2*tke_s(1)),'-');

figure (3)
hold on;
% plot(Jinv,R(:,1,2)./(2*tke),'--');
plot(Jinv,Rs(:,1,2)./(2*tke_s),'-');

else
figure (1)
hold on;
plot(time,2*tke/(2*tke_init),'--k');

figure (2)
hold on;
plot(time(2:end-1),(tke(3:end)-tke(1:end-2))./(S*tke(2:end-1)*2*dt(1)),'--k');

figure (3)
hold on;
plot(time,-R(:,1,2)./tke,'--k');
if mm==1; plot(time,-Rs(:,1,2)./tke_s,'-k'); end %plot solenoidal limit for lowest gradient Mach number

figure (4)
[mag_k,idx]=sort(sum(K,1)/n_particles);
loglog(mag_k,E(idx));
ylabel('$E(k,t_f)$','Interpreter','latex','FontSize',14);
xlabel('$k$','Interpreter','latex','FontSize',14);

figure (5)
hold on;
plot(time,R(:,1,1)./(2*tke),'--r',time,R(:,2,2)./(2*tke),'--b',time,R(:,3,3)./(2*tke),'--k');
% plot(time,R(:,1,2)./(2*tke),'--','color','#7E2F8E');
if mm==1 
    plot(time,Rs(:,1,1)./(2*tke_s),'-r',time,Rs(:,2,2)./(2*tke_s),'-b',time,Rs(:,3,3)./(2*tke_s),'-k'); 
    % plot(time,Rs(:,1,3)./(2*tke_s),'-','color','#7E2F8E');
end %plot solenoidal limit for lowest gradient Mach number
box on;

figure (5+mm)
[x,y,z] = sphere(100);
subplot(1,2,1)
ray_integrand = griddata(k(:,1,1),k(:,1,2),k(:,1,3),(E_ray(:,1,1)+E_ray(:,2,2)+E_ray(:,3,3)),x,y,z,'nearest');
surf(x,y,z,ray_integrand);
axis square
colormap jet;
shading interp;
cb = colorbar('southoutside');
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\Phi^{\nearrow k}_{kk}$';
cb.Label.FontSize = 14;
subplot(1,2,2)
ray_integrand = griddata(k(:,1,1),k(:,1,2),k(:,1,3),E_ray(:,1,2),x,y,z,'nearest');
surf(x,y,z,ray_integrand);
axis square
colormap jet;
shading interp;
cb = colorbar('southoutside');
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\Phi^{\nearrow k}_{12}$';
cb.Label.FontSize = 14;
cb.Ruler.Exponent = 3;

end

% % DNS data of Blaisdell
% close all;
% B_dat=readmatrix('blaisdell_data.dat');
% figure (1)
% hold on;
% plot(time,R(:,1,1)./(2*tke),time,R(:,2,2)./(2*tke),time,R(:,3,3)./(2*tke),time,R(:,1,2)./(2*tke));
% plot(S*B_dat(:,1),B_dat(:,2),'o');
% ylabel('$\overline{\rho u_i'' u_j''}$','Interpreter','latex','FontSize',14);
% xlabel('$St$','Interpreter','latex','FontSize',14);
% legend('$r_{11}$','$r_{22}$','$r_{33}$','$r_{12}$','Interpreter','latex','FontSize',14,'location','eastoutside');

end
% plot outside loop
if strcmp(deftype,'axc')
    axc_dat=readmatrix('axc_data.dat');
    LIA=readmatrix('ribnerLIA.dat');
    figure (1)
    plot(Jinv,(2+Jinv.^2)/3,'-k');%pressure-released limit
    plot(Jinv,0.5*(1+Jinv.^2.*atan(sqrt(Jinv.^2-1))./(sqrt(Jinv.^2-1))),'-k');% solenoidal
    plot(axc_dat(1:9,1),axc_dat(1:9,2),'ok');
    plot(axc_dat(10:19,1),axc_dat(10:19,2),'ob');
    plot(axc_dat(20:28,1),axc_dat(20:28,2),'or');
    % h = get(gca, 'Children');
    % set(h(8), 'Color', 'b','LineStyle',':'); set(h(7), 'Color', 'b','LineStyle','-.'); set(h(6), 'Color', 'b','LineStyle','--');
    % set(h(11), 'Color', 'r','LineStyle',':'); set(h(10), 'Color', 'r','LineStyle','-.'); set(h(9), 'Color', 'r','LineStyle','--');
    % set(h(14), 'Color', 'k','LineStyle',':'); set(h(13), 'Color', 'k','LineStyle','-.'); set(h(12), 'Color', 'k','LineStyle','--');
    % plot(M0,(2+Jinv.^2)/3,'-k');%pressure-released limit
    % plot(M0,0.5*(1+Jinv.^2.*atan(sqrt(Jinv.^2-1))./(sqrt(Jinv.^2-1))),'-k');% solenoidal
    % plot(sqrt(2*axc_dat(1:9,1)./(2.4+axc_dat(1:9,1)-1.4*axc_dat(1:9,1))),axc_dat(1:9,2),'ok');
    % plot(sqrt(2*axc_dat(10:19,1)./(2.4+axc_dat(10:19,1)-1.4*axc_dat(10:19,1))),axc_dat(10:19,2),'ob');
    % plot(sqrt(2*axc_dat(20:28,1)./(2.4+axc_dat(20:28,1)-1.4*axc_dat(20:28,1))),axc_dat(20:28,2),'or')
    % h = get(gca, 'Children');
    % set(h(6), 'Color', 'b'); set(h(7), 'Color', 'r');set(h(8), 'Color', 'k');
    % plot(LIA(1:33,1),LIA(1:33,2),'-.k');
    % plot(LIA(34:68,1),LIA(34:68,2),':k');
    ylim([0,12]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$q^2(t)/q^2(0)$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);
    % xlabel('$M_0$','Interpreter','latex','FontSize',14);
    % legend([h(14) h(13) h(12)],'$N_c=200$','$N_c=500$','$N_c=800$','Interpreter','latex','FontSize',14,'location','northwest');
    % h = get(gca, 'Children');
    % legend([h(1) h(2) h(5) h(6) h(10)],'Nearfield LIA','Farfield LIA','Homogeneous DNS','Limits','CPRM','Interpreter','latex','FontSize',13,'location','northwest');

    figure (2)
    h = get(gca, 'Children');
    set(h(1), 'Color', 'b'); set(h(2), 'Color', 'b');
    set(h(3), 'Color', 'r'); set(h(4), 'Color', 'r');
    set(h(6), 'Color', 'k'); set(h(5), 'Color', 'k');
    plot(axc_dat(29:42,1),axc_dat(29:42,2),'ok');
    plot(axc_dat(43:53,1),axc_dat(43:53,2),'or');
    plot(axc_dat(54:79,1),axc_dat(54:79,2),'ob');
    ylim([0,4]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$q_{s,d}^2(t)/q^2(0)$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);
     legend('solenoidal','dilatational','Interpreter','latex','FontSize',14,'location','northwest');

     figure (3)
    h = get(gca, 'Children');
    set(h(1), 'Color', 'b'); set(h(2), 'Color', 'b');
    set(h(3), 'Color', 'r'); set(h(4), 'Color', 'r');
    set(h(6), 'Color', 'k'); set(h(5), 'Color', 'k');
    plot(axc_dat(80:95,1),axc_dat(80:95,2),'ok');
    plot(axc_dat(96:111,1),axc_dat(96:111,2),'or');
    plot(axc_dat(112:127,1),axc_dat(112:127,2),'ob');
    ylim([0,1]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$D_{11}^{s,d}/q_{s,d}^2$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);
     legend('solenoidal','dilatational','Interpreter','latex','FontSize',14,'location','northwest');
elseif strcmp(deftype,'sc')
    sc_data=readmatrix('sc_data.dat');
    figure (1)
    h = get(gca, 'Children');
    set(h(1), 'Color', 'b'); set(h(2), 'Color', 'k'); set(h(3), 'Color', 'r');
    plot(sc_data(1:9,1),sc_data(1:9,2),'ob');
    plot(sc_data(10:16,1),sc_data(10:16,2),'ok');
    plot(sc_data(17:22,1),sc_data(17:22,2),'or');
    ylim([1,6]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$q_s^2(t)/q^2(0)$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);

    figure (2)
    hold on;
    % plot(Jinv,R(:,1,1)./(2*tke),'-b');
    % plot(Jinv,R(:,2,2)./(2*tke),'-k');
    % plot(Jinv,R(:,3,3)./(2*tke),'-r');
    plot(Jinv,Rs(:,1,1)./(2*tke_s),'-b');
    plot(Jinv,Rs(:,2,2)./(2*tke_s),'-k');
    plot(Jinv,Rs(:,3,3)./(2*tke_s),'-r');
    plot(sc_data(23:29,1),sc_data(23:29,2),'ob')
    plot(sc_data(30:36,1),sc_data(30:36,2),'or')
    plot(sc_data(37:43,1),sc_data(37:43,2),'ok')
    ylim([0,1]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$b^s_{ij}$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);
    % legend('$b_{11}$','$b_{22}$','$b_{33}$','Interpreter','latex','FontSize',14,'location','northeast');

    figure (3)
    h = get(gca, 'Children');
    set(h(1), 'Color', 'b'); set(h(2), 'Color', 'k'); set(h(3), 'Color', 'r');
    plot(sc_data(44:53,1),sc_data(44:53,2),'ok')
    plot(sc_data(54:61,1),sc_data(54:61,2),'ob')
    plot(sc_data(62:end,1),sc_data(62:end,2),'or')
    ylim([-0.2,0.1]);
    xlim([1,maxCompRatio]);
    box on;
    ylabel('$b^s_{12}$','Interpreter','latex','FontSize',16);
    xlabel('$\rho(t)/\rho(0)$','Interpreter','latex','FontSize',14);

elseif strcmp(deftype,'shear')
    figure (1)
    hold on;
    plot(time,1+time.^2/3,'-k');
    box on;
    pbaspect([1 1 1])
    ax = gca;
    xticks(ax, 0:5:15);
    xlim(ax, [0 St])
    ax.XMinorTick = 'on';
    yticks(ax, 0:2:10);
    ylim(ax, [0 10])
    ax.YMinorTick = 'on';
    ax.FontSize = 8;
    ylabel('$q^2(t)/q^2(0)$','Interpreter','latex','FontSize',12);
    xlabel('$St$','Interpreter','latex','FontSize',10);
    
    figure (2)
    plot(time,2*time./(3+time.^2),'-k');
    ylabel('$\Lambda$','Interpreter','latex','FontSize',16);
    xlabel('$St$','Interpreter','latex','FontSize',14);
    box on;

    % run pressure-released assuming last set is highest Mach number
    R_pr=zeros(nts,3,3);
    [Rhat]=initClusters(nshells,n_particles,n_init,tke_init,spectrum); %initalize clusters of energy rays
    R_pr(1,:,:)=sum(Rhat,1)/n_particles;
    for tt=2:length(time)
        [Rhat]=RK4PressureReleased(Rhat,Gij,dt(tt-1),n_particles);
        R_pr(tt,:,:)=sum(Rhat,1)/n_particles;
    end
    tke_pr=0.5*(R_pr(:,1,1)+R_pr(:,2,2)+R_pr(:,3,3));
    
    figure (3)
    hold on;
    plot(time,-R_pr(:,1,2)./tke_pr,'-k');
    % GR=readmatrix('girimajiResults.dat');
    % for ii=1:11
    %     plot(GR(:,2*ii-1),GR(:,2*ii),'--b');
    % end
    box on;
    pbaspect([1 1 1])
    ax = gca;
    xticks(ax, 0:5:15);
    xlim(ax, [0 St])
    ax.XMinorTick = 'on';
    yticks(ax, 0:0.2:0.6);
    ylim(ax, [0 0.6])
    ax.YMinorTick = 'on';
    ax.FontSize = 8;
    ylabel('$-2b_{12}$','Interpreter','latex','FontSize',12);
    xlabel('$St$','Interpreter','latex','FontSize',10);

    figure (5)
    hold on;
    plot(time,R_pr(:,1,1)./(2*tke_pr),'-r',time,R_pr(:,2,2)./(2*tke_pr),'-b',time,R_pr(:,3,3)./(2*tke_pr),'-k');
    % plot(time,R_pr(:,1,2)./(2*tke_pr),'-', 'Color','#7E2F8E');
    % plot(time,R_b(:,1,1)./tke_b(1),'-r',time,R_b(:,2,2)./tke_b(1),'-b',time,R_b(:,3,3)./tke_b(1),'-k');
    ylabel('$b_{ij}$','Interpreter','latex','FontSize',16);
    % ylabel('$2R_{ij}(t)/q^2(t_0)$','Interpreter','latex','FontSize',16);
    xlabel('$St$','Interpreter','latex','FontSize',14);
    % legend('$b_{11}$','$b_{22}$','$b_{33}$','Interpreter','latex','FontSize',14,'location','eastoutside');
    box on;


end

%%--------------------END OF MAIN CODE-----------------------------------%%

%% Time integration functions

function [Rhat,Bhat,Phihat,k,K]=RK4CPRM(Rhat,Bhat,Phihat,k,K,Gij,a,dt,n_particles,nshells)
    Rhat_0=Rhat; Bhat_0=Bhat; k_0=k; Phihat_0=Phihat; %K_0=K;

    % calculate RHS
    [rhsR1,rhsk1, rhsK1, rhsB1, rhsPhi1] = calcRHSTerms(Gij,Rhat,k,K,Bhat,Phihat,a,n_particles,nshells);
    % K=K+dt*rhsK1/2;
    % update next time steps
    Rhat=Rhat+dt*rhsR1/2;
    k=k+dt*rhsk1/2;
    Bhat=Bhat+dt*rhsB1/2;
    Phihat=Phihat+dt*rhsPhi1/2;

    % calculate RHS2
    [rhsR2,rhsk2, rhsK2, rhsB2, rhsPhi2] = calcRHSTerms(Gij,Rhat,k,K,Bhat,Phihat,a,n_particles,nshells);
    % K=K_0+dt*rhsK2/2;
    % update next time steps
    Rhat=Rhat_0+dt*rhsR2/2;
    k=k_0+dt*rhsk2/2;
    Bhat=Bhat_0+dt*rhsB2/2;
    Phihat=Phihat_0+dt*rhsPhi2/2;

    % calculate RHS3
    [rhsR3,rhsk3, rhsK3, rhsB3, rhsPhi3] = calcRHSTerms(Gij,Rhat,k,K,Bhat,Phihat,a,n_particles,nshells);
    % K=K_0+dt*rhsK3;
    % update next time steps
    Rhat=Rhat_0+dt*rhsR3;
    k=k_0+dt*rhsk3;
    Bhat=Bhat_0+dt*rhsB3;
    Phihat=Phihat_0+dt*rhsPhi3;

    % calculate RHS4
    [rhsR4,rhsk4, rhsK4, rhsB4, rhsPhi4] = calcRHSTerms(Gij,Rhat,k,K,Bhat,Phihat,a,n_particles,nshells);
    % K=K_0+dt*(rhsK1+2*rhsK2+2*rhsK3+rhsK4)/6;
    % update next time steps
    Rhat=Rhat_0+dt*(rhsR1+2*rhsR2+2*rhsR3+rhsR4)/6;
    k=k_0+dt*(rhsk1+2*rhsk2+2*rhsk3+rhsk4)/6;
    Bhat=Bhat_0+dt*(rhsB1+2*rhsB2+2*rhsB3+rhsB4)/6;
    Phihat=Phihat_0+dt*(rhsPhi1+2*rhsPhi2+2*rhsPhi3+rhsPhi4)/6;
end

function [Rhat,n]=RK4Solenoidal_cluster(Rhat,n,Gij,dt,n_particles)
    Rhat_0=Rhat;n_0=n;

    [rhsR1,rhsn1]=calcRhsSolenoidal_cluster(Gij,Rhat,n,n_particles);
    Rhat=Rhat_0+dt*rhsR1/2;
    n=n_0+dt*rhsn1/2;

    [rhsR2,rhsn2]=calcRhsSolenoidal_cluster(Gij,Rhat,n,n_particles);
    Rhat=Rhat_0+dt*rhsR2/2;
    n=n_0+dt*rhsn2/2;

    [rhsR3,rhsn3]=calcRhsSolenoidal_cluster(Gij,Rhat,n,n_particles);
    Rhat=Rhat_0+dt*rhsR3;
    n=n_0+dt*rhsn3;

    [rhsR4,rhsn4]=calcRhsSolenoidal_cluster(Gij,Rhat,n,n_particles);
    Rhat=Rhat_0+dt*(rhsR1+2*rhsR2+2*rhsR3+rhsR4)/6;
    n=n_0+dt*(rhsn1+2*rhsn2+2*rhsn3+rhsn4)/6;

end

function [Rhat,n]=RK4Solenoidal(Rhat,n,Gij,dt,n_particles,nshells)
    Rhat_0=Rhat;n_0=n;

    [rhsR1,rhsn1]=calcRhsSolenoidal(Gij,Rhat,n,n_particles,nshells);
    Rhat=Rhat_0+dt*rhsR1/2;
    n=n_0+dt*rhsn1/2;

    [rhsR2,rhsn2]=calcRhsSolenoidal(Gij,Rhat,n,n_particles,nshells);
    Rhat=Rhat_0+dt*rhsR2/2;
    n=n_0+dt*rhsn2/2;

    [rhsR3,rhsn3]=calcRhsSolenoidal(Gij,Rhat,n,n_particles,nshells);
    Rhat=Rhat_0+dt*rhsR3;
    n=n_0+dt*rhsn3;

    [rhsR4,rhsn4]=calcRhsSolenoidal(Gij,Rhat,n,n_particles,nshells);
    Rhat=Rhat_0+dt*(rhsR1+2*rhsR2+2*rhsR3+rhsR4)/6;
    n=n_0+dt*(rhsn1+2*rhsn2+2*rhsn3+rhsn4)/6;
end

function [Rhat]=RK4PressureReleased(Rhat,Gij,dt,n_particles)
    Rhat_0=Rhat;

    rhsR1 = calcRHSPressureReleased(Gij,Rhat,n_particles);
    Rhat=Rhat_0+dt*rhsR1/2;

    rhsR2 = calcRHSPressureReleased(Gij,Rhat,n_particles);
    Rhat=Rhat_0+dt*rhsR2/2;

    rhsR3 = calcRHSPressureReleased(Gij,Rhat,n_particles);
    Rhat=Rhat_0+dt*rhsR3;

    rhsR4 = calcRHSPressureReleased(Gij,Rhat,n_particles);
    Rhat=Rhat_0+dt*(rhsR1+2*rhsR2+2*rhsR3+rhsR4)/6;
end

function [Rshat,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k]=RK4Decomposed(Rshat,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,K,k,Gij,a,dt,n_particles,nshells)
    Rshat_0=Rshat;Rdhat_0=Rdhat;Rsdhat_0=Rsdhat;Bshat_0=Bshat;Bdhat_0=Bdhat;Phihat_0=Phihat;k_0=k;

    [rhsRs1,rhsk1]=calcRhsSolenoidal(Gij,Rshat,k,n_particles,nshells);
    [rhsRd1]=calcRhsDilatationalRij(Gij,Rdhat,Rsdhat,Bdhat,k,K,a,n_particles,nshells);
    [rhsRsd1]=calcRhsCrossRij(Gij,Rshat,Rsdhat,Bshat,k,K,a,n_particles,nshells);
    [rhsBs1]=calcRhsBs(Gij,Rshat,Rsdhat,Bshat,k,K,n_particles,nshells);
    [rhsBd1]=calcRhsBd(Gij,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k,K,a,n_particles,nshells);
    [rhsPhi1, ~]=calcRhsPhiK(Gij,Bshat+Bdhat,Phihat,k,K,n_particles,nshells);
    Rshat=Rshat_0+dt*rhsRs1/2;
    Rdhat=Rdhat_0+dt*rhsRd1/2;
    Rsdhat=Rsdhat_0+dt*rhsRsd1/2;
    Bshat=Bshat_0+dt*rhsBs1/2;
    Bdhat=Bdhat_0+dt*rhsBd1/2;
    Phihat=Phihat_0+dt*rhsPhi1/2;
    k=k_0+dt*rhsk1/2;

    [rhsRs2,rhsk2]=calcRhsSolenoidal(Gij,Rshat,k,n_particles,nshells);
    [rhsRd2]=calcRhsDilatationalRij(Gij,Rdhat,Rsdhat,Bdhat,k,K,a,n_particles,nshells);
    [rhsRsd2]=calcRhsCrossRij(Gij,Rshat,Rsdhat,Bshat,k,K,a,n_particles,nshells);
    [rhsBs2]=calcRhsBs(Gij,Rshat,Rsdhat,Bshat,k,K,n_particles,nshells);
    [rhsBd2]=calcRhsBd(Gij,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k,K,a,n_particles,nshells);
    [rhsPhi2, ~]=calcRhsPhiK(Gij,Bshat+Bdhat,Phihat,k,K,n_particles,nshells);
    Rshat=Rshat_0+dt*rhsRs2/2;
    Rdhat=Rdhat_0+dt*rhsRd2/2;
    Rsdhat=Rsdhat_0+dt*rhsRsd2/2;
    Bshat=Bshat_0+dt*rhsBs2/2;
    Bdhat=Bdhat_0+dt*rhsBd2/2;
    Phihat=Phihat_0+dt*rhsPhi2/2;
    k=k_0+dt*rhsk2/2;

    [rhsRs3,rhsk3]=calcRhsSolenoidal(Gij,Rshat,k,n_particles,nshells);
    [rhsRd3]=calcRhsDilatationalRij(Gij,Rdhat,Rsdhat,Bdhat,k,K,a,n_particles,nshells);
    [rhsRsd3]=calcRhsCrossRij(Gij,Rshat,Rsdhat,Bshat,k,K,a,n_particles,nshells);
    [rhsBs3]=calcRhsBs(Gij,Rshat,Rsdhat,Bshat,k,K,n_particles,nshells);
    [rhsBd3]=calcRhsBd(Gij,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k,K,a,n_particles,nshells);
    [rhsPhi3, ~]=calcRhsPhiK(Gij,Bshat+Bdhat,Phihat,k,K,n_particles,nshells);
    Rshat=Rshat_0+dt*rhsRs3;
    Rdhat=Rdhat_0+dt*rhsRd3;
    Rsdhat=Rsdhat_0+dt*rhsRsd3;
    Bshat=Bshat_0+dt*rhsBs3;
    Bdhat=Bdhat_0+dt*rhsBd3;
    Phihat=Phihat_0+dt*rhsPhi3;
    k=k_0+dt*rhsk3;

    [rhsRs4,rhsk4]=calcRhsSolenoidal(Gij,Rshat,k,n_particles,nshells);
    [rhsRd4]=calcRhsDilatationalRij(Gij,Rdhat,Rsdhat,Bdhat,k,K,a,n_particles,nshells);
    [rhsRsd4]=calcRhsCrossRij(Gij,Rshat,Rsdhat,Bshat,k,K,a,n_particles,nshells);
    [rhsBs4]=calcRhsBs(Gij,Rshat,Rsdhat,Bshat,k,K,n_particles,nshells);
    [rhsBd4]=calcRhsBd(Gij,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k,K,a,n_particles,nshells);
    [rhsPhi4, ~]=calcRhsPhiK(Gij,Bshat+Bdhat,Phihat,k,K,n_particles,nshells);
    Rshat=Rshat_0+dt*(rhsRs1+2*rhsRs2+2*rhsRs3+rhsRs4)/6;
    Rdhat=Rdhat_0+dt*(rhsRd1+2*rhsRd2+2*rhsRd3+rhsRd4)/6;
    Rsdhat=Rsdhat_0+dt*(rhsRsd1+2*rhsRsd2+2*rhsRsd3+rhsRsd4)/6;
    Bshat=Bshat_0+dt*(rhsBs1+2*rhsBs2+2*rhsBs3+rhsBs4)/6;
    Bdhat=Bdhat_0+dt*(rhsBd1+2*rhsBd2+2*rhsBd3+rhsBd4)/6;
    Phihat=Phihat_0+dt*(rhsPhi1+2*rhsPhi2+2*rhsPhi3+rhsPhi4)/6;
    k=k_0+dt*(rhsk1+2*rhsk2+2*rhsk3+rhsk4)/6;
end

%% RHS functions

function [rhsR,rhsk, rhsK, rhsB, rhsPhi] = calcRHSTerms(Gij,Rhat,k,K,Bhat,Phihat,a0,n_particles,nshells)
rhsR=zeros(n_particles,nshells,3,3);
rhsk=zeros(n_particles,nshells,3);
rhsB=zeros(n_particles,nshells,3);
[rhsPhi,rhsK]=calcRhsPhiK(Gij,Bhat,Phihat,k,K,n_particles,nshells);
     for ii=1:3
        [rhsk(:,:,ii), rhsB(:,:,ii)]=calcRhsBi(ii,Gij,Rhat,Bhat,Phihat,k,K,a0,n_particles,nshells);
        for jj=1:3
            rhsR(:,:,ii,jj)=calcRhsRij(ii,jj,Gij,Rhat,Bhat,k,K,a0,n_particles,nshells);
        end
     end
     % symmetrize
    rhsR(:,:,2,1)=rhsR(:,:,1,2);
    rhsR(:,:,3,1)=rhsR(:,:,1,3);
    rhsR(:,:,3,2)=rhsR(:,:,2,3);
end


function [rhsR]=calcRhsRij(ii,jj,Gij,Rhat,Bhat,k,K,a,n_particles,nshells)
    rhsR=zeros(n_particles,nshells,1); 
    for kk=1:3
        rhsR=rhsR-Gij(ii,kk)*Rhat(:,:,kk,jj)-Gij(jj,kk)*Rhat(:,:,kk,ii);%+2*Gij(kk,kk)*Rhat(:,ii,jj);
    end
    rhsR=rhsR+K.*k(:,:,ii).*Bhat(:,:,jj)*a^2+K.*k(:,:,jj).*Bhat(:,:,ii)*a^2; %negative OG
    %solenoidal component
    % for kk=1:3
    %     rhsRs=rhsRs-Gij(ii,kk)*Rshat(:,:,kk,jj)-Gij(jj,kk)*Rshat(:,:,kk,ii);
    %     for mm=1:3
    %         rhsRs=rhsRs+2*Gij(kk,mm)*(Rshat(:,:,ii,mm).*k(:,:,kk).*k(:,:,jj)+Rshat(:,:,jj,mm).*k(:,:,kk).*k(:,:,ii)); 
    %     end
    % end
end

function [rhsk, rhsB]=calcRhsBi(ii,Gij,Rhat,Bhat,Phihat,k,K,a,n_particles,nshells)
    rhsB=zeros(n_particles,nshells,1);
    rhsk=zeros(n_particles,nshells,1);
    for kk=1:3
        rhsB=rhsB-K.*k(:,:,kk).*Rhat(:,:,kk,ii)-Gij(ii,kk)*Bhat(:,:,kk);%first term positive OG
        rhsk=rhsk-Gij(kk,ii)*k(:,:,kk);
        for mm=1:3
            rhsk=rhsk+Gij(kk,mm)*k(:,:,kk).*k(:,:,mm).*k(:,:,ii);
        end
    end
    rhsB=rhsB+K.*k(:,:,ii).*Phihat*a^2; %negative OG
end

function [rhsPhi, rhsK]=calcRhsPhiK(Gij,Bhat,Phihat,k,K,n_particles,nshells)
    rhsPhi=zeros(n_particles,nshells); 
    rhsK=zeros(n_particles,nshells,1);
    for kk=1:3
        rhsPhi=rhsPhi-2*K.*k(:,:,kk).*Bhat(:,:,kk);%positive OG
        for nn=1:3
            rhsK = rhsK-K.*Gij(kk,nn).*k(:,:,kk).*k(:,:,nn);
        end
    end

end

%% Special cases
function rhsR = calcRHSPressureReleased(Gij,Rhat,n_particles)
rhsR=zeros(n_particles,3,3);
    for ii=1:3
        for jj=1:3
            rhs=zeros(n_particles,1);
            for kk=1:3
                rhs=rhs-Gij(ii,kk)*Rhat(:,kk,jj)-Gij(jj,kk)*Rhat(:,kk,ii);%+2*Gij(kk,kk)*Rhat(:,ii,jj);
            end
            rhsR(:,ii,jj)=rhs;
        end
    end
end

function [rhsR,rhsn]=calcRhsSolenoidal_cluster(Gij,Rhat,n,n_particles)
    rhsR=zeros(n_particles,3,3); 
    rhsn=zeros(n_particles,3); 
     for ii=1:3
       for kk=1:3
            rhsn(:,ii)=rhsn(:,ii)-Gij(kk,ii)*n(:,kk);
            for mm=1:3
                rhsn(:,ii)=rhsn(:,ii)+Gij(kk,mm)*n(:,kk).*n(:,mm).*n(:,ii);
            end
        end
        for jj=1:3
            rhs=zeros(n_particles,1);
            for kk=1:3
                rhs=rhs-Gij(ii,kk)*Rhat(:,kk,jj)-Gij(jj,kk)*Rhat(:,kk,ii);%+2*Gij(kk,kk)*Rhat(:,ii,jj);
                for mm=1:3
                    rhs=rhs+2*Gij(kk,mm)*(Rhat(:,ii,mm).*n(:,kk).*n(:,jj)+Rhat(:,jj,mm).*n(:,kk).*n(:,ii)); 
                end
            end
            rhsR(:,ii,jj)=rhs;
        end
    end
end

function [rhsR,rhsn]=calcRhsSolenoidal(Gij,Rhat,n,n_particles,nshells)
    rhsR=zeros(n_particles,nshells,3,3); 
    rhsn=zeros(n_particles,nshells,3); 
     for ii=1:3
       for kk=1:3
            rhsn(:,:,ii)=rhsn(:,:,ii)-Gij(kk,ii)*n(:,:,kk);
            for mm=1:3
                rhsn(:,:,ii)=rhsn(:,:,ii)+Gij(kk,mm)*n(:,:,kk).*n(:,:,mm).*n(:,:,ii);
            end
        end
        for jj=1:3
            for kk=1:3
                rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)-Gij(ii,kk)*Rhat(:,:,kk,jj)-Gij(jj,kk)*Rhat(:,:,kk,ii);
                for mm=1:3
                    rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)+2*Gij(kk,mm)*(Rhat(:,:,ii,mm).*n(:,:,kk).*n(:,:,jj)+Rhat(:,:,jj,mm).*n(:,:,kk).*n(:,:,ii)); 
                end
            end
        end
    end
end

function [rhsR]=calcRhsDilatationalRij(Gij,Rdhat,Rsdhat,Bdhat,k,K,a,n_particles,nshells)
    rhsR=zeros(n_particles,nshells,3,3); 
     for ii=1:3
        for jj=1:3
            rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)+a^2.*(k(:,:,ii).*K(:,:).*Bdhat(:,:,jj)+k(:,:,jj).*K(:,:).*Bdhat(:,:,ii));
            for kk=1:3
                rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)-Gij(ii,kk).*Rdhat(:,:,kk,jj)-Gij(jj,kk).*Rdhat(:,:,kk,ii);
                for mm=1:3
                    rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)-2*Gij(kk,mm).*(Rsdhat(:,:,mm,ii).*k(:,:,kk).*k(:,:,jj)+Rsdhat(:,:,mm,jj).*k(:,:,kk).*k(:,:,ii)); 
                end
            end
        end
    end
end

function [rhsR]=calcRhsCrossRij(Gij,Rshat,Rsdhat,Bshat,k,K,a,n_particles,nshells)
    rhsR=zeros(n_particles,nshells,3,3); 
     for ii=1:3
        for jj=1:3
            rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)+a^2.*k(:,:,jj).*K(:,:).*Bshat(:,:,ii);
            for kk=1:3
                rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)-Gij(ii,kk).*Rsdhat(:,:,kk,jj)-Gij(jj,kk).*Rsdhat(:,:,ii,kk);
                for mm=1:3
                    rhsR(:,:,ii,jj)=rhsR(:,:,ii,jj)+2*Gij(kk,mm)*(-Rshat(:,:,mm,ii).*k(:,:,kk).*k(:,:,jj)+Rsdhat(:,:,mm,jj).*k(:,:,kk).*k(:,:,ii)); 
                end
            end
        end
    end
end

function [rhsB]=calcRhsBs(Gij,Rshat,Rsdhat,Bshat,k,K,n_particles,nshells)
    rhsB=zeros(n_particles,nshells,3);
    for ii=1:3
        for kk=1:3
            rhsB(:,:,ii)=rhsB(:,:,ii)-K.*k(:,:,kk).*Rshat(:,:,kk,ii)-K.*k(:,:,kk).*Rsdhat(:,:,ii,kk)-Gij(ii,kk).*Bshat(:,:,kk);
            for mm=1:3
                rhsB(:,:,ii)=rhsB(:,:,ii)+2*Gij(kk,mm).*Bshat(:,:,mm).*k(:,:,kk).*k(:,:,ii);
            end
        end
    end
end

function [rhsB]=calcRhsBd(Gij,Rdhat,Rsdhat,Bshat,Bdhat,Phihat,k,K,a,n_particles,nshells)
    rhsB=zeros(n_particles,nshells,3);
    for ii=1:3
    for kk=1:3
        rhsB(:,:,ii)=rhsB(:,:,ii)-K.*k(:,:,kk).*Rsdhat(:,:,kk,ii)-K.*k(:,:,kk).*Rdhat(:,:,kk,ii)-Gij(ii,kk).*Bdhat(:,:,kk);
        for mm=1:3
            rhsB(:,:,ii)=rhsB(:,:,ii)-2*Gij(kk,mm).*Bshat(:,:,mm).*k(:,:,kk).*k(:,:,ii);
        end
    end
    rhsB(:,:,ii)=rhsB(:,:,ii)+K.*k(:,:,ii).*a^2.*Phihat;
    end
end


%% initialization and conversion functions

function [Rhat]=initClusters(nshells,n_particles,k,tke_init,spectrum)
    Rhat=zeros(n_particles,3,3);
    mag_k=logspace(-3,5,nshells)';
    if strcmp(spectrum,'VK')
        % Von Karman Spectrum
        Cnuk=0.4843;
        E=2.0*tke_init*Cnuk*(mag_k.^4)./((1.0 + mag_k.^2).^(17/6));
    elseif strcmp(spectrum,'SC97')
        % Simone Cambon 1997 Spectrum
        kp=8;
        A = 0.00025907;
        E=A*tke_init*exp(-2*mag_k.^2/(kp^2)).*mag_k.^4;
    else
        fprintf('spectrum specified not available')
        return;
    end
    Eint = trapz(mag_k,E);
    for pp=1:n_particles
      Rhat(pp,1,1) = Eint*(1-k(pp,1)^2);
      Rhat(pp,2,2) = Eint*(1-k(pp,2)^2);
      Rhat(pp,3,3) = Eint*(1-k(pp,3)^2);
      Rhat(pp,1,2) = Eint*(-k(pp,1)*k(pp,2));Rhat(pp,2,1)=Rhat(pp,1,2);
      Rhat(pp,1,3) = Eint*(-k(pp,1)*k(pp,3));Rhat(pp,3,1)=Rhat(pp,1,3);
      Rhat(pp,2,3) = Eint*(-k(pp,2)*k(pp,3));Rhat(pp,3,2)=Rhat(pp,2,3);

    end
end

function [Rhat,Rshat,Rdhat,Rsdhat,Bhat,Bshat,Bdhat,Phihat,k,mag_k] = initVelSpecTen(n_particles,nshells,max_k_exp,tke_init,tke_d_init,a0,gamma,spectrum,intType)
    Rhat=zeros(n_particles,nshells,3,3);
    Rshat=zeros(n_particles,nshells,3,3);
    Rdhat=zeros(n_particles,nshells,3,3);
    Rsdhat=zeros(n_particles,nshells,3,3);
    Bhat=zeros(n_particles,nshells,3);
    Bshat=zeros(n_particles,nshells,3);
    Bdhat=zeros(n_particles,nshells,3);
    Phihat=zeros(n_particles,nshells,1);
    k=zeros(n_particles,3);
    tke_s_init=tke_init-tke_d_init;
    % generate wave-vectors randomly
    pp=1;
    while (pp<n_particles+1)
        u1 = rand*2.0 - 1.0;
        u2 = rand*2.0 - 1.0;
        if (u1*u1 + u2*u2 < 1.0)
	        k(pp,1) = 2.0*u1*sqrt(1.0-u1*u1-u2*u2);
	        k(pp,2) = 2.0*u2*sqrt(1.0-u1*u1-u2*u2);
	        k(pp,3) = (1.0-2.0*(u1*u1+u2*u2));
	        pp = pp+1;
        end
    end
    if strcmp(intType,'quad')
        GL_pts_wts=readmatrix('GL_pts_wts.dat');
        X=GL_pts_wts(:,1)';
        mag_k=repmat(X,n_particles,1);
    else
        mag_k=repmat(logspace(-2,max_k_exp,nshells),n_particles,1);
    % mag_k=repmat(linspace(10^-3,10^max_k_exp,nshells),n_particles,1);
    end
    if strcmp(spectrum,'VK')
        % Von Karman Spectrum
        Cnuk=0.4843;
        E=2.0*tke_init*Cnuk*(mag_k.^4)./((1.0 + mag_k.^2).^(17/6));
        E_d=2.0*tke_init_d*Cnuk*(mag_k.^4)./((1.0 + mag_k.^2).^(17/6));
    elseif strcmp(spectrum,'SC97')
        % Simone Cambon 1997 Spectrum
        kp=8;
        A = 0.00025907;
        E=A*tke_init*exp(-2*mag_k.^2/(kp^2)).*mag_k.^4;
        E_d=A*tke_d_init*exp(-2*mag_k.^2/(kp^2)).*mag_k.^4;
    else
        fprintf('spectrum specified not available')
        return;
    end
    for ss=1:nshells
        for pp=1:n_particles
          Rhat(pp,ss,1,1) = E(pp,ss)*(1-(k(pp,1)^2))/(4*pi*mag_k(pp,ss)^2);
          Rhat(pp,ss,2,2) = E(pp,ss)*(1-(k(pp,2)^2))/(4*pi*mag_k(pp,ss)^2);
          Rhat(pp,ss,3,3) = E(pp,ss)*(1-(k(pp,3)^2))/(4*pi*mag_k(pp,ss)^2);
          Rhat(pp,ss,1,2) = E(pp,ss)*(-(k(pp,1)*k(pp,2)))/(4*pi*mag_k(pp,ss)^2);Rhat(pp,ss,2,1)=Rhat(pp,ss,1,2);
          Rhat(pp,ss,1,3) = E(pp,ss)*(-(k(pp,1)*k(pp,3)))/(4*pi*mag_k(pp,ss)^2);Rhat(pp,ss,3,1)=Rhat(pp,ss,1,3);
          Rhat(pp,ss,2,3) = E(pp,ss)*(-(k(pp,2)*k(pp,3)))/(4*pi*mag_k(pp,ss)^2);Rhat(pp,ss,3,2)=Rhat(pp,ss,2,3);
          Phihat(pp,ss)= -E_d(pp,ss)/(2*pi*mag_k(pp,ss)^2*a0^2);%*(1-(k(pp,3)^2))/(2*pi*mag_k(pp,ss)^2*gamma*a0^2)
          % Bhat(pp,ss,1)= sqrt(Rhat(pp,ss,1,1))*sqrt(abs(Phihat(pp,ss)));
          % Bhat(pp,ss,2)= sqrt(Rhat(pp,ss,2,2))*sqrt(abs(Phihat(pp,ss)));
          % Bhat(pp,ss,3)= sqrt(Rhat(pp,ss,3,3))*sqrt(abs(Phihat(pp,ss)));
        end
    end
    Rshat=Rhat.*(tke_s_init/tke_init);
    Rdhat=Rhat.*(tke_d_init/tke_init);
    k = repmat(reshape(k,n_particles,1,3),1,nshells,1);
end



%% OLD CODE
% function [Rhat,k,mag_k] = initVelSpecTen(n_particles,max_k_exp,k,tke_init)
%     Rhat=zeros(n_particles,3,3);
%     Cnuk=0.4843;
%     % mag_k=logspace(-2,max_k_exp,n_particles)';
%     mag_k=linspace(10^-2,10^max_k_exp,n_particles)';
%     % Von karman Spectrum
%     E=2.0*tke_init*Cnuk*(mag_k.^4)./((1.0 + mag_k.^2).^(17/6));
%     % % Simone Cambon 1997 Spectrum
%     % kp=8;
%     % A = 0.00025907;
%     % E=A*tke_init*exp(-2*mag_k.^2/(kp^2)).*mag_k.^4;
%     for pp=1:n_particles
%       ss=pp;
%       Rhat(pp,1,1) = E(ss)*(1-(k(pp,1)^2))/(4*pi*mag_k(ss)^2);
%       Rhat(pp,2,2) = E(ss)*(1-(k(pp,2)^2))/(4*pi*mag_k(ss)^2);
%       Rhat(pp,3,3) = E(ss)*(1-(k(pp,3)^2))/(4*pi*mag_k(ss)^2);
%       Rhat(pp,1,2) = E(ss)*(-(k(pp,1)*k(pp,2)))/(4*pi*mag_k(ss)^2);Rhat(pp,2,1)=Rhat(pp,1,2);
%       Rhat(pp,1,3) = E(ss)*(-(k(pp,1)*k(pp,3)))/(4*pi*mag_k(ss)^2);Rhat(pp,3,1)=Rhat(pp,1,3);
%       Rhat(pp,2,3) = E(ss)*(-(k(pp,2)*k(pp,3)))/(4*pi*mag_k(ss)^2);Rhat(pp,3,2)=Rhat(pp,2,3);
%       k(pp,:)=k(pp,:)*mag_k(ss);
%     end
% end
% 
% function [rhsR,rhsk, rhsB, rhsPhi] = calcRHSTerms(Gij,Rhat,k,Bhat,Phihat,a0,n_particles)
% rhsR=zeros(n_particles,3,3);
% rhsk=zeros(n_particles,3);
% rhsB=zeros(n_particles,3);
% 
% rhsPhi=calcRhsPhi(Gij,Bhat,Phihat,k,n_particles);
%      for ii=1:3
%         [rhsk(:,ii), rhsB(:,ii)]=calcRhsBi(ii,Gij,Rhat,Bhat,Phihat,k,a0,n_particles);
%         for jj=1:3
%             rhsR(:,ii,jj)=calcRhsRij(ii,jj,Gij,Rhat,Bhat,k,a0,n_particles);
%         end
%      end
%      % symmetrize
%     rhsR(:,2,1)=rhsR(:,1,2);
%     rhsR(:,3,1)=rhsR(:,1,3);
%     rhsR(:,3,2)=rhsR(:,2,3);
% end
% 
% 
% function rhsR=calcRhsRij(ii,jj,Gij,Rhat,Bhat,k,a,n_particles)
%     rhsR=zeros(n_particles,1); 
%     for kk=1:3
%         rhsR=rhsR-Gij(ii,kk)*Rhat(:,kk,jj)-Gij(jj,kk)*Rhat(:,kk,ii)+2*Gij(kk,kk)*Rhat(:,ii,jj);
%     end
%     rhsR=rhsR+k(:,ii).*Bhat(:,jj)*a^2+k(:,jj).*Bhat(:,ii)*a^2;
%     % rhsR=rhsR-mag_k.*k(:,ii).*Bhat(:,jj)*a^2-mag_k.*k(:,jj).*Bhat(:,ii)*a^2;
%     % %incompressible
%     % n= k./sqrt(sum(k.^2,2));
%     % % n=k;
%     % for kk=1:3
%     %     for mm=1:3
%     %         rhsR=rhsR+2*Gij(kk,mm)*(Rhat(:,ii,mm).*n(:,kk).*n(:,jj)+Rhat(:,jj,mm).*n(:,kk).*n(:,ii)); 
%     %     end
%     % end
% end
% 
% 
% function [rhsPhi]=calcRhsPhi(Gij,Bhat,Phihat,k,n_particles)
%     rhsPhi=zeros(n_particles,1); 
%     % rhsK2=zeros(n_particles,1);
%     for kk=1:3
%         rhsPhi=rhsPhi+2*k(:,kk).*Bhat(:,kk)+2*Gij(kk,kk)*Phihat;
%         % rhsPhi=rhsPhi+2*mag_k.*k(:,kk).*Bhat(:,kk)+2*Gij(kk,kk)*Phihat;
%         % for nn=1:3
%         %     rhsK2 = rhsK2-2*Gij(kk,nn)*k(:,kk).*k(:,nn);
%         % end
%     end
% 
% end
% 
% function [rhsk, rhsB]=calcRhsBi(ii,Gij,Rhat,Bhat,Phihat,k,a,n_particles)
%     rhsB=zeros(n_particles,1);
%     rhsk=zeros(n_particles,1);
%     for kk=1:3
%         rhsB=rhsB+k(:,kk).*Rhat(:,kk,ii)-Gij(ii,kk)*Bhat(:,kk)+2*Gij(kk,kk)*Bhat(:,ii);
%         % rhsB=rhsB+mag_k.*k(:,kk).*Rhat(:,kk,ii)-Gij(ii,kk)*Bhat(:,kk)+2*Gij(kk,kk)*Bhat(:,ii);
%         rhsk=rhsk-Gij(kk,ii)*k(:,kk);
%         % for mm=1:3
%         %     rhsk=rhsk+Gij(kk,mm)*k(:,kk).*k(:,mm).*k(:,ii);
%         % end
%     end
%     rhsB=rhsB-k(:,ii).*Phihat*a^2;
%     % rhsB=rhsB-mag_k.*k(:,ii).*Phihat*a^2;
% end