%% Fig1_topological_pump_foundation.m
% Manuscript Fig. 1: intrinsic two-tone C_B=2 wave-action pump
clear; clc; close all;

Omega=1; mu=1; epsb=9/2; g=3/2; a=1/10; Rs=3/4; Rr=1/2;
Nh=8; ell=(-Nh:Nh).'; M=numel(ell);

figure('Color','w','Position',[80 80 1200 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% (a) pump loop
phi=linspace(0,2*pi,600);
s=Rs*cos(phi); r=-Rr*sin(phi);
nexttile;
plot(s,r,'LineWidth',2); hold on;
plot(a,0,'o','MarkerSize',7,'LineWidth',1.3);
plot(-a,0,'o','MarkerSize',7,'LineWidth',1.3);
axis equal; grid on;
xlabel('s'); ylabel('r');
title('(a) Two-tone C_B=2 pump loop');
text(a+0.035,0.035,'N_1: q/\Omega=1/2, s=a');
text(-a+0.035,-0.075,'N_2: q/\Omega=0, s=-a');

%% (b) instantaneous bands at phi=0
q=linspace(-0.5,0.5,401); K=zeros(numel(q),4);
for iq=1:numel(q)
    L=makeL(q(iq),0,ell,mu,epsb,g,a,Rs,Rr,Omega);
    ev=sort(real(eig(L)),'ascend');
    K(iq,:)=sqrt(max(ev(1:4),0));
end
nexttile;
plot(q/Omega,K(:,1)/Omega,'LineWidth',1.1); hold on;
plot(q/Omega,K(:,2)/Omega,'LineWidth',2.6);
plot(q/Omega,K(:,3)/Omega,'LineWidth',1.1);
plot(q/Omega,K(:,4)/Omega,'--','LineWidth',0.9);
grid on;
xlabel('q/\Omega'); ylabel('k/\Omega');
title('(b) Instantaneous longitudinal bands at \phi=0');
legend('lower band','target B','upper band','next band','Box','off','Location','north');

%% (c) protecting gaps
phiG=linspace(0,2*pi,121);
qG=linspace(-0.5,0.5,302); qG(end)=[];
gapLo=zeros(size(phiG)); gapHi=zeros(size(phiG));
for ip=1:numel(phiG)
    dl=zeros(size(qG)); du=zeros(size(qG));
    for iq=1:numel(qG)
        L=makeL(qG(iq),phiG(ip),ell,mu,epsb,g,a,Rs,Rr,Omega);
        ev=sort(real(eig(L)),'ascend');
        kk=sqrt(max(ev(1:3),0));
        dl(iq)=kk(2)-kk(1); du(iq)=kk(3)-kk(2);
    end
    gapLo(ip)=min(dl); gapHi(ip)=min(du);
end
nexttile;
plot(phiG/(2*pi),gapLo/Omega,'LineWidth',2); hold on;
plot(phiG/(2*pi),gapHi/Omega,'LineWidth',2);
grid on;
xlabel('\phi/2\pi'); ylabel('minimum gap / \Omega');
title('(c) Protecting gaps remain open over the cycle');
legend('\Delta k_-','\Delta k_+','Box','off','Location','best');

%% (d) Wilson-loop action polarization
NqW=401; qW=linspace(-0.5,0.5,NqW+1); qW(end)=[];
phiW=linspace(0,2*pi,101);
S=zeros(M);
for m=1:M-1, S(m,m+1)=1; end
Pmod=zeros(size(phiW));

for ip=1:numel(phiW)
    U=zeros(M,NqW);
    for iq=1:NqW
        L=makeL(qW(iq),phiW(ip),ell,mu,epsb,g,a,Rs,Rr,Omega);
        [V,D]=eig((L+L')/2,'vector');
        [~,idx]=sort(real(D),'ascend');
        U(:,iq)=V(:,idx(2)); % target band = second k^2 eigenvector
    end
    Wprod=1;
    for iq=1:NqW-1
        z=U(:,iq)'*U(:,iq+1);
        Wprod=Wprod*z/abs(z);
    end
    z=U(:,end)'*(S*U(:,1));
    Wprod=Wprod*z/abs(z);
    Pmod(ip)=-angle(Wprod)/(2*pi);
end
Pun=unwrap(2*pi*Pmod)/(2*pi);
Pun=Pun-Pun(1);

nexttile;
plot(phiW/(2*pi),Pun,'LineWidth',2.5); hold on;
plot([0 1],[0 2],'--','LineWidth',1);
grid on; xlim([0 1]);
xlabel('\phi/2\pi');
ylabel('[\bar t_A(\phi)-\bar t_A(0)]/T');
title('(d) Adiabatic Wannier-center flow: C_B=2');
legend('action polarization','two-cell reference','Box','off','Location','best');

exportgraphics(gcf,'Fig1_topological_pump_foundation_MATLAB.png','Resolution',300);
exportgraphics(gcf,'Fig1_topological_pump_foundation_MATLAB.pdf','ContentType','vector');

function L=makeL(q,phi,ell,mu,epsb,g,a,Rs,Rr,Omega)
    s=Rs*cos(phi); r=-Rr*sin(phi);
    ep1=g/2*((s-a)-1i*r); em1=conj(ep1);
    ep2=g/2*((s+a)+1i*r); em2=conj(ep2);
    M=numel(ell); E=zeros(M,M);
    for m=1:M
        for n=1:M
            d=ell(m)-ell(n);
            if d==0, E(m,n)=epsb;
            elseif d==1, E(m,n)=ep1;
            elseif d==-1, E(m,n)=em1;
            elseif d==2, E(m,n)=ep2;
            elseif d==-2, E(m,n)=em2;
            end
        end
    end
    Q=q+ell*Omega;
    L=mu*(Q*Q.').*E;
    L=(L+L')/2;
end
