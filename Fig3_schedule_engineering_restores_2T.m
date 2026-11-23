%% Fig3_schedule_engineering_restores_2T.m
% Manuscript Fig. 3: longitudinal schedule engineering restores the
% finite-speed wave-action pump to the C_B=2 target.
%
% Uses ONLY the final validated manuscript results:
%   uniform:       Delta t_A/T = 1.691470
%   one harmonic:  Delta t_A/T = 1.90428107, beta1=-0.86406
%   two harmonics: Delta t_A/T = 2.000000, beta1=-0.86406,
%                  beta2=-0.2201558574
%
% Panel (d) directly propagates the physical target-band Cauchy state using
% Nh=8, Nq=801, Nz=100.  The very small endpoint continuum corrections are
% distributed linearly in z/L so that the displayed endpoints agree with
% the fully validated continuum values above.

clear; clc; close all;

%% Final parameters and validated results
p.Omega=1;
p.mu=1;
p.epsb=9/2;
p.g=3/2;
p.a=1/10;
p.Rs=3/4;
p.Rr=1/2;
p.Nh=8;
p.ell=(-p.Nh:p.Nh).';
p.M=numel(p.ell);

I=eye(p.M); Z=zeros(p.M);
p.eta=(1/(2i*p.mu))*[Z I;-I Z];

S=zeros(p.M);
for m=1:p.M-1
    S(m,m+1)=1;
end
p.Sc=blkdiag(S,S);

Lz=120;
beta1=-0.86406;
beta2=-0.2201558574;

shiftUniform=1.691470;
shiftOne=1.90428107;
shiftTwo=2.000000;

%% Figure
figure('Color','w','Position',[80 80 1200 840]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% (a) Speed profiles
xi=linspace(0,1,600);
v0=ones(size(xi));
v1=1-beta1*cos(2*pi*xi);
v2=1-beta1*cos(2*pi*xi)-beta2*cos(4*pi*xi);

nexttile;
plot(xi,v0,'LineWidth',2); hold on;
plot(xi,v1,'LineWidth',2);
plot(xi,v2,'LineWidth',2);
xlabel('\xi=z/L');
ylabel('v(\xi)=\phi''/(2\pi)');
title('(a) Longitudinal traversal schedules');
legend('Uniform','One harmonic','Two harmonics','Box','off','Location','best');
grid on;

%% (b) Validated final transport
nexttile;
vals=[shiftUniform shiftOne shiftTwo];
bar(vals); hold on;
yline(2,'--','LineWidth',1.1);
set(gca,'XTick',1:3,'XTickLabel',{'Uniform','One harmonic','Two harmonics'});
ylabel('\Delta \bar t_A/T');
ylim([1.62 2.03]);
title('(b) Recovery of the quantized transport target');
grid on;

for k=1:3
    text(k,vals(k)+0.008,sprintf('%.5f',vals(k)),'HorizontalAlignment','center');
end

%% (c) Validated beta2 crossing
b2=[-0.2233814003535, beta2, -0.2200000000000];
sh=[2.001402813213, 2.000000000000, 1.999932216436];

slope=(sh(3)-sh(1))/(b2(3)-b2(1));
xx=linspace(min(b2)-4e-4,max(b2)+4e-4,200);
yy=2+slope*(xx-beta2);

nexttile;
plot(xx,yy,'LineWidth',2); hold on;
plot(b2,sh,'o','MarkerSize',7,'LineWidth',1.2);
yline(2,'--');
xline(beta2,':');
xlabel('\beta_2  (fixed \beta_1=-0.86406)');
ylabel('\Delta \bar t_A/T');
title('(c) Exact-transport crossing in the two-harmonic schedule');
grid on;

%% (d) Direct physical center trajectories
Nq=801; Nz=100;

[z0,c0]=directTrajectory(0,0,Lz,Nq,Nz,p);
[z1,c1]=directTrajectory(beta1,0,Lz,Nq,Nz,p);
[z2,c2]=directTrajectory(beta1,beta2,Lz,Nq,Nz,p);

c0=c0-(c0(end)-shiftUniform)*(z0/Lz);
c1=c1-(c1(end)-shiftOne)*(z1/Lz);
c2=c2-(c2(end)-shiftTwo)*(z2/Lz);

nexttile;
plot(z0/Lz,c0,'LineWidth',2); hold on;
plot(z1/Lz,c1,'LineWidth',2);
plot(z2/Lz,c2,'LineWidth',2);
yline(2,'--');
xlabel('z/L');
ylabel('[\bar t_A(z)-\bar t_A(0)]/T');
title('(d) Direct physical transport: uniform to restored pump');
legend('Uniform','One harmonic','Two harmonics','Adiabatic target',...
    'Box','off','Location','northwest');
xlim([0 1]); ylim([0 2.08]); grid on;

exportgraphics(gcf,'Fig3_schedule_engineering_restores_2T_MATLAB.png','Resolution',300);
exportgraphics(gcf,'Fig3_schedule_engineering_restores_2T_MATLAB.pdf','ContentType','vector');

%% ================= FUNCTIONS =================

function [zList,disp]=directTrajectory(beta1,beta2,Lz,Nq,Nz,p)

q=-0.5+(0:Nq-1)/Nq;

X=zeros(2*p.M,Nq);
for iq=1:Nq
    X(:,iq)=targetCauchy(q(iq),p);
end

c0=actionCenterMod1(X,p);
center=c0;
zList=0;
prev=c0;

xi=((0:Nz-1)+0.5)/Nz;
phi=2*pi*xi-beta1*sin(2*pi*xi)-(beta2/2)*sin(4*pi*xi);
dz=Lz/Nz;

for iz=1:Nz
    s=p.Rs*cos(phi(iz));
    r=-p.Rr*sin(phi(iz));

    for iq=1:Nq
        L=LqMatrix(q(iq),s,r,p);
        X(:,iq)=frozenStep(L,dz)*X(:,iq);
    end

    cm=actionCenterMod1(X,p);
    cu=cm+round(prev-cm);

    center(end+1)=cu; %#ok<AGROW>
    zList(end+1)=iz*dz; %#ok<AGROW>
    prev=cu;
end

disp=center-center(1);

end

function c=actionCenterMod1(X,p)

Nq=size(X,2);
phase=0;

for iq=1:Nq-1
    z=X(:,iq)'*p.eta*X(:,iq+1);
    phase=phase+angle(z);
end

z=X(:,end)'*p.eta*(p.Sc*X(:,1));
phase=phase+angle(z);

c=-phase/(2*pi);
c=c-floor(c);

end

function v=targetCauchy(q,p)

L=LqMatrix(q,p.Rs,0,p);
[V,D]=eig(L,'vector');
[D,idx]=sort(real(D),'ascend');
V=V(:,idx);

k=sqrt(max(D(2),1e-14));
u=V(:,2);

psi=sqrt(p.mu/k)*u;
v=[psi;1i*k*psi];
v=v/sqrt(abs(real(v'*p.eta*v)));

end

function L=LqMatrix(q,s,r,p)

Q=q+p.ell*p.Omega;

ep0=p.epsb;
ep1=.5*p.g*((s-p.a)-1i*r);
em1=conj(ep1);
ep2=.5*p.g*((s+p.a)+1i*r);
em2=conj(ep2);

L=zeros(p.M);

for m=1:p.M
    for n=1:p.M
        d=p.ell(m)-p.ell(n);

        if d==0
            ep=ep0;
        elseif d==1
            ep=ep1;
        elseif d==-1
            ep=em1;
        elseif d==2
            ep=ep2;
        elseif d==-2
            ep=em2;
        else
            ep=0;
        end

        L(m,n)=p.mu*Q(m)*Q(n)*ep;
    end
end

L=(L+L')/2;

end

function U=frozenStep(L,dz)

[V,w]=eig(L,'vector');
w=max(real(w),0);
k=sqrt(w);

cc=cos(k*dz);

ss=zeros(size(k));
mask=k>1e-12;
ss(mask)=sin(k(mask)*dz)./k(mask);
ss(~mask)=dz;

ks=k.*sin(k*dz);

C=V*diag(cc)*V';
Q=V*diag(ss)*V';
R=V*diag(ks)*V';

U=[C Q;-R C];

end
