%% code for the Model of mgrB safA interactions
%clear
% low, high
m = [0,0.4,4,40,40]; % [Mg2+]/K_M, the low value is calibrated on Fig 2d DmgrB
% PhoP - small proteins binding constants - normalized
f = 1;
%Variable binding constants
% lb1 = 0.1*f; % mgrB high salt 0.100  
% lb0 = 0.325*f; % mgrB low salt  0.325
% la1 = 0.525*f; % safA high salt 0.525
% la0 = 0.24*f; % safA low salt 0.240

%Constant binding constants
% lb1 = 0.3*f; % mgrB high salt 
% lb0 = 0.3*f; % mgrB low salt  
% la1 = 0.3*f; % safA high salt 
% la0 = 0.3*f; % safA low salt 

%Only mgrB binding constant varies
lb1 = 0.1*f; % mgrB high salt 
lb0 = 0.325*f; % mgrB low salt  
la1 = 0.3*f; % safA high salt 
la0 = 0.3*f; % safA low salt 

% fraction of active Q = of phosphorylated PhoP
act = @(a,b,m) (1+a./la0 + m.*a./la1)./(1+a./la0+b./lb0 + m.*(1+a./la1+b./lb1));

% Growth rate parameters
% PhoP-P binding constants to 
ca = 0.4; % p-safA % arbitrarily taken to 0.4
cq = 1.2; % p-PhoPQ %1.2 
cb = cq/4; % p-mgrB % from Fig 4b of Miyashiro and Goulian (2007), assuming their b_PQ = 1, we need (cq/cb)^2 = 16 -> cq/cb=4 
% Basal level of PhoP/Q and mgrB
qo=0.1;
bo=0.01;

% growth/time params
tlag = 3.4; % lagtime duration in units of 1/exponential growth rate
dt = 0.01; % time step in units of 1/expo growth rate
Tsimu = (tlag+3.4)/dt; % total simulation duration in steps

% evgs on?
s = 1; %1   set to 0 for pH=7
S = s*[1,1,1,0.65,1]; % [1,1,1,1,1] for no effect of Mg on safA induction

% variables
a = zeros(Tsimu+1,length(m)); % normalized safA cc
q = zeros(Tsimu+1,length(m)); % normalized PhoP/Q cc 
b = zeros(Tsimu+1,length(m)); % normalized mgrB cc
% initial values 
a(1,:)=0.*ones(1,length(m));
q(1,:)=qo.*ones(1,length(m));
b(1,:)=bo.*ones(1,length(m));

aq = act(a,b,m); % fraction of phosphorylating PhoQ 
pp = q.*aq; % PhoP-P concentration (normalized)

od = zeros(Tsimu+1,1); % optical density
od(1)=0.01;

% Euler integration
for i =1:Tsimu
    
    h=(i*dt>=tlag); % are we still in the lag phase?
    % h=(i*dt/tlag).^20./(1+(i*dt/tlag).^20); % for testing a smoother transition
    
    g =  1;% to modulate growth rate - deprecated g = 1 bc normalization
    
    % compute increments
    dod = dt.*(h*od(i,:));
    % additive effects of PhoP-P and EvgA binding, PhoP-P does not shut
    % down p-safA induction completely. Factor 0.2 from Burton et al jmb (2010)
    da = dt.*( g*S.*( (1+0.2.*(pp(i,:)/ca).^2) ./ (1+(pp(i,:)/ca).^2 ) ) - h*a(i,:) ); 
    % p-PhoP/Q: basal expression of PhoP/Q + positive feedback by PhoP-P
    dq =  dt.*( g*(qo + (1-qo).*(pp(i,:)/cq).^2./(1+(pp(i,:)/cq).^2)) - h*q(i,:) );
    % p-mgrB: Induction by PhoP-P  
    db = dt.*( g*(pp(i,:)/cb).^2./(1+(pp(i,:)/cb).^2) - h*b(i,:) );
    
    % update cc
    a(i+1,:) = a(i,:)+da;
    q(i+1,:) = q(i,:)+dq;
    b(i+1,:) = b(i,:)+db;
    od(i+1,:) = od(i,:)+dod;
    
    % compute new activity and cc of PhoP-P
    aq(i+1,:) = act(a(i+1,:),b(i+1,:),m); %act(0,0,m);%
    pp(i+1,:) = q(i+1,:).*aq(i+1,:);
end

% plots all outputs
t0 = 1;
t = (t0:Tsimu)*dt/log(2); % t in units of doubling time

figure
subplot(2,4,1)
plot(t,a(t0:Tsimu,:))
legend(num2str(m'))
title('safA')
%set(gca,'YScale','log')

subplot(2,4,2)
plot(t,b(t0:Tsimu,:))
legend(num2str(m'))
title('mgrB')

subplot(2,4,3)
plot(t,pp(t0:Tsimu,:))
legend(num2str(m'))
title('PhoP-P')

subplot(2,4,4)
plot(t,q(t0:Tsimu,:))
legend(num2str(m'))
title('PhoP-tot')


subplot(2,4,5)
semilogx(od(t0:Tsimu,:),a(t0:Tsimu,:))
legend(num2str(m'))
title('safA')
%set(gca,'YScale','log')

subplot(2,4,6)
semilogx(od(t0:Tsimu,:),b(t0:Tsimu,:))
legend(num2str(m'))
title('mgrB')

subplot(2,4,7)
semilogx(od(t0:Tsimu,:),pp(t0:Tsimu,:))
legend(num2str(m'))
title('PhoP-P')
%hold on

subplot(2,4,8)
plot(od(t0:Tsimu,:),q(t0:Tsimu,:))
legend(num2str(m'))
title('PhoP-P tot')

%% figures formated for paper

figure
subplot(2,2,1)
plot(t,a(2:end,:))
%legend('low Mg','high Mg, EvgSA deg','high Mg')
legend('0','0.4','4','40 + EvgSA red','40')
xlabel 'time / t1/2'
ylabel 'Normalized concentration'
title('safA')
set(gca,'Ylim',[0 3.5],'PlotboxAspectRatio',[1 1 1],'Fontsize',14)

subplot(2,2,2)
plot(t,b(2:end,:))
%legend('low Mg','high Mg, EvgSA deg','high Mg')
legend('0','0.4','4','40 + EvgSA red','40')
xlabel 'time / t1/2'
ylabel 'Normalized concentration'
title('mgrB')
set(gca,'Ylim',[0 1.5],'PlotboxAspectRatio',[1 1 1],'Fontsize',14)

subplot(2,2,3)
semilogx(od,a)
%legend('low Mg','high Mg, EvgSA deg','high Mg')
legend('0','0.4','4','40 + EvgSA red','40')
xlabel 'OD'
ylabel 'Normalized concentration'
title('safA')
set(gca,'Ylim',[0 3.5],'PlotboxAspectRatio',[1 1 1],'Fontsize',14)

subplot(2,2,4)
semilogx(od,b)
%legend('low Mg','high Mg, EvgSA deg','high Mg')
legend('0','0.4','4','40 + EvgSA red','40')
xlabel 'OD'
ylabel 'Normalized concentration'
title('mgrB')
set(gca,'Ylim',[0 1.5],'PlotboxAspectRatio',[1 1 1],'Fontsize',14)

%% to save max mgrB in case of constant binding constant
b_CstB = max(b);
%% to save max mgrB in case of variable binding constant
b_var = max(b);
%% to save max mgrB in case of only mgrB binding constant varies
b_BvarOnly = max(b);
%% to save max mgrB in case pH=7
b_pH7 = max(b);

%% Effect of having variable binding constants
cats = categorical(["0";"0.4";"4";"40";"40, EvgA red."]);
maxbs = [b_var',b_CstB',b_BvarOnly',b_pH7'];
maxbs = maxbs([1:3,5,4],:);
figure
bar(cats,maxbs) 
legend('pH5, var Ki','pH5, fixed Ki','pH5, var mgrB Ki only','pH 7')
xlabel '[Mg2+] / KM'
ylabel 'MgrB max'
set(gca,'Ylim',[0 1.5],'PlotboxAspectRatio',[1.5 1 1],'Fontsize',14)
%% Effect of having variable binding constants, plot v2
cats = categorical(["0";"0.4";"4";"40"]);
maxbs = [b_CstB',b_var',b_BvarOnly',b_pH7'];
maxbs = maxbs([1:3,5,4],:);
figure
plot(cats(1:4),maxbs(1:4,:),'-o') %bar
hold on
plot(cats(4),maxbs(5,2),'o')
legend('pH5, fixed Ki','pH5, var Ki','pH5, var mgrB Ki only','pH 7','pH5, var Ki, red EvgA')
xlabel '[Mg2+] / KM'
ylabel '[MgrB] at peak'
set(gca,'Ylim',[0 1.5],'PlotboxAspectRatio',[1.5 1 1],'Fontsize',14)