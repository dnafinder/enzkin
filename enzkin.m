function enzkinout=enzkin(S,v)
% ENZyme KINetics is the study of the chemical reactions that are catalysed by
% enzymes. In enzyme kinetics the reaction rate is measured and the effects of
% varying the conditions of the reaction investigated. Studying an enzyme
% kinetics in this way can reveal the catalytic mechanism of this enzyme, its
% role in metabolism, how its activity is controlled, and how a drug or a poison
% might inhibit the enzyme.
% Michaelis–Menten kinetics approximately describes the kinetics of many
% enzymes. It is named after Leonor Michaelis and Maud Menten. This kinetic
% model is relevant to situations where very simple kinetics can be assumed,
% (i.e. there is no intermediate or product inhibition, and there is no
% allostericity or cooperativity).
% The Michaelis–Menten equation relates the initial reaction rate v0  to the
% substrate concentration [S]. The corresponding graph is a rectangular
% hyperbolic function; the maximum rate is described as Vmax (asymptote); the
% concentration of substrate where the v0 is the half of Vmax is the
% Michaelis-Menten costant (Km).
% To determine the maximum rate of an enzyme mediated reaction, a series of
% experiments is carried out with varying substrate concentration ([S]) and the
% initial rate of product formation is measured. 'Initial' here is taken to mean
% that the reaction rate is measured after a relatively short time period,
% during which complex builds up but the substrate concentration remains
% approximately constant and the quasi-steady-state assumption will hold.
% Accurate values for Km and Vmax can only be determined by non-linear
% regression of Michaelis-Menten data.
% The Michaelis-Menten equation can be linearized using several techniques.
% ENZKIN uses 6 regression models (2 non-linear and 4 linear) to obtain the
% kinetic parameters.
%
% Syntax: 	enzkinout=enzkin(S,v)
%      
%     Inputs:
%           S - data array of substrate concentrations
%           v - data array of measured initial velocity
%     Outputs:
%           - Vmax and Km estimation by:
%                ° Michaelis-Menten non linear regression
%                ° loglog non linear regression
%                ° Lineweaver-Burk linear regression
%                ° Hanes-Woolf linear regression
%                ° Eadie-Hofstee linear regression
%                ° Scatchard linear regression
%           - for the linear regressions, all regression data are summarized
%           - Plots
%
% The function requires another function of mine MYREGR. If it is not present on
% the computer, enzkin will try to download it from FEX
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2010). Enzkin: a tool to estimate Michaelis-Menten kinetic
% parameters
% http://www.mathworks.com/matlabcentral/fileexchange/26653

%Input Error handling
p=inputParser;
addRequired(p,'S',@(x) validateattributes(x,{'numeric'},{'vector','real','finite','nonnan','nonnegative','row','increasing'}));
addRequired(p,'v',@(x) validateattributes(x,{'numeric'},{'vector','real','finite','nonnan','nonnegative','row','increasing'}));
parse(p,S,v)
assert(length(S)==length(v))
clear p

if exist('myregr.m','file')==0
    filename=unzip('https://it.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/15473/versions/10/download/zip','prova');
    Index = contains(filename,'myregr.m');
    current=cd;
    copyfile(filename{Index},current)
    rmdir('prova','s')
    clear filename Index current 
end

%set the costants
n=length(S); %number of data
vc=tinv(0.95,n-2); %critical value for confidence intervals
tr=repmat('-',1,100); %divisor
txtlbl={'Michaelis & Menten non linear fit' '[S]' 'v'; 
    sprintf('Lineweaver-Burk\n(x=1/S; y=1/v)') '1/[S]' '1/V'; ...
    sprintf('Hanes-Woolf\n(x=s; y=s/v)') '[S]' '[S]/v'; ...
    sprintf('Logarithmic non linear fit') 'log(S)' 'log(v)'; ...
    sprintf('Eadie-Hofstee\n(x=v/s; y=v)') 'v/[S]' 'v'; ...
    sprintf('Scatchard\n(x=v; y=v/s)') 'v' 'v/[S]'}; %labels
KM=NaN(6,4); VMAX=KM;

disp(tr)
fprintf('Lineweaver-Burk linearization (x=1/S; y=1/v => Vmax=1/q; Km=m/q)...\n')
disp(tr)
[x,Idx]=sort(1./S); y=1./v(Idx);
[slope, intercept, stat]=myregr(x,y,0);
VMAX(1,1:2)=[1/intercept.value intercept.se/intercept.value^2]; VMAX(1,3:4)=VMAX(1,1)+[-1 1].*vc*VMAX(1,2);
KM(1,1:2)=[slope.value/intercept.value realsqrt((intercept.value*slope.se)^2+(slope.value*intercept.se)^2)/intercept.value^2]; KM(1,3:4)=KM(1,1)+[-1 1].*vc*KM(1,2);
disppar(1)
clear x y slope intercept stat 
disp(' ')

disp(tr)
fprintf('Hanes-Woolf linearization (x=S; y=S/v => Vmax=1/m; Km=q/m)...\n')
disp(tr)
[slope, intercept, stat]=myregr(S,S./v,0);
VMAX(2,1:2)=[1/slope.value slope.se/slope.value^2]; VMAX(2,3:4)=VMAX(2,1)+[-1 1].*vc*VMAX(2,2);
KM(2,1:2)=[intercept.value/slope.value realsqrt((slope.value*intercept.se)^2+(intercept.value*slope.se)^2)/slope.value^2]; KM(2,3:4)=KM(2,1)+[-1 1].*vc*KM(2,2);
disppar(2)
clear slope intercept stat 
disp(' ')

disp(tr)
[x,Idx]=sort(v./S); y=v(Idx);
fprintf('Eadie-Hofstee linearization (x=v/S; y=v => Vmax=q; Km=-m)...\n')
disp(tr)
[slope, intercept, stat]=myregr(x,y,0);
VMAX(3,:)=[intercept.value intercept.se intercept.lv intercept.uv]; 
KM(3,1:2)=[-slope.value slope.se]; KM(3,3:4)=KM(3,1)+[-1 1].*vc*KM(3,2);
disppar(3)
clear x y slope intercept stat 
disp(' ')

disp(tr)
[x,Idx]=sort(v); y=v(Idx)./S(Idx);
fprintf('Scatchard linearization (x=v; y=v/s => Vmax=-q/m; Km=-1/m)...\n')
disp(tr)
[slope, intercept, stat]=myregr(x,y,0);
KM(4,1:2)=[-1/slope.value slope.se/slope.value^2]; KM(4,3:4)=KM(4,1)+[-1 1].*vc*KM(4,2);
VMAX(4,1:2)=[-intercept.value/slope.value realsqrt((-intercept.se/slope.value)^2+(-intercept.value/slope.value^2*slope.se)^2)]; VMAX(4,3:4)=VMAX(4,1)+[-1 1].*vc*VMAX(4,2);
disppar(4)
clear x y slope intercept stat 
disp(' ')

%Michaelis and Menten non linear fit
%log-log non linear fit
xfit = S(:); yfit = v(:); lyfit = log(v(:));
%check if x and y are ok
ok = isfinite(log(xfit)) & isfinite(lyfit);
if ~all(ok)
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
%set the start point: 
%The same of Michaelis & Menten fit but using logs
%fitting option:
fo = fitoptions('method','NonlinearLeastSquares','Robust','On','Lower',[0 0]);
set(fo,'Startpoint',log([mean(KM(1:4,1)) mean(VMAX(1:4,1))]));
%Set the log(Michaelis & Menten) equation
ft = fittype('log(Vmax)+log(x)-log(Km+x)',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'Km', 'Vmax'});
% Fit this model using data
[cfl, goodness] = fit(xfit(ok),lyfit(ok),ft,fo);
%Display the results
disp(tr)
fprintf('Logarithm nonlinear fit...\n')
disp(tr)
disp(cfl)
fprintf('\tR = \t%0.4f\n',realsqrt(goodness.rsquare))
disp(tr)
disp(' ')
p=coeffvalues(cfl); KM(5,1)=p(1); VMAX(5,1)=p(2); 
p=confint(cfl); KM(5,3:4)=p(1,:)'; VMAX(5,3:4)=p(2,:)';
KM(5,2)=(KM(5,4)-KM(5,1))/vc;
VMAX(5,2)=(VMAX(5,4)-VMAX(5,1))/vc;
clear p st fo ft cf goodness

%Hyperbolic fit
%check if x and y are ok
ok = isfinite(xfit) & isfinite(yfit);
if ~all(ok)
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
%set fitting options: 
fo = fitoptions('method','NonlinearLeastSquares','Robust','On','Lower',[min(KM(1:5,1)) min(VMAX(1:5,1))],'Upper',[max(KM(1:5,1)) max(VMAX(1:5,1))],'DiffMaxChange',9.9999999999999995475e-07);
set(fo,'Startpoint',[mean(KM(1:5,1)) mean(VMAX(1:5,1))]);
%Set the Michaelis & Menten equation
ft = fittype('(Vmax*x)/(Km+x)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'Km', 'Vmax'});
% Fit this model using the data
[cf, goodness] = fit(xfit(ok),yfit(ok),ft,fo);
%Display the results
disp(tr)
fprintf('Michaelis & Menten nonlinear fit...\n')
disp(tr)
disp(cf)
fprintf('\tR = \t%0.4f\n',sqrt(goodness.rsquare))
disp(tr)
disp(' ')
p=coeffvalues(cf); KM(6,1)=p(1); VMAX(6,1)=p(2); 
p=confint(cf); KM(6,3:4)=p(1,:)'; VMAX(6,3:4)=p(2,:)';
KM(6,2)=(KM(6,4)-KM(6,1))/vc;
VMAX(6,2)=(VMAX(6,4)-VMAX(6,1))/vc;
clear p st fo ft cf goodness

disp(array2table(KM,'VariableNames',{'KM','Std_Err','Lower_bound','Upper_bound'},...
    'Rownames',{'Lineweaver_Burk','Hanes-Woolf','Eadie-Hofstee','Scatchard','Log_fit','Hyperbolic_fit'}));

disp(array2table(VMAX,'VariableNames',{'VMAX','Std_Err','Lower_bound','Upper_bound'},...
    'Rownames',{'Lineweaver_Burk','Hanes-Woolf','Eadie-Hofstee','Scatchard','Log_fit','Hyperbolic_fit'}));

if nargout
    enzkinout.KM=KM;
    enzkinout.VMAX=VMAX;
end

figure('Color',[1 1 1],'outerposition',get(groot,'ScreenSize'));

%Lineweaver-Burk plot 
m=KM(1,1)/VMAX(1,1); q=1/VMAX(1,1);
x=[-1/KM(1,1) 0 1./S];
y=[0 q 1./v];
stringa={['Vmax = ' num2str(VMAX(1,1)) '    Km = ' num2str(KM(1,1))];  ...
    'y = (Km/Vmax) * x + (1/Vmax)'; ...
    ['x=0; y=1/Vmax= ' num2str(y(2))]; ...
    ['y=0; x=-1/Km= ' num2str(x(1))]}; %labels
kingraph(2,m,q,stringa)

%Hanes-Woolf plot 
m=1/VMAX(2,1); q=KM(2,1)/VMAX(2,1);
x=[-KM(2,1) 0 S];
y=[0 q S./v];
stringa={['Vmax = ' num2str(VMAX(2,1)) '    Km = ' num2str(KM(2,1))];  ...
    'y = (1/Vmax) * x + (Km/Vmax)'; ...
    ['x=0; y=Km/Vmax=' num2str(y(2))]; ...
    ['y=0; x=-Km=' num2str(x(1))]}; %labels
kingraph(3,m,q,stringa)

%Eadie-Hofstee plot
m=-KM(3,1); q=VMAX(3,1);
x=[VMAX(3,1)/KM(3,1) 0 v./S];
y=[0 VMAX(3,1) v];
stringa={['Vmax = ' num2str(VMAX(3,1)) '    Km = ' num2str(KM(3,1))];  ...
    'y = -Km * x + Vmax'; ...
    ['x=0; y=Vmax=' num2str(y(2))]; ...
    ['y=0; x=Vmax/Km=' num2str(x(1))]}; %labels
kingraph(5,m,q,stringa)

%Scatchard plot
m=-1/KM(4,1); q=VMAX(4,1)/KM(4,1);
x=[VMAX(4,1) 0 v];
y=[0 VMAX(4,1)/KM(4,1) v./S];
stringa={['Vmax = ' num2str(VMAX(4,1)) '    Km = ' num2str(KM(4,1))]; ...
    'y = -1/Km * x + Vmax/Km'; ...
    ['x=0; y=Vmax/Km=' num2str(y(2))]; ...
    ['y=0; x=Vmax=' num2str(x(1))]}; %labels
kingraph(6,m,q,stringa)

%Michaelis & Menten plot 
x=S; y=v;
[~,xtick,ytick]=kingraph(1,VMAX(6,1),KM(6,1));
H=text(5*xtick,VMAX(6,1)+3*ytick,['Vmax = ' num2str(VMAX(6,1))]); set(H,'Color','m')
H=text(KM(6,1)+5*xtick,VMAX(6,1)/2,['Km = ' num2str(KM(6,1))]); set(H,'Color','m')

%log(Michaelis & Menten) plot 
[as,xtick,ytick]=kingraph(4,VMAX(5,1),KM(5,1));
H=text(as(1)+5*xtick,log(VMAX(5,1))+2*ytick,'log(v)=log(Vmax)+log(S)-log(Km+S)'); set(H,'Color','m')
H=text(as(1)+5*xtick,log(VMAX(5,1))+5*ytick,['log(Vmax) = ' num2str(log(VMAX(5,1))) '; Vmax = ' num2str(VMAX(5,1))]); set(H,'Color','m')
H=text(log(KM(5,1))+5*xtick,log(VMAX(5,1)/2),['log(Km) = ' num2str(log(KM(5,1))) '; Km = ' num2str(KM(5,1))]); set(H,'Color','m')

function disppar(I) %display parameters
    z=[struct2array(slope),(1-tcdf(abs(slope.value/slope.se),n-2))*2; ...
        struct2array(intercept),(1-tcdf(abs(intercept.value/intercept.se),n-2))*2; ...
        stat.r(1:4),(1-tcdf(abs(stat.r(1)/stat.r(2)),n-2))*2; ...
        VMAX(I,:),(1-tcdf(abs(VMAX(I,1)/VMAX(I,2)),n-2))*2; ...
        KM(I,:),(1-tcdf(abs(KM(I,1)/KM(I,2)),n-2))*2; ...
        ];
disp(array2table(z,'Rownames',{'Slope' 'Intercept' 'R' 'Km' 'Vmax'},...
    'VariableNames',{'Value','SE','Lower_bound','Upper_bound','p_value'}))
disp(tr)
end

function [out1,out2,out3]=kingraph(I,m,q,stringa) % main plot function
if I==4
    xtick=range(log(x))/50; ytick=range(log(y))/50;
else
    xtick=range(x)/50; ytick=range(y)/50; 
end
subplot(2,3,I)
hold on
if I==1 || I==4
    start=1;
else
    start=3;
end
%plot the data points (blue circles)
if I==4
    plot(log(x(start:end)),log(y(start:end)),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4,'LineStyle','none')
else
    plot(x(start:end),y(start:end),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4,'LineStyle','none')
end
if I~=1 && I~=4
   %plot the x and y intersections (red circles)
   plot(x(1:2),y(1:2),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4,'LineStyle','none')
end

%plot the regression line (orange line)
if I==1
    xs=linspace(0,max(x),500); ys=(m.*xs)./(q+xs);
    plot(xs,ys,'LineStyle','-','Color',[0.87 0.49 0])
elseif I==4
    xs=linspace(min(x),max(x),500); ys=log(m)+log(xs)-log(q+xs);
    plot(log(xs),ys,'LineStyle','-','Color',[0.87 0.49 0])
else
    xs=linspace(min(x),max(x),500);
    plot(xs,polyval([m q],xs),'LineStyle','-','Color',[0.87 0.49 0])
end

as=axis;
if nargout
    out1=as;
    out2=xtick;
    out3=ytick;
end
if I==1
    plot(as(1:2),[m m],'k--') %plot the Vmax asymptote
    plot([0 q],[m m]./2,'k--',[q q],[0 m/2],'k--') %plot the Km coordinates
elseif I==4
    as(4)=as(4)+10*ytick;
    axis(as)
    plot(as(1:2),log([m m]),'k--') %plot the Vmax asymptote
    ax=axis;
    plot([as(1) log(q)],log([m m]./2),'k--',log([q q]),[log(m/2) ax(3)],'k--') %plot the Km coordinates
else    
    %Mark the axis (black line')
    axis([as(1)-10*xtick as(2)+10*xtick as(3)-10*ytick as(4)+10*ytick])
    xk1=[0 0]; yk1=as(3:4)+[-10 10].*ytick;
    xk2=as(1:2)+[-10 10].*xtick; yk2=[0 0];
    plot(xk1,yk1,'k',xk2,yk2,'k')
end
hold off
axis square
%set the labels
title(txtlbl(I,1))
xlabel(txtlbl(I,2))
ylabel(txtlbl(I,3))
if  I~=1 && I~=4
    H=text(as(1)-5*xtick,as(4)+5*ytick,stringa(1)); set(H,'Color','m')
    H=text(as(1)-5*xtick,as(4)-ytick,stringa(2)); set(H,'Color','m')
    H=text(5*xtick,y(2),stringa(3)); set(H,'Color','m')
    if I==2 || I==3
        H=text(x(1)+0.5*xtick,-5*ytick,stringa(4)); 
    elseif I==5 || I==6
        H=text(x(1)-12*xtick,-5*ytick,stringa(4));
    end
    set(H,'Color','m')
end
end

end
