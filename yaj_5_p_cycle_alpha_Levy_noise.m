clear; clc;

a=0.55;
q=8;
b=1;
rr=1;
m=1;


alpha = 1.5;
beta = -0.6;
D = 0.1; 

T = 20000; N = 50*2^10; dt = T/N;
alpha_1 = D*(dt)^(1/alpha);
dL = zeros(1,N);  
L = zeros(1,N);
dL(1) = sqrt(dt)*randn; 
L(1) = dL(1); 

X_0 = 0.55; 
Xtemps=X_0;
Xstroges=zeros(1,N);Xstroges(1)=Xtemps;
Xstroges_1 = [];
time_first = [];
tep = [];


for j = 2:N
    
    if alpha ~= 1 
        V = unifrnd(-pi/2,pi/2);
        W = exprnd(1);
        const = beta * tan(pi*alpha/2);
        B = atan( const );
        S = (1 + const * const).^(1/(2*alpha));
        dL(j) = ( S * sin( alpha*V + B ) ./ ( cos(V) ).^(1/alpha) .* ...
            ( cos( (1-alpha) * V - B ) ./ W ).^((1-alpha)/alpha) )*alpha_1;
    else                             % General case, alpha = 1
        V = unifrnd(-pi/2,pi/2);
        W = exprnd(1);
        piover2 = pi/2;
        sclshftV =  piover2 + beta * V ;
        dL(j) = ( 1/piover2 * ( sclshftV .* tan(V) - ...
            beta * log( (piover2 * W .* cos(V) ) ./ sclshftV ) ) )*alpha_1;
    end
    
    L(j) = L(j-1) + dL(j);
    
    Xtemps = Xtemps + (a - b*Xtemps + rr*Xtemps^q/(m^q + Xtemps^q))*dt + dL(j); 
    Xstroges(j)=Xtemps;
    
    
    if Xtemps >= 0 & Xtemps <=3
        Xstroges_1 = [Xstroges_1 Xtemps];
    end
    
    
    if abs(Xtemps-1.5)<=0.01
        tep = [tep j];
        time_first = [time_first dt*(j)];   
    end
end

% subplot(2,1,1);
plot([0:dt:(T-dt)],Xstroges,'r')
xlabel('t','FontSize',16)
ylabel('X(t)','FontSize',16,'Rotation',0)
axis([0 T -10 10]); 


% subplot(2,1,2);
[f,xi] = ksdensity(Xstroges_1); 
plot(xi,f,'b--o');
xlabel('P Concentration','FontSize',13);
ylabel('P_s','FontSize',13);
xlim([0,2]);







