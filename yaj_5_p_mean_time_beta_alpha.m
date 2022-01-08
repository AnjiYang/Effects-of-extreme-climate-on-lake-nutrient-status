
clear; clc;

a=0.55;
q=8;
b=1;
rr=1;
m=1;


alpha = 1.7;
beta = -1:0.1:1;
D=0.1; 
mean_time=[];

for ii = 1:length(beta)
mean_time_total = 0;
for i = 1:1000
T = 2000; N = 50*2^10; dt = T/N;
alpha_1 = D*(dt)^(1/alpha);
dL = zeros(1,N);  
L = zeros(1,N);
dL(1) = sqrt(dt)*randn; 
L(1) = dL(1); 

X_0=0.55; 
Xtemps=X_0;
Xstroges=zeros(1,N);Xstroges(1)=Xtemps;
Xstroges_1 = [];
time_first = [];
tep = [];
f_funtion = @(x) a - b*x + rr*x^q/(m^q + x^q); 

for j = 2:N
    
    if alpha ~= 1 
        V = unifrnd(-pi/2,pi/2);
        W = exprnd(1);
        const = beta(ii) * tan(pi*alpha/2);
        B = atan( const );
        S = (1 + const * const).^(1/(2*alpha));
        dL(j) = ( S * sin( alpha*V + B ) ./ ( cos(V) ).^(1/alpha) .* ...
            ( cos( (1-alpha) * V - B ) ./ W ).^((1-alpha)/alpha) )*alpha_1;
    else                             % General case, alpha = 1
        V = unifrnd(-pi/2,pi/2);
        W = exprnd(1);
        piover2 = pi/2;
        sclshftV =  piover2 + beta(ii) * V ;
        dL(j) = ( 1/piover2 * ( sclshftV .* tan(V) - ...
            beta(ii) * log( (piover2 * W .* cos(V) ) ./ sclshftV ) ) )*alpha_1;
    end
    
    L(j) = L(j-1) + dL(j);
    
    k1 = dt*f_funtion(Xtemps);
    k2 = dt*f_funtion(Xtemps + 1/2*(k1));
    k3 = dt*f_funtion(Xtemps + 1/2*(k2));
    k4 = dt*f_funtion(Xtemps + k3);
 
    Xtemps = Xtemps + (k1 + 2*k2 + 2*k3 + k4)/6 + dL(j); 
    

    Xstroges(j)=Xtemps;
    
    
    if Xtemps >= 0 && Xtemps <=2
        Xstroges_1 = [Xstroges_1 Xtemps];
    end
    
    
    if Xstroges(j) >= 1.5
        tep = [tep j];
        time_first = [time_first dt*(j)];
    end
    
end

mean_time_total = mean_time_total + time_first(1);

end
mean_time = [mean_time log10(mean_time_total/i)]; 

end


plot(beta,mean_time,'-*b'); 
xlabel('\beta','FontSize',13)
ylabel('log10(MFPT)','FontSize',13)








