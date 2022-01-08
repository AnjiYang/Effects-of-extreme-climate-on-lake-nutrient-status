clear; clc; close;

a=0.55;
q=8;
b=1;
rr=1;
m=1;

P_mean_temp = [];
total_mean = [];
for D = 0.01:0.01:0.1
    for J = 1:10
        
        alpha = 1.5; 
        beta = -0.6;
        
        T = 2000; N = 50*2^10; dt = T/N;
        alpha_1 = (D*dt)^(1/alpha);
        dL = zeros(1,N);  
        L = zeros(1,N);
        dL(1) = sqrt(dt)*randn; 
        L(1) = dL(1); 
        
        X_0=0.55; 
        Xtemps=X_0;
        Xstroges=zeros(1,N);Xstroges(1)=Xtemps;
        Xstroges_1 = [];
        
        f_funtion = @(x) a - b*x + rr*x^q/(m^q + x^q); 
        
        for j = 2:N
            
            if alpha ~= 1 % General case, alpha not 1 i.e. X ~ S_alpha(1,beta,0) 
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
            
            k1 = dt*f_funtion(Xtemps);
            k2 = dt*f_funtion(Xtemps + 1/2*(k1));
            k3 = dt*f_funtion(Xtemps + 1/2*(k2));
            k4 = dt*f_funtion(Xtemps + k3);
            
            Xtemps = Xtemps + (k1 + 2*k2 + 2*k3 + k4)/6 + dL(j);
            Xstroges(j)=Xtemps;
            
            if Xtemps >= 0 & Xtemps <=2 
                Xstroges_1 = [Xstroges_1 Xtemps];
            end
        end
        
        
        [f,xi] = ksdensity(Xstroges_1);
        f = f*(xi(2) - xi(1));
        P_mean = sum(xi.*f);
        
        P_mean_temp = [P_mean_temp P_mean];
    end
    total_mean_temp = mean(P_mean_temp);
    total_mean_temp
    total_mean = [total_mean total_mean_temp];
    
end


D = 0.01:0.01:0.1;
% load('P_mean_vs_noise_indensity.mat','total_mean')
plot(D,total_mean,'m','LineWidth',2,'Marker','p','MarkerSize',10);
xlabel('Noise intensity \sigma','FontSize',13);
ylabel('The average P concentration','FontSize',13);
xlim([0.01,0.1]);


