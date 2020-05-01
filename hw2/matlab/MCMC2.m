clear;
clc;

load data ;


%% package
addpath('../demo/jplv7/distrib')
addpath('../demo/jplv7/gibbs')
addpath('../demo/jplv7/util')

%%
% N=30;
G=50;
T=5000;    % number of discarded samples during Markov process
R=2;     % number of replication for Monte Carlo experiment
k = 13;    % number of covariates

lambda_R=zeros(1,R);
beta_R=zeros(2*k,R);
sige_R=zeros(1,R);

%%
parpool(2);
tic;
parfor r=1:R
    disp('r'); disp(r);
    
%     tic;
        
    % assign parameter in prior distributions %
    
    beta_0=zeros([2*k, 1]);
    B_0=eye(2*k)*3;
    ALPHA_0=1;
    rho_0=2.2;
    eta_0=0.1;
    
    acc_1=0; acc_2=0; acc_3=0; acc_4=zeros(G,1);
    
    acc_rate1=zeros(T,1);
    acc_rate2=zeros(T,1);
    acc_rate3=zeros(T,1);
    acc_rate4=zeros(G,T);
    
    lambda_T =zeros([1,T]);  % save for lambda
    beta_T   =zeros([2*k,T]);  % save for beta
    alpha_T  =zeros([G,T]);  % save for alpha
    sige_T   =zeros([1,T]);  % save for Sigma_e^2
    
    % starting value of draw
    sige_T(1)=1;
        
    for t=2:T
        
        accept=0;
        %while accept==0;                                        % propose lambda^*
            if t<2
                lambda_1=mvnrnd(lambda_T(t-1),eye(1)*0.1^2);
            else
                lambda_1=mvnrnd(lambda_T(t-1),cov(lambda_T(1:t-1))*2.38^2)*0.95+mvnrnd(lambda_T(t-1),eye(1)*0.1^2)*0.05;
            end
            %if lambda_1>-1/(min(max(max(sum(w(:,:,:,r),1))),max(max(sum(w(:,:,:,r),2))))) && lambda_1<1/(min(max(max(sum(w(:,:,:,r),1))),max(max(sum(w(:,:,:,r),2)))));
            %    accept=1;
            %end;
        %end;
        
        pp_l = 1;
        
        
        for g=1:G
            
            Ng = length(W{g});
            V = sige_T(t-1)*eye(Ng);
            
            S_1 = eye(Ng) - lambda_1 * W{g};
            S_2 = eye(Ng) - lambda_T(t-1) * W{g};
            
            ep_1 = S_1*Y{g} - [X{g}, W{g}*X{g}] * beta_T(:,t-1) - ones([Ng,1])*alpha_T(g,t-1);
            ep_2 = S_2*Y{g} - [X{g}, W{g}*X{g}] * beta_T(:,t-1) - ones([Ng,1])*alpha_T(g,t-1);
            
            like_1 = det(S_1)*exp(-(1/2)*(ep_1)'/V*(ep_1));
            like_2 = det(S_2)*exp(-(1/2)*(ep_2)'/V*(ep_2));
            
            pp_l=pp_l*(like_1/like_2);
        end
        
        pp_l=min(pp_l,1);
        
        if rand(1) <= pp_l
            
            lambda_T(t) = lambda_1;
            acc_2 = acc_2+1;
        else
            
            lambda_T(t) = lambda_T(t-1);
        end
        
        acc_rate2(t,1) = acc_2/t;
        
        
        % THE SAMPLING OF BETA FROM PROSTERIOR DISTRIBUTION %
        ZVY=0;
        ZVX=0;
        
        for g=1:G
            
            Ng = size(W{g}, 1); % row number
            
            V = sige_T(t-1)*eye(Ng);
            
            SS = eye(Ng) - lambda_T(t)*W{g};
            YY = SS*Y{g} - ones([Ng,1])*alpha_T(g,t-1);
            ZZ = [X{g} W{g}*X{g}];
            
            ZVX = ZVX + ZZ'/V*ZZ;
            ZVY = ZVY + ZZ'/V*YY;
        end
        beta_T(:,t) = norm_rnd(inv(inv(B_0)+ZVX)) + ((inv(B_0)+ZVX)\((B_0)\beta_0+ZVY));
        
        
        % THE SAMPLING OF SIGMA_E^2 FROM PROSTERIOR DISTRIBUTION %
        
%         ep_v=zeros(2020, 1); % known total observation 2020
        
        ee = 0;
        total_n = 0;
        
        for g=1:G
            
            Ng = length(W{g});
            total_n = total_n + Ng ;
            
            
            SS = eye(Ng) - lambda_T(t)*W{g};
            ep = SS*Y{g} - [X{g} W{g}*X{g}]*beta_T(:,t) - ones([Ng,1])*alpha_T(g,t-1);
            
%             ep_v((g-1)*N+1:g*N)=ep;

            ee = ee + ep' * ep ;
            
            
        end
        
        rho_1 = rho_0 + total_n ;
        sige_T(t) = (ee + eta_0)/chis_rnd(1,rho_1);
        
        % THE SAMPLING OF ALPHA_G FROM PROSTERIOR DISTRIBUTION %
        
        for g=1:G
            
            Ng = length(W{g}) ;
            
            dd = (ALPHA_0^(-1)+(sige_T(t))^(-1)*ones([1,Ng])/(eye(Ng))*ones([Ng,1]))^(-1);
            
            SS = eye(Ng) - lambda_T(t)*W{g};
            YY = SS*Y{g};
            XX = [X{g} W{g}*X{g}];
            
            alpha_T(g,t) = sige_T(t)^(-1)*dd*ones([1,Ng])/(eye(Ng))*(YY-XX*beta_T(:,t))+randn(1)*sqrt(dd);
            
        end
                
%         if (t/100)-round(t/100)==0
%             disp(t);
%             disp(lambda_T(t));
%             disp(beta_T(:,t)');
%             disp(sige_T(t));
%             
%             subplot(3,1,1)
%             plot(lambda_T(1:t));figure(gcf);
%             title('\lambda')
%             subplot(3,1,2)
%             plot(beta_T(1,1:t));figure(gcf);
%             title('\beta_1')
%             subplot(3,1,3)
%             plot(beta_T(2,1:t));figure(gcf);
%             title('\beta_2')
%             drawnow;
%         end
                
    end
%     time=toc;
    
%     disp('lambda');disp(mean(lambda_T(1000:t)));
%     disp('beta');disp(mean(beta_T(:,1000:t),2));
%     disp('sigma_e^2');disp(mean(sige_T(1000:t)));
%     disp('time'); disp(time);
%     
%     lambda_R(r)=mean(lambda_T(1000:t));
%     beta_R(:,r)=mean(beta_T(:,1000:t),2);
%     sige_R(r)=mean(sige_T(1000:t));
%     
%     info.q = 0.025;
%     info.r = 0.005;
%     info.s = 0.95;
%     info.p1 = 0.3;
%     info.p2 = 0.5;
%     vnames = strvcat('lambda');
%     coda(lambda_T(1000:5:t)',vnames,info);
%     plot(lambda_T(1000:5:t));
%     
%     pause;
     
end
toc;

disp('lambda');disp(mean(lambda_R,2));
disp(std(lambda_R,0,2));
disp('beta');disp(mean(beta_R,2));
disp(std(beta_R,0,2));
disp('sigma_e^2');disp(mean(sige_R,2));
disp(std(sige_R,0,2));
