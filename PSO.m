function [Lmresult Diresult PPresult Massresult] = PSO(iternumber) 

Low_bound = [0.0000001 1 1];
Upp_bound = [10 100 200];
 
var = 3; % # of variables
particles = 100; % population size
w_max = 0.9; % inertia weight
w_min = 0.4;
c_1 = 2 ; % acceleration factor
c_2 = 2;

max_iter = iternumber;
max_run = 15;

for run = 1:max_run

    for i=1:particles
        for j=1:var
            x_0(i,j) = round(Low_bound(j) + rand()*(Upp_bound(j)-Low_bound(j)));
        end
    end
    x=x_0; %initial population
    V=0.1*x_0; %initial velocity
    for i=1:particles
        ftemp(i,1) = RF_mass_NoStrc(x_0(i,:));
    end
    [fmin0, index0] = min(ftemp);
    pBEST=x_0;
    gBEST=x_0(index0,:);
    
    iter=1;
    tolerance=1;
    
    while iter<=max_iter && tolerance>10^-11
        w=w_max-(w_max-w_min)*iter/max_iter; %inertial weight update
        
        for i = 1:particles
            for j=1:var
                
                V(i,j) = w*V(i,j) + c_1*rand()*(pBEST(i,j) - x(i,j))+ c_2*rand()*(gBEST(1,j)-x(i,j)); %velocity update
            end
        end
        
        x = x + V; %position update
        % check for boundaries
        for i=1:particles
            for j=1:var
                if x(i,j)<Low_bound(j)
                    x(i,j) = Low_bound(j);
                elseif x(i,j)>Upp_bound(j)
                    x(i,j)=Upp_bound(j);
                end
            end
        end
        
        % evaluate cost function
        for i=1:particles
            f(i,1)=RF_mass_NoStrc(x(i,:));
        end 
        % update pbest and function
        for i=1:particles
            if f(i,1)<ftemp(i,1)
                pBEST(i,:)=x(i,:);
                ftemp(i,1)=f(i,1);
            end
        end
        
        %find the best particle and store function value and corresponding
        %iteration
        [minmin,index] = min(ftemp);
        ffmin(iter,run) = minmin;
        ffiter(run) = iter;
        
        % update gbest
        if minmin<fmin0
            gBEST = pBEST(index,:);
            fmin0=minmin;
        end
        
        %calculate tolerance after a number of iterations
        if iter>100;
            tolerance=abs(ffmin(iter-100,run)-fmin0);
        end
        
        % displaying iterative results
        if iter==1
            disp(sprintf('iteration   best particle   cost func'));
            
        end
        disp(sprintf('%8g   %8g      %8.4f',iter,index,fmin0));
        iter=iter + 1;
    end
    fvalue= RF_mass_NoStrc(gBEST);
    fff(run)=fvalue;
    rgbest(run,:)=gBEST;
    disp(sprintf(' --------------------------------'));
end
disp(sprintf('final results\n'));
[bestfun,bestrun]=min(fff)
best_vars=rgbest(bestrun,:)


Lmresult = rgbest(bestrun,1);
Diresult = rgbest(bestrun,2);
PPresult = rgbest(bestrun,3);
Massresult = bestfun;
%ploting
plot(ffmin(1:ffiter(bestrun),bestrun),'-k');

    
      
         
        
    
    
   