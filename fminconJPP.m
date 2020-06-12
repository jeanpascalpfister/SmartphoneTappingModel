function [x,f,flag,xlist] = fminconJPP(fun,x0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
% JPP 31.8.2017


n = length(x0);
converged = false;
x = x0;
i = 0;
imax = OPTIONS.MaxIter;
tolx = OPTIONS.TolX;     
tolfun = OPTIONS.TolFun;
dis = OPTIONS.Display;
eta = 1;
xlist(:,1) = x0;
Dx = 0;
df = 10; % any initial value larger than tolfun

while ~converged
    [f,g] = fun(x);
    if strcmp(dis,'iter')
        disp(['i = ' num2str(i) ', dx = ' num2str(norm(Dx)/sqrt(n)) ', C = ' num2str(f) ', theta = ' num2str(x')])        
    end
           
    dx = -eta*g;
    xnew = x+dx;
            
    % enforce the bounds
    xnew = xnew.*((LB<=xnew).*(xnew<= UB)) + LB.*(xnew<LB) + UB.*(UB<xnew);
    
    % implement constraints (inequalities)
    if (isempty(A) && isempty(B)) % no inequality constraints
        Dx = xnew-x;
        xCOK = xnew;
    else
        if all(A*xnew<=B) % all constraints are OK
            Dx = xnew-x;
            xCOK = xnew;
        else        
            [ma, indma] = max(A*x-B); % find the worst constraint violation;
            if strcmp(dis,'iter')
                disp(['constraint violation: returning in the feasible region. ma = ' num2str(ma)])
            end
            w = A(indma,:);
            xB = xnew -(w*xnew-1)*w'/(w*w');    
            Dx = xB-x;
            xCOK = xB;
        end        
    end
    
    % checks if the new solution is better than the previous one
    fCOK = fun(xCOK);
    if fCOK < f
        df = f-fCOK;
        x = xCOK;
       
    else
        if strcmp(dis,'iter')
            disp('decrease learning rate')
        end
        eta = eta/2;
    end
   
    i = i+1;
    if nargout > 3
        xlist(:,i+1) = x;
    end
    if i==imax 
       converged = true;
       disp('JPP: max number of iterations exceeded')
        if nargout > 2
            flag = 0;
        end
    end
    if norm(Dx)/sqrt(n) < tolx
       converged = true;
       disp('JPP: |dx|<tolx')
       if nargout > 2
            flag = 1;
       end
    end
    if df < tolfun
       converged = true;
       disp('JPP: |dL|<tolfun')
       if nargout > 2
            flag = 2;
       end
    end

end

end
