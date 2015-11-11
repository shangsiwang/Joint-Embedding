function [H, lambda, objective,iter] = oned_embed(A,maxiter,group)

%     A=randn(20,20,200);
%     A(:,:,1:100)=randn(20,20,100)>0+0;
%     A(:,:,101:200)=randn(20,20,100)>1+0;
%     A(1:10,1:10,1:100)=randn(10,10,100)>-0.5+0;
    
    [m,~,n]=size(A);
    
    Aver=zeros(m);
        for i=1:n
        Aver=A(:,:,i)+Aver;
        end
    Aver=Aver+Aver';
    [V, ~] = eigs(Aver, 1, 'LA');

    H=V;
    %H=randn(m,1);
    
    H=H/norm(H);
    
    lambda=zeros(n,1);
    

    if nargin==1
        maxiter=20;
    end
    
    
    if nargin<=2
        group=1:n;
    end
    ugroup=unique(group);
    ugroup=reshape(ugroup,1,length(ugroup));
    
    iter=0;
    while iter<maxiter
    iter=iter+1;
    h=reshape(H*H',m^2,1);
    
    for g=ugroup
        hth=0;
        hty=0;
        for i=1:n
            if(group(i)==g)
            hth=hth+h'*h;
            v1=reshape(A(:,:,i),m^2,1);
            hty=hty+h'*v1;
            end
        end
        lambda(group==g)=1/hth*hty;
    end
    
    
    objective=cmp_obj(A,H,lambda);
   
    
    %disp(objective)
    
    Gradiant=zeros(m,1);
    for i=1:n
        Gradiant=Gradiant+4*lambda(i)*(lambda(i)*(H*H')-A(:,:,i))*H;
    end
    Gradiant=Gradiant/n;
    %disp(Gradiant)
    
    if iter==1
        GradiantP=Gradiant;
    end
    
    if iter>1 && norm(Gradiant)<norm(GradiantP)/(10^5)
        %disp(1);
        break;
    end
    
    stepsize=1;
    alpha=0.1;
    beta=0.1;
   
    Ht=H-Gradiant*stepsize;
    Ht=Ht/norm(H);
    
    tempobjective=cmp_obj(A,Ht,lambda);
    
    while tempobjective-objective>= -alpha*stepsize*(Gradiant'*Gradiant)
        stepsize=stepsize*beta;
        if stepsize<=2.2204e-016
            break;
        end
        
        Ht=H-Gradiant*stepsize;
        Ht=Ht/norm(Ht);
    tempobjective=cmp_obj(A,Ht,lambda);
    end
    
    
    
        if stepsize<=2.2204e-016 || abs(tempobjective-objective)<(1+abs(objective))/10^5
%             disp('2')
%             disp(stepsize)
%             disp(tempobjective)
%             disp(objective)
            break;
        end
    
    H=Ht;
    
    
    end
    

end


function obj = cmp_obj(A,H,lambda)
 [~,~,n]=size(A);
    tempobjective=0;
    for i=1:n
        tempobjective=tempobjective+norm(A(:,:,i)-H*H'*lambda(i))^2;
    end
    obj=tempobjective;
    
end