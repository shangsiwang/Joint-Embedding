function [H, lambda, objective] = tensor_embed(A,d,maxiter,group)

%     A=randn(100,100,200);
%     A(:,:,1:100)=randn(100,100,100)>0+0;
%     A(:,:,101:200)=randn(100,100,100)>0+0;
%     A(1:50,1:50,1:100)=randn(50,50,100)>-0.5+0;
    
    


    [m,~,n]=size(A);
    H=zeros(m,d);
    lambda=zeros(n,d);
    if nargin<4
        group=1:n;
    end
    Atemp=A;
    for k=1:d
        
        if nargin<=2
            maxiter=20;
        end
        [H1,lambda1,objective1,iter1]=oned_embed(Atemp,maxiter,group);
        H(:,k)=H1;
        lambda(:,k)=lambda1;
        objective=objective1;
        for i=1:n
            Atemp(:,:,i)=Atemp(:,:,i)-lambda1(i)*(H1*H1');
        end

        disp(['Dimension ' num2str(k) ' is finished in iteration ' num2str(iter1) '. ' 'Objective is ' num2str(objective1) '.']);
        
    end
    
end