function [sources] = sources_tracking(lambda,X,amps)

% source tracking to plot the regularization path for SFW

nblamb = length(lambda);
sources = cell(0,1);


init = 1;
while sum(X{init})==0
    init=init+1;
end


for l = init:nblamb
    A=pdist2(X{l},X{l-1});
    

    temp=zeros(length(X{l}(:,1)),1);
    for a = 1:length(X{l-1}(:,1))
        [~, n] = min(A(:));
        [x, y] = ind2sub(size(A),n);
        temp(y)=1;
        sources{x}=[sources{x};[lambda(l),X{l}(y,:),amps{l}(y,:)]];

        A(x,:)=Inf;
        A(:,y)=Inf;
    end

    for i = 1:length(temp)
        if temp(i)==0
            sources{length(sources)+1}=[lambda(l),X{l}(i,:),amps{l}(i,:)];
        end
    end
end


end

