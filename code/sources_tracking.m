function [sources] = sources_tracking(lambda,X,amps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inititalisation

nblamb = length(lambda);
sources = cell(0,1);

%on initialise le dictionnaire sources avec les premières sources non
%nulles suivant lambda décroissant
init = 1;
while sum(X{init})==0
    init=init+1;
end

% 
% for i = 1:length(X{init}(:,1))
%     sources{i}=[lambda(init),X{init}(i,:),amps{init}(i,:)];
% end

for l = init:nblamb
    %A is the array containing the distances between each pair of point
    A=pdist2(X{l},X{l-1});
    
    %temp permet de retrouver les nouvelles sources qui n'ont pas été
    %apairées. La valeur de temp(i) est assoicée à la source i, elle vaut 0
    %si la source n'a pas été aprairée et 1 sinon
    temp=zeros(length(X{l}(:,1)),1);
    for a = 1:length(X{l-1}(:,1))
        %We take the minimal distance and we link the two points
        [~, n] = min(A(:));
        [x, y] = ind2sub(size(A),n);
        %On note que la source y a été associée à une ancienne source
        temp(y)=1;
        sources{x}=[sources{x};[lambda(l),X{l}(y,:),amps{l}(y,:)]];
        %On retire artificiellement les sources x et y de la liste des
        %sources à apairer
        A(x,:)=Inf;
        A(:,y)=Inf;
    end
    %Pour toutes les nouvelles sources non apairées, on crée une cellule
    %dans le dictionnaire sources
    for i = 1:length(temp)
        if temp(i)==0
            sources{length(sources)+1}=[lambda(l),X{l}(i,:),amps{l}(i,:)];
        end
    end
end


end

