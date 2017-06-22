function [J,IX]  = maximalMarginalDiversity(Data,Label)

J = zeros(1,size(Data,2));
nb_classe = length(unique(Label));
PC = zeros(1,nb_classe);
for k = 1:nb_classe
    
    eval(['indices_' num2str(k) '= find(Label==' num2str(k) ');'])
    eval(['classe_' num2str(k) '=Data(indices_' num2str(k) ',:);'])
    eval(['PC(' num2str(k) ')= length(classe_' num2str(k) ');'])
    
end
taille = max(PC);
PC = PC/length(Label);


for p = 1:size(Data,2)
    
    [~,b] = hist(Data(:,p),ceil(sqrt(taille)));clear dummy
    prob_MCK = zeros(1,ceil(sqrt(taille)));
    for k = 1:nb_classe
        eval(['a = hist(classe_' num2str(k) '(:,' num2str(p) '),b);'])
        eval(['deltaC = (max(classe_' num2str(k) '(:,' num2str(p) '))-min(classe_' num2str(k) '(:,' num2str(p) ')))/length(b);'])
        eval(['alpha = 1/(length(classe_' num2str(k) '(:,' num2str(p) '))*deltaC);'])
        eval(['prob_' num2str(k) '= a*alpha;'])
        eval(['prob_MCK= prob_MCK+prob_' num2str(k) ';'])
    end
    prob_MCK = prob_MCK/nb_classe;
    
    for k = 1:nb_classe
    eval(['temp(' num2str(k) ') = PC(' num2str(k) ')*Distance(prob_' num2str(k) ',prob_MCK);'])
    end
    
    J(p) = sum(temp);
end

[~,IX] = sort(J,'descend');clear dummy