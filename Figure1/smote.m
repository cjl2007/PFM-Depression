function [X,Y] = smote(X,Y)
% cjl

    % group counts;
    [gc,g] = groupcounts(Y);
    
    % define the classes;
    min_c = g(gc==min(gc)); % minority
    maj_c = g(gc~=min(gc)); % majority 

    % extract 
    % classes;
    X1 = X(Y==min_c,:);
    X2 = X(Y~=min_c,:);
    
    % preserve 
    % original 
    % subjects 
    X1_Orig = X1;

    % over-sample the minority group ;
    for ii = 1:abs((size(X1,1)-size(X2,1)))
        Idx = 1:size(X1_Orig,1)==randi(size(X1_Orig,1)); % random starting point
        d = corr(X1_Orig(Idx==1,:)',X1_Orig','type','Spearman'); % distance
        [~,v] = sort(d,'descend'); % sort from most similar, to least similar
        n = v(2:6); % five nearest neighbors
        N = n(randi(length(n))); % random neighbor
        X1(end+1,:) = ((X1_Orig(N,:)- X1_Orig(Idx,:))* rand) + X1_Orig(Idx,:); % simulated minority class observation;
    end
    
    X = [X1; X2]; % features & labels (post-smote)
    Y = [ones(size(X1,1),1) * min_c ; ones(size(X2,1),1) * maj_c];

end