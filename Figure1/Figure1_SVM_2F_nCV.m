
X = [HC; MDD]; % HC and MDD are Subject x Number of functional networks arrays of network size
Y = [ones(size(HC,1),1)*1; ones(size(MDD,1),1)*2];

% Define the hyperparameter ranges
boxConstraintValues = [0.1 1 10 100]; 
kernelScaleValues = [0.1 1 10 100];  

% preallocate some variables;
AllPredictedLabels = zeros(100,length(Y));
FinalConfusionmatrix = [];
FinalWeights = [];

rng(1); % for reproducibility

% 100 iterations;
for nperm = 1:100
    
    % Create partition for
    % outer cross-validation
    outerCV = cvpartition(Y,'KFold',2);
    
    % Initialize array for outer fold accuracies
    outerFoldAccuracies = zeros(outerCV.NumTestSets, 1);
    
    % preallocate;
    OuterWeights = [];
    OuterConfusionMatrix = [];
    
    % this is the outer loop
    for i = 1:outerCV.NumTestSets
        
        % Separate data into training &
        % testing sets for the outer loop
        outerTrainData = X(outerCV.training(i), :);
        outerTrainLabels = Y(outerCV.training(i));
        outerTestData = X(outerCV.test(i), :);
        outerTestLabels = Y(outerCV.test(i));
        
        % over-sample the minority class in the training data;
        [outerTrainData,outerTrainLabels] = smote(outerTrainData,outerTrainLabels);
        
        % Create partition for inner cross-validation
        innerCV = cvpartition(outerTrainLabels,'KFold',2);
        
        % Initialize variables;
        bestBoxConstraint = NaN;
        bestKernelScale = NaN;
        bestAccuracy = 0;
        
        % Now, inner loop for hyperparameter tuning
        
        % sweep the box constraint values
        for boxConstraint = boxConstraintValues
            
            % sweep the kernel scale values
            for kernelScale = kernelScaleValues
                
                % preallocate accuracy values;
                acc = zeros(innerCV.NumTestSets,1);
                
                % sweep the test tests;
                for j = 1:innerCV.NumTestSets
                    
                    % Training and validation data for this inner fold
                    innerTrainData = outerTrainData(innerCV.training(j), :);
                    innerTrainLabels = outerTrainLabels(innerCV.training(j));
                    innerValidationData = outerTrainData(innerCV.test(j), :);
                    innerValidationLabels = outerTrainLabels(innerCV.test(j));
                    
                    % Train the SVM model
                    model = fitcsvm(innerTrainData, innerTrainLabels, ...
                        'BoxConstraint', boxConstraint, ...
                        'KernelScale', kernelScale);
                    
                    % Evaluate the model
                    predictions = predict(model,innerValidationData);
                    acc(j) = mean(predictions == innerValidationLabels);
                    
                end
                
                % Average accuracy for this
                % combination of hyperparameters
                avgAccuracy = mean(acc);
                
                % log best combination
                % of hyperparameters;
                if avgAccuracy > bestAccuracy
                    bestAccuracy = avgAccuracy;
                    bestBoxConstraint = boxConstraint;
                    bestKernelScale = kernelScale;
                end
                
            end
            
        end
        
        % Train the model on the entire outer training set w/ best hyperparameters
        finalModel = fitcsvm(outerTrainData, outerTrainLabels,'BoxConstraint',...
        bestBoxConstraint,'KernelScale', bestKernelScale);
        
        % Evaluate the model on the outer test set
        finalPredictions = predict(finalModel, outerTestData);
        outerFoldAccuracies(i) = mean(finalPredictions == outerTestLabels);
        
        % log the model weights
        % and confusion matrix;
        OuterWeights(i,:) = finalModel.Beta;
        OuterConfusionMatrix(:,:,i) = confmat(finalPredictions,outerTestLabels);
        AllPredictedLabels(nperm,outerCV.test(i)) = finalPredictions;
        
    end
    
    % Calculate the overall accuracy across all outer folds
    overallAccuracy(nperm) = mean(outerFoldAccuracies);
    
    % log the average model 
    % weights and confusion matrix
    FinalWeights(nperm,:) = mean(OuterWeights);
    FinalConfusionmatrix(:,:,nperm) = mean(OuterConfusionMatrix,3);
    
end
