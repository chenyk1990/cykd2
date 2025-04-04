function [Dictionary,output] =KSVDnoiseS(...
    Data,... % an nXN matrix that contins N signals (Y), each of dimension n.
    param)
% =========================================================================
%                          K-SVD algorithm
% =========================================================================
% The K-SVD algorithm finds a dictionary for linear representation of
% signals. Given a set of signals, it searches for the best dictionary that
% can sparsely represent each signal. Detailed discussion on the algorithm
% and possible applications can be found in "The K-SVD: An Algorithm for
% Designing of Overcomplete Dictionaries for Sparse Representation", written
% by M. Aharon, M. Elad, and A.M. Bruckstein and appeared in the IEEE Trans.
% On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006.
% =========================================================================
% INPUT ARGUMENTS:
% Data                         an nXN matrix that contins N signals (Y), each of dimension n.
% param                        structure that includes all required
%                                 parameters for the K-SVD execution.
%                                 Required fields are:
%    K, ...                    the number of dictionary elements to train
%    numIteration,...          number of iterations to perform.
%    errorFlag...              if =0, a fix number of coefficients is
%                                 used for representation of each signal. If so, param.m must be
%                                 specified as the number of representing atom. if =1, arbitrary number
%                                 of atoms represent each signal, until a specific representation error
%                                 is reached. If so, param.errorGoal must be specified as the allowed
%                                 error.
%    preserveDCAtom...         if =1 then the first atom in the dictionary
%                                 is set to be constant, and does not ever change. This
%                                 might be useful for working with natural
%                                 images (in this case, only param.K-1
%                                 atoms are trained).
%    (optional, see errorFlag) L,...                 % maximum coefficients to use in OMP coefficient calculations.
%    (optional, see errorFlag) errorGoal, ...        % allowed representation error in representing each signal.
%    InitializationMethod,...  mehtod to initialize the dictionary, can
%                                 be one of the following arguments:
%                                 * 'DataElements' (initialization by the signals themselves), or:
%                                 * 'GivenMatrix' (initialization by a given matrix param.initialDictionary).
%    (optional, see InitializationMethod) initialDictionary,...      % if the initialization method
%                                 is 'GivenMatrix', this is the matrix that will be used.
%    (optional) TrueDictionary, ...        % if specified, in each
%                                 iteration the difference between this dictionary and the trained one
%                                 is measured and displayed.
%    displayProgress, ...      if =1 progress information is displyed. If param.errorFlag==0,
%                                 the average repersentation error (RMSE) is displayed, while if
%                                 param.errorFlag==1, the average number of required coefficients for
%                                 representation of each signal is displayed.
% =========================================================================
% OUTPUT ARGUMENTS:
%  Dictionary                  The extracted dictionary of size nX(param.K).
%  output                      Struct that contains information about the current run. It may include the following fields:
%    CoefMatrix                  The final coefficients matrix (it should hold that Data equals approximately Dictionary*output.CoefMatrix.
%    ratio                       If the true dictionary was defined (in
%                                synthetic experiments), this parameter holds a vector of length
%                                param.numIteration that includes the detection ratios in each
%                                iteration).
%    totalerr                    The total representation error after each
%                                iteration (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 0)
%    numCoef                     A vector of length param.numIteration that
%                                include the average number of coefficients required for representation
%                                of each signal (in each iteration) (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 1)
%    time                        Average time taken for one iteration of
%                                the dictionary update algorithm.
% =========================================================================

if (~isfield(param,'displayProgress'))
    param.displayProgress = 0;
end
totalerr(1) = 99999;
if (isfield(param,'errorFlag')==0)
    param.errorFlag = 0;
end

if (isfield(param,'TrueDictionary'))
    displayErrorWithTrueDictionary = 1;
    ErrorBetweenDictionaries = zeros(param.numIteration+1,1);
    ratio = zeros(param.numIteration+1,1);
else
    displayErrorWithTrueDictionary = 0;
    ratio = 0;
end
if (param.preserveDCAtom>0)
    FixedDictionaryElement(1:size(Data,1),1) = 1/sqrt(size(Data,1));
else
    FixedDictionaryElement = [];
end
% coefficient calculation method is OMP with fixed number of coefficients

if (size(Data,2) < param.K)
    disp('Size of data is smaller than the dictionary size. Trivial solution...');
    Dictionary = Data(:,1:size(Data,2));
    return;
elseif (strcmp(param.InitializationMethod,'DataElements'))
    Dictionary(:,1:param.K-param.preserveDCAtom) = Data(:,1:param.K-param.preserveDCAtom);
elseif (strcmp(param.InitializationMethod,'GivenMatrix'))
    Dictionary = param.initialDictionary;
end
% reduce the components in Dictionary that are spanned by the fixed
% elements
if (param.preserveDCAtom)
    tmpMat = FixedDictionaryElement \ Dictionary;
    Dictionary = Dictionary - FixedDictionaryElement*tmpMat;
    Dictionary = Dictionary(:,sum(Dictionary.^2)>10^-2);
end
%normalize the dictionary.
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.

% the K-SVD algorithm starts here.
time = 0;
for iterNum = 1:param.numIteration
    
    D = [FixedDictionaryElement,Dictionary];
    G = D'*D; epsilon = sqrt(size(D,1))*param.errorGoal;
    
    % find the coefficients
    if (param.errorFlag==0)
        %%% Faster computation of OMP
        CoefMatrix = OMPa(D,Data,param.m,0,0);
    else
        %%% Faster computation of OMP
        CoefMatrix = OMPa(D,Data,param.m,inf,param.errorGoal^2);
        param.m = 1;
    end
    
    tic;
    replacedVectorCounter = 0;
    rPerm = randperm(size(Dictionary,2));
%     rPerm = 1:size(Dictionary,2);    
    for j = rPerm
        [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...
            [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...
            CoefMatrix, 0.1*(param.errorGoal/1.15)^2);
        Dictionary(:,j) = betterDictionaryElement;
        if (param.preserveDCAtom)
            tmpCoef = FixedDictionaryElement\betterDictionaryElement;
            Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;
            Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));
        end
        replacedVectorCounter = replacedVectorCounter+addedNewVector;
    end
    time = time +toc;
    if (iterNum>1 && param.displayProgress)
        if (param.errorFlag==0)
            output.totalerr(iterNum-1) = sqrt(sum(sum((Data-[FixedDictionaryElement,Dictionary]*CoefMatrix).^2))/numel(Data));
            disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);
        else
            output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2);
            disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);
        end
    end
    if (displayErrorWithTrueDictionary )
        [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));
    end
    Dictionary = I_clearDictionary(Dictionary);
    
    if (isfield(param,'waitBarHandle'))
        waitbar(iterNum/param.counterForWaitBar);
    end
end

if (displayErrorWithTrueDictionary )
    output.ratio = ratio;
end
output.CoefMatrix = CoefMatrix;
Dictionary = [FixedDictionaryElement,Dictionary];
output.time = time/iterNum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findBetterDictionaryElement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,errT)
if (length(who('numCoefUsed'))==0)
    numCoefUsed = 1;
end
relevantDataIndices = find(CoefMatrix(j,:)); % the data indices that uses the j'th dictionary element.
if (length(relevantDataIndices)<2*errT) %(length(relevantDataIndices)==0)
    betterDictionaryElement = Dictionary(:,j);%ErrorMat(:,i); %
    NewVectorAdded = 0;
    return;
end

NewVectorAdded = 0;
tmpCoefMatrix = CoefMatrix(:,relevantDataIndices);
tmpCoefMatrix(j,:) = 0;% the coeffitients of the element we now improve are not relevant.
errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix); % vector of errors that we want to minimize with the new element
% % the better dictionary element and the values of beta are found using svd.
% % This is because we would like to minimize || errors - beta*element ||_F^2.
% % that is, to approximate the matrix 'errors' with a one-rank matrix. This
% % is done using the largest singular value.
%%%%%%------exact SVD computation-------
[betterDictionaryElement,singularValue,betaVector] = svds(errors,1);
CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';% *signOfFirstElem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findDistanseBetweenDictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)
% first, all the column in oiginal starts with positive values.
catchCounter = 0;
totalDistances = 0;
for i = 1:size(new,2)
    new(:,i) = sign(new(1,i))*new(:,i);
end
for i = 1:size(original,2)
    d = sign(original(1,i))*original(:,i);
    distances =sum ( (new-repmat(d,1,size(new,2))).^2);
    [minValue,index] = min(distances);
    errorOfElement = 1-abs(new(:,index)'*d);
    totalDistances = totalDistances+errorOfElement;
    catchCounter = catchCounter+(errorOfElement<0.01);
end
ratio = 100*catchCounter/size(original,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  I_clearDictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D_Trim,idx] = I_clearDictionary(D)
T=0.99;
K= size(D,2);

G=D'*D; G(eye(K)==1)= 0;

IDX = 1:K;
mark = false(1,K);

count = 1;
while sum(~mark)>count
    idx = IDX(~mark);
    mark = mark|(G(idx(count),:)>T);
    count = count + 1;
end
idx = IDX(~mark);
D_Trim = D(:,idx);
return