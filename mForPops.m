% %////////////////////////// M MATRIX PLASTIC /////////////////////////////%
% %Evaluate M matrix for static population of 2 output units
%
% mNum = 100;
%
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_1_Pop.mat');
%
% numY = numel(Pop(1).Y);
% numZ = numel(Pop(1).Z);
% numy = numel(Pop(1).y);
% numz = numel(Pop(1).z);
%
% rand_a = rand(mNum,1);
% rand_m = randi(numY + numZ, mNum, 1);
%
% for i=1:30
%     i
%     %load Pop%
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_Pop.mat']);
%     %evaluate unmutated and mutated phenotypes
%     [uPhens,mPhens,recChange] = computeM(Pop,2,rand_a,rand_m);
%     %save unmutated and mutated phenotypes
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_mutational_effects.mat'],'uPhens','mPhens','recChange');
% end
% %////////////////////////////////// END /////////////////////////////////%

% %////////////////////////// M MATRIX STATIC /////////////////////////////%
% %Evaluate M matrix for static population of 2 output units
%
% %Get a random set of values using Sobol sequences
%
% mNum = 100;
%
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_1_Pop.mat');
%
% numY = numel(Pop(1).Y);
% numZ = numel(Pop(1).Z);
% numy = numel(Pop(1).y);
% numz = numel(Pop(1).z);
%
% rand_a = rand(mNum,1);
% rand_m = randi(numY + numZ, mNum, 1);
%
% for i=1:36
%     i
%     %load Pop%
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_Pop.mat']);
%     %evaluate unmutated and mutated phenotypes
%     [uPhens,mPhens,recChange] = computeM(Pop,0,rand_a,rand_m);
%     %save unmutated and mutated phenotypes
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_mutational_effects.mat'],'uPhens','mPhens','recChange');
% end
% %////////////////////////////////// END /////////////////////////////////%

% %////////////////////////// M MATRIX STATIC WITH CUES /////////////////////////////%
% %Evaluate M matrix for static population of 2 output units when they
% %receive environmental cues
%
% %Get a random set of values using Sobol sequences
%
% mNum = 100;
%
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Static_Cue\10000\replicate_1_Pop.mat');
%
% numY = numel(Pop(1).Y);
% numZ = numel(Pop(1).Z);
% numy = numel(Pop(1).y);
% numz = numel(Pop(1).z);
%
% rand_a = rand(mNum,1);
% rand_m = randi(numY + numZ, mNum, 1);
%
% for i=1:10
%     %load Pop%
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static_Cue\10000\replicate_' num2str(i) '_Pop.mat']);
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static_Cue\10000\replicate_' num2str(i) '_phi.mat']);
%     %evaluate unmutated and mutated phenotypes
%     [uPhens,mPhens,recChange] = computeM(Pop,5,rand_a,rand_m,phi);
%     %save unmutated and mutated phenotypes
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static_Cue\10000\replicate_' num2str(i) '_mutational_effects.mat'],'uPhens','mPhens','recChange');
% end
% %////////////////////////////////// END /////////////////////////////////%

% %////////////////////////// M MATRIX INFERENCE /////////////////////////////%
% %Evaluate M matrix for static population of 2 output units when they
% %receive environmental cues
% 
% %Get a random set of values using Sobol sequences
% 
% mNum = 300;
% 
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_1_Pop.mat');
% 
% numY = numel(Pop(1).Y);
% numZ = numel(Pop(1).Z);
% numy = numel(Pop(1).y);
% numz = numel(Pop(1).z);
% 
% rand_a = rand(mNum,1);
% rand_m = randi(numY + numZ, mNum, 1);
% 
% for i=1:20
%     %load Pop%
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_' num2str(i) '_Pop.mat']);
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_' num2str(i) '_phi.mat']);
%     %evaluate unmutated and mutated phenotypes
%     [uPhens,mPhens,recChange] = computeM(Pop,2,rand_a,rand_m,phi);
%     %save unmutated and mutated phenotypes
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_' num2str(i) '_mutational_effects_C.mat'],'uPhens','mPhens','recChange');
% end
% %////////////////////////////////// END /////////////////////////////////%

% 
% %/////////////////////////////DISPLAY M - INFERENECE //////////////////////////////
% mNum = 300;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=1:20
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_' num2str(i) '_mutational_effects_C.mat']);
%     B = mPhens - uPhens;
%     tempCorr = NaN(N,3);
%     
%     for k=1:N
%         tempB = B((k-1)*mNum+1:k*mNum,:); %foreach individual
%         tempB = tempB(sum(~isinf(tempB),2)==3 & sum(~isnan(tempB),2)==3,:);
%         tempCorr(k,1) = corr(tempB(:,1),tempB(:,2));
%         tempCorr(k,2) = corr(tempB(:,2),tempB(:,3));
%         tempCorr(k,3) = corr(tempB(:,1),tempB(:,3));
%     end
%     %mean(tempCorr(~isnan(tempCorr)))
%     tempCorr = mean(tempCorr);
%     meanCorr = [meanCorr; tempCorr]
%     
%     tempB = B(sum(~isinf(B),2)==3 & sum(~isnan(B),2)==3,:);
%     tempVar = sum(tempB.^2)/size(tempB,1);
%     meanVar = [meanVar; mean(tempVar)];
%     
%     fig=figure;
%     subplot(1,3,1);
%     hold on;
%     title(['Var: ' num2str(tempVar(1)) '  Corr:' num2str(tempCorr(1))]);
%     plot(B(logical(~recChange),1),B(logical(~recChange),2),'.');
%     plot(B(logical(recChange),1),B(logical(recChange),2),'.','color',[1 .5 0]);
%     hold off;
%     xlim([-1 1]); ylim([-1 1]);
%     xlabel('Trait A'); ylabel('Trait B');
%     axis square;
%     subplot(1,3,2);
%     hold on;
%     title(['Var: ' num2str(tempVar(2)) '  Corr:' num2str(tempCorr(2))]);
%     plot(B(logical(~recChange),2),B(logical(~recChange),3),'.');
%     plot(B(logical(recChange),2),B(logical(recChange),3),'.','color',[1 .5 0]);
%     hold off;
%     xlim([-1 1]); ylim([-1 1]);
%     xlabel('Trait B'); ylabel('Trait C');
%     axis square;
%     subplot(1,3,3);
%     hold on;
%     title(['Var: ' num2str(tempVar(3)) '  Corr:' num2str(tempCorr(3))]);
%     plot(B(logical(~recChange),1),B(logical(~recChange),3),'.');
%     plot(B(logical(recChange),1),B(logical(recChange),3),'.','color',[1 .5 0]);
%     hold off;
%     xlim([-1 1]); ylim([-1 1]);
%     xlabel('Trait A'); ylabel('Trait C');
%     axis square;
%     
%     print(fig,['Inference_mut_distribution_C_' num2str(i)],'-dpng');
%     
% end
% [mean(meanVar) mean(meanCorr)]




% for i=1:36
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_Pop.mat']);
%     [A,B] = computeG(Pop,0);
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_parents_traits.mat'],'A');
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_child_traits.mat'],'B');
% end

% for i=1:10
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50k\2outputs\replicate_' num2str(i) '_Pop.mat']);
%     [A,B] = computeG(Pop,2);
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50k\2outputs\replicate_' num2str(i) '_parents_traits.mat'],'A');
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50k\2outputs\replicate_' num2str(i) '_child_traits.mat'],'B');
% end
%
% for i=21:30
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_Pop.mat']);
%     [pTraits,cTraits] = computeG(Pop,2);
%     save(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_genetic_effects.mat'],'pTraits','cTraits');
% end

% %/////////////////////////////////////////////////////////////////////////
% mNum = 100;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=1:30
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_mutational_effects.mat']);
%     B = mPhens - uPhens;
% 	tempCorr = [];
%
%     for k=1:N
%         B1 = B((k-1)*mNum+1:k*mNum,:); %foreach individual
%         B1 = B1(sum(~isinf(B1),2)==2 & sum(~isnan(B1),2)==2,:);
%         tempCorr = [tempCorr corr(B1(:,1),B1(:,2))];
%     end
%     meanCorr = [meanCorr mean(tempCorr(~isnan(tempCorr)))];
%
%     B = B(sum(~isinf(B),2)==2 & sum(~isnan(B),2)==2,:);
%     meanVar = [meanVar mean(sum(B.^2)/size(B,1))];
%
% end
% [mean(meanVar) mean(meanCorr)]


% %/////////////////////////////////////////////////////////////////////////
% mNum = 100;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=1:10
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static_Cue\10000\replicate_' num2str(i) '_mutational_effects.mat']);
%     B = mPhens - uPhens;
% 	tempCorr = [];
%
%     for k=1:N
%         B1 = B((k-1)*mNum+1:k*mNum,:); %foreach individual
%         B1 = B1(sum(~isinf(B1),2)==2 & sum(~isnan(B1),2)==2,:);
%         tempCorr = [tempCorr corr(B1(:,1),B1(:,2))];
%     end
%     meanCorr = [meanCorr mean(tempCorr(~isnan(tempCorr)))];
%
%     B = B(sum(~isinf(B),2)==2 & sum(~isnan(B),2)==2,:);
%     meanVar = [meanVar mean(sum(B.^2)/size(B,1))];
%
% end
% [mean(meanVar) mean(meanCorr)]


% %/////////////////////////////////////////////////////////////////////////
% mNum = 100;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=1:36
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_mutational_effects.mat']);
%     B = mPhens - uPhens;
% 	tempCorr = [];
%
%     for k=1:N
%         B1 = B((k-1)*mNum+1:k*mNum,:); %foreach individual
%         B1 = B1(sum(~isinf(B1),2)==2 & sum(~isnan(B1),2)==2,:);
%         tempCorr = [tempCorr corr(B1(:,1),B1(:,2))];
%     end
%     meanCorr = [meanCorr mean(tempCorr(~isnan(tempCorr)))];
%
%     B = B(sum(~isinf(B),2)==2 & sum(~isnan(B),2)==2,:);
%     meanVar = [meanVar mean(sum(B.^2)/size(B,1))];
%
% end
% [mean(meanVar) mean(meanCorr)]

% % Evolvability - Mean Initial Fitness
% N_Static = 36;
% initFit_Collection_Static = NaN(N_Static,7);
% for i=1:N_Static
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_Pop.mat']);
%     [~,initFit_Collection_Static(i,:)] = evalEvolvability(Pop,0);
% end
% N_Hetero = 10;
% initFit_Collection_Hetero = NaN(N_Hetero,7);
% for i=1:N_Hetero
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Heterogeneous\replicate_' num2str(i) '_Pop.mat']);
%     [~,initFit_Collection_Hetero(i,:)] = evalEvolvability(Pop,0);
% end
% N_Plastic = 30;
% initFit_Collection_Plastic = NaN(N_Plastic,7);
% for i=1:N_Plastic
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_Pop.mat']);
%     [~,initFit_Collection_Plastic(i,:)] = evalEvolvability(Pop,2);
% end
% fig = figure; hold on; plot(mean(initFit_Collection_Static,1),'x'); plot(mean(initFit_Collection_Plastic,1),'o'); plot(mean(initFit_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; print(fig,['evol_s_' num2str(50)],'-dpng');





% % % %%%%%%%%%%%%%%Evolvability%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeats = 1;
% N_Static = 36;
% initFit_Collection_Static = NaN(Repeats * N_Static,7);
% fitDiff_Collection_Static = NaN(Repeats * N_Static,4,7);
% for i=1:N_Static
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_Pop.mat']);
%     for j = 1: Repeats
%         [fitDiff_Collection_Static(i,:,:),initFit_Collection_Static(i,:)] = evalEvolvability(Pop,0);
%     end
% end
% N_Hetero = 10;
% initFit_Collection_Hetero = NaN(Repeats * N_Hetero,7);
% fitDiff_Collection_Hetero = NaN(Repeats * N_Hetero,4,7);
% for i=1:N_Hetero
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Heterogeneous\replicate_' num2str(i) '_Pop.mat']);
%     for j = 1: Repeats
%         [fitDiff_Collection_Hetero(i,:,:),initFit_Collection_Hetero(i,:)] = evalEvolvability(Pop,0);
%     end
% end
% N_Plastic = 30;
% initFit_Collection_Plastic = NaN(Repeats * N_Plastic,7);
% fitDiff_Collection_Plastic = NaN(Repeats * N_Plastic,4,7);
% for i=1:N_Plastic
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_Pop.mat']);
%     for j = 1: Repeats
%     [fitDiff_Collection_Plastic(i,:,:),initFit_Collection_Plastic(i,:)] = evalEvolvability(Pop,2); %%%%%%%% 222 %%%%%%%
%     end
% end
% %
% %fig = figure; hold on; plot(mean(initFit_Collection_Static,1),'x'); plot(mean(initFit_Collection_Plastic,1),'o'); plot(mean(initFit_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; print(fig,['evol_s_' num2str(50)],'-dpng');
% %fig = figure; hold on; plot(mean(initFit_Collection_Plastic,1),'o'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};
% %
% temp_plastic = mean(fitDiff_Collection_Plastic(:,4,:),1);
% % temp_hetero = mean(fitDiff_Collection_Hetero(:,1,:),1);
% % temp_static = mean(fitDiff_Collection_Static(:,1,:),1);
% plastic_vector_1 = zeros(7,1); for i=1:7, plastic_vector_1(i) = temp_plastic(:,:,i); end
% % static_vector_1 = zeros(7,1); for i=1:7, static_vector_1(i) = temp_static(:,:,i); end
% % hetero_vector_1 = zeros(7,1); for i=1:7, hetero_vector_1(i) = temp_hetero(:,:,i); end
% %
% % fig = figure; hold on; plot(static_vector_1,'x'); plot(plastic_vector_1,'o'); plot(hetero_vector_1,'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};
%
% fig = figure; hold on; plot(plastic_vector_1,'o'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};

%fig = figure; hold on; plot(mean(fitDiff_Collection_Static,1),'x'); plot(mean(fitDiff_Collection_Plastic,1),'o'); plot(mean(fitDiff_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};
% fig = figure; hold on; plot(mean(fitDiff_Collection_Static,1)-mean(initFit_Collection_Static,1),'x'); plot(mean(fitDiff_Collection_Plastic,1)-mean(initFit_Collection_Plastic,1),'o'); plot(mean(fitDiff_Collection_Hetero,1)-mean(initFit_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; print(fig,['evol_s_' num2str(50)],'-dpng');

%fig = figure; hold on; plot(mean(initFit_Collection_Plastic),'o'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};


% initFit_Collection_Plastic = []
% for k=1:1
%     k
% for i=1:1
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_Pop.mat']);
%     [fitDiff, initFit] = evalEvolvability(Pop,2);
%     initFit_Collection_Plastic = [initFit_Collection_Plastic; initFit];
% end
% end

% %/////////////////////////////////////////////////////////////////////////
% C = 499500;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=1:36
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_child_traits.mat']);
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Static\50000\replicate_' num2str(i) '_parents_traits.mat']);
%     D = A-B;
% 	tempCorr = [];
%     D = D(sum(~isinf(D),2)==2 & sum(~isnan(D),2)==2,:);
%     tempCorr = corr(D(:,1),D(:,2));
%     %tempCorr = corr(A(:,2),B(:,1));
%     meanCorr = [meanCorr mean(tempCorr(~isnan(tempCorr)))];
%     meanVar = [meanVar mean(sum(D.^2)/size(D,1))];
% end
% [mean(meanVar) mean(meanCorr)]

% %/////////////////////////////////////////////////////////////////////////
% C = 499500;
% N = 1000;
% meanVar = [];
% meanCorr = [];
% for i=[1:5 11:30]
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Plasticity\50000\replicate_' num2str(i) '_genetic_effects.mat']);
%     B = cTraits-pTraits;
%     ind = sum(~isinf(B),2)==2 & sum(~isnan(B),2)==2;
%     B = B(ind,:);
%     tempCorr = corr(B(:,1),B(:,2));
%     %tempCorr = [corr(cTraits(ind,1)-pTraits(ind,2))];
%     meanCorr = [meanCorr mean(tempCorr(~isnan(tempCorr)))];
%     meanVar = [meanVar mean(sum(B.^2)/size(B,1))];
% end
% [mean(meanVar) mean(meanCorr)]


%fig=figure; hold on; h1=histfit(meanVar_Static); h2=histfit(meanVar_Plastic'); hold off; xlabel('Mutational Variance'); ylabel('Probability Density'); set(gca,'YTickLabel',[]); set(h1(1),'facecolor',[.17 .17 .17]); set(h1(2),'color','r'); set(h2(1),'facecolor',[.57 .57 .57],'facealpha',.5); set(h2(2),'color','g'); legend([h1(1) h2(1)],{'Static','Plastic'});
%phenDiffs = mPhens - uPhens; fig=figure; hold on; plot(phenDiffs(~logical(recChange),1),phenDiffs(~logical(recChange),2),'.','color','blue'); plot(phenDiffs(logical(recChange),1),phenDiffs(logical(recChange),2),'.','color',[1 .5 0]); hold off; axis square; ylim([-1 1]); xlim([-1 1]); xlabel('Trait 1'); ylabel('Trait 2'); title('Static'); print(fig,['mut_distr' num2str(3)],'-dpng');

% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\InitFit_Static.mat');
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\fitDiff_Static.mat');
% temp_static = mean(fitDiff_Static(:,1,:),1);
% static_vector_1 = zeros(7,1); for i=1:7, static_vector_1(i) = temp_static(:,:,i)-mean(initFit_Collection_Static(:,i)); end
% figure; plot(log10(static_vector_1),'x');

% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\InitFit_Plastic.mat');
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\fitDiff_Plastic.mat');
% temp_plastic = mean(fitDiff_Plastic(:,1,:),1);
% plastic_vector_1 = zeros(7,1); for i=1:7, plastic_vector_1(i) = temp_plastic(:,:,i); end
% figure; hold on; plot(plastic_vector_1,'x'); plot(mean(initFit_Collection_Plastic,1),'o'); hold off;

% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\fitDiff_Plastic.mat');
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\fitDiff_Static.mat');
% load('C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\fitDiff_Hetero.mat');
%
% plastic_vector = NaN(4,7);
% static_vector = NaN(4,7);
% hetero_vector = NaN(4,7);
%
% for i=1:4
%     temp_plastic = mean(fitDiff_Plastic(:,i,:),1);
%     temp_hetero = mean(fitDiff_Hetero(:,i,:),1);
%     temp_static = mean(fitDiff_Static(:,i,:),1);
%
%     for j=1:7
%         plastic_vector(i,j) = temp_plastic(:,:,j);
%         static_vector(i,j) = temp_static(:,:,j);
%         hetero_vector(i,j) = temp_hetero(:,:,j);
%     end
% end
%
% fig = figure;
% for i=1:4
%     plastic_vector(i,:) = log(plastic_vector(i,:));%-log(mean(initFit_Collection_Plastic,1));
%     static_vector(i,:) = log(static_vector(i,:));%-log(mean(initFit_Collection_Static,1));
%     hetero_vector(i,:) = log(hetero_vector(i,:));%-log(mean(initFit_Collection_Hetero,1));
%
%     subplot(2,2,i);
%     if i==1, title('t = 1'); end
%     if i==2, title('t = 10'); end
%     if i==3, title('t = 30'); end
%     if i==4, title('t = 100'); end
%     hold on; plot(plastic_vector(i,:),'ko'); plot(static_vector(i,:),'kx'); plot(hetero_vector(i,:),'k*'); hold off; xlabel('Angle of selection'); ylabel('log-fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; legend('Plastic','Static','Hetero');
% end



% % % %%%%%%%%%%%%%%Evolvability Inference%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeats = 1;
% N_Plastic = 20;
% initFit_Collection_Plastic = NaN(Repeats * N_Plastic,10);
% fitDiff_Collection_Plastic = NaN(Repeats * N_Plastic,4,10);
% for i=1:N_Plastic
%     i
%     load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Inference\50000\3outputs\replicate_' num2str(i) '_Pop.mat']);
%     for j = 1: Repeats
%         j
%     [fitDiff_Collection_Plastic((i-1)*Repeats+j,:,:),initFit_Collection_Plastic((i-1)*Repeats+j,:)] = evalEvolvability_Inference(Pop);
%     end
% end
% % %
% %fig = figure; plot(mean(initFit_Collection_Plastic,1),'o'); xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; 
% 
% temp_plastic = mean(fitDiff_Collection_Plastic(:,1,:),1);
% plastic_vector_1 = zeros(10,1); for i=1:10, plastic_vector_1(i) = temp_plastic(:,:,i); end
% % % fig = figure; hold on; plot(static_vector_1,'x'); plot(plastic_vector_1,'o'); plot(hetero_vector_1,'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};
% 
%  fig = figure; hold on; plot(plastic_vector_1,'o'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};

%fig = figure; hold on; plot(mean(fitDiff_Collection_Static,1),'x'); plot(mean(fitDiff_Collection_Plastic,1),'o'); plot(mean(fitDiff_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};
% fig = figure; hold on; plot(mean(fitDiff_Collection_Static,1)-mean(initFit_Collection_Static,1),'x'); plot(mean(fitDiff_Collection_Plastic,1)-mean(initFit_Collection_Plastic,1),'o'); plot(mean(fitDiff_Collection_Hetero,1)-mean(initFit_Collection_Hetero,1),'*'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'}; print(fig,['evol_s_' num2str(50)],'-dpng');

%fig = figure; hold on; plot(mean(initFit_Collection_Plastic),'o'); hold off; xlabel('Angle of selection'); ylabel('Initial Fitness'); ax = gca; ax.XTickLabel = {'45','30','15','0','-15','-30','-45'};

% fitDiffs = [];
% for i=31:56
% load(['C:\Users\Sardokius\Documents\MATLAB\Draghi\Evolvability Assays\Inference\fitDiff_' num2str(i) '.mat']);
% fitDiffs = [fitDiffs; fitDiff_Collection_Plastic];
% end



plastic_vector = NaN(4,10);

for i=1:4
    temp_plastic = mean(fitDiffs(:,i,:),1);

    for j=1:10
        plastic_vector(i,j) = temp_plastic(:,:,j);
    end
end

fig = figure;
for i=1:4
    plastic_vector(i,:) = (plastic_vector(i,:))-(mean(initFit_Collection_Plastic,1));

    subplot(2,2,i);
    if i==1, title('t = 1'); end
    if i==2, title('t = 10'); end
    if i==3, title('t = 30'); end
    if i==4, title('t = 100'); end
    hold on; plot(plastic_vector(i,:),'ko'); xlabel('Directions'); ylabel('log-fitness difference'); xlim([0 11]);
end

