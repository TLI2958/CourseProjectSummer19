% In this script, I will implement Treisman's visual search paradigm, prune
% data, visualize relationships between variables, and save data.

% Reasons for choosing it:
% a) I made a raw version of this paradigm at the beginning of the class.
% Now I would like to update it with the codes learned throughout the
% 1-month class;
% b) Coding skills for Cognitive Psychology, data analysis(permutation), and graphing 
% will be combined in this script, though GUI is not indicated (I did assignment5).
% c) I would like to use it as a template for experiments in CogPsy field;
% d) I have no access to psychophysical data. 

% There are 3 independent variables(IV),
% a) condition: popout vs. conjunction (blocked);
% Popout: target and non-targets in each trial 
% b) set size:4,8,12,16 (blocked); referring to total number of symbols in one trial;
% c) target: present vs. absent (randomized);
% The dependent variable (DV) is the reaction time(RT);
% The hypothesis goes as: 
% a) The condition under which stimuli are presented will modulate speed of
% reaction. Specifically, participants will respond quicker under popout
% condition;
% b) Whether target is present will modulate reaction speed. Specifically, when
% target is present, participants will respond quicker.
% c) The condition under which stimuli are presented
% will modulate the relationship between reaction time and set size. Specifically,
% under popout condition, reaction time will be independent of set size,
% while under conjunction condition, reaction time will increase with the
% increase of set size.
% NB: only trials with correct response will enter data analysis process.

% Input: N/A; The reviewers need to go through the experiment on their own
% and collect their own data to implement processes following the experiment.
% Output: a) reaction time
%         b) statistics: permutation will be employed to test hypothesis a)
%         & b)
%         c) plots visualizing the correlation between reaction time and set
%         size under popout and conjunction conditions, respectively.

% author: Tiancheng Li; contact: tl2546@nyu.edu
% June ~ July,2019

%% 0 Initialization
close all % close all the figures
clear all % clear workspace
clc % clear command window

%% 1 Priors
% In this section, we will define constants, and create variables to store
% values which will be used/collected later.

% 1a) Constants
n_group = 30; % number of trials for each size;
n_total = 120; % total number of trials;
size_all = [4,8,12,16]; % 4 different levels of set size;
Fsize = 22; % fontsize;
hair = .5; % the room between axis and data;

% 1b) Matrices and cells to store data
CorrectAnswer = cell(120,2); % create an empty matrix to store correct answers
ActualResponse = cell(120,2); % create empty cell array to store actual responses of participant;
RT = zeros(120,2); % Empty matrix to store response time;

% 1c) Below are stimuli
target_sym = 'o'; % target symbol 'o';
non_target_sym = 'x'; % non target symbol 'x';
symbol = ['x';'o']; % collection of all the symbols;
t_color = 'b';% target color blue;
non_t_color = 'g'; % non-target color green;

%% 1-1 Setting backgroup and Instructions
h = figure; % call a figure
h.Color = 'w'; % set background color to white
set(gca,'Visible','off'); % remove axes
ScreenSize = get(0,'screensize');
set(h,'Position',ScreenSize);
initialInstruction = text(.5,.5,sprintf('Welcome to the "Blue Circle" experiment.\n\nIn each trial, one image will be presented at the center of the screen.\nYou are expected to judge whether there exists a "Blue Circle".\nPress Z if there is; M if there is not.\n\n Please respond as quick as possible but avoid random guessing.\nPress any key to continue.\n'),'Color','k','FontSize',Fsize);
initialInstruction.HorizontalAlignment = 'center'; % center the instructions
pause % wait for response
close all % close the instruction page


%% 2 Presenting Stimuli
% In this section, we will present all the trials of the experiment.
% Design: 2(condition: Popout vs. Conjunction) x 4(number of stimuli: 4, 8, 12, 16) x 2(Target: Presence vs. Absent)

rng('shuffle'); % seed

% 2a) Target: Present vs. Absent
for ii = 1:n_group % all the 120 trials
    for jj = 1:2
        if rand > .5 % the condition where target is present
            CorrectAnswer{ii,jj} = 'z'; % z corresponds to present
        else
            CorrectAnswer{ii,jj} = 'm'; % m corresponds to absent
        end % end of if loop
    end % end of loop for columns
end % end of loop for rows

% 2b) Set size: randomized
randSetSize_popout = randperm(length(size_all));
randSetSize_conjunction = randperm(length(size_all));
% 2c) Presenting stimuli in blocks
countTrial = 0; % set the initial value for number of trials

if randperm(2,1) == 1 % first, conterbalancing sequence of Popout and Conjunction Conditions
    % if randperm(2,1) = 1, display Popout first; if randperm(2,1) = 2,
    % display Conjunction first.

    for ss = 1:length(size_all)
        rs = randSetSize_popout(ss); % select corresponding set size generated randomly
        location_stim = rand(n_group,size_all(rs),2); %location of stimuli between 0 and 1 for each number of stimuli
        for ii = 1:n_group % all the 120 trials
            for jj = 1:2
                if rand > .5 % the condition where target is present
                    CorrectAnswer{ii,jj} = 'z'; % z corresponds to present
                else
                    CorrectAnswer{ii,jj} = 'm'; % m corresponds to absent
                end % end of if loop
            end % end of loop for columns
        end % end of loop for rows
        for ii = 1:n_group % 30 trials per condition
            figure
            set(gcf,'Color','w'); 
            axis off % call a figure to present each trial, set features.
            tic % start timing
            if strcmp(CorrectAnswer{ii,1},'z') == 1 % if target is present
                text(location_stim(ii,1,1),location_stim(ii,1,2),target_sym,...
                    'fontsize',Fsize,'color',t_color); % present stimuli
                symbol_all = [repmat('x',((size(location_stim,2)/2) - 1),1);...
                    repmat('o',size(location_stim,2)/2,1)]; % store all the remaining symbols
                symbol_index = randperm(size(location_stim,2) - 1); % index for the symbols. Target has used one symbol, so we need to minus 1.
                symbol_all = symbol_all(symbol_index); % convert index to symbols
                for jj = 2:size(location_stim,2) % 2nd to the last column
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),symbol_all(jj-1),...
                        'fontsize',Fsize,'color',non_t_color);
                    % present all the nontarget stimuli;
                end
            else % if target is absent
                symbol_all = [repmat('x',size(location_stim,2)/2,1);...
                    repmat('o',size(location_stim,2)/2,1)]; % replicate 'x' and 'o' with (half of set size) times
                symbol_index = randperm(size(location_stim,2)); % randomize location of symbols by randomizing their sequence
                symbol_all = symbol_all(symbol_index); % convert index to symbols
                for jj = 1:size(location_stim,2) 
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),symbol_all(jj),...
                        'fontsize',Fsize,'color',non_t_color); % present all stimuli;
                end
            end
            pause(2) % pause for 2000 ms for presentation of stimuli
            RT(ii,1) = toc; % collect reaction time
            ActualResponse{ii} = get(gcf,'CurrentCharacter'); % collect actual response
            pause(.1); % pause for 100 ms between trials;
            countTrial = countTrial + 1; % count number of trials to ensure the break of loop after 120 trials
            close all % close all the figures
        end
        interval = figure; % short break between blocks with different numbers of stimuli;
        interval.Color = 'w'; % set background color to white;
        hold on
        ttinter = text(.5,.5,'Hold on! Take a short break!','FontSize',Fsize); 
        ttinter.HorizontalAlignment = 'center'; 
        axis off
        pause(30) % pause for 30 s between blocks
        close all
    end
    for ss = 1:length(size_all) % 4 different set sizes
        rs = randSetSize_conjunction(ss);
        location_stim = rand(n_group,size_all(rs),2); %location of stimuli between 0 and 1 for each set size
        for ii = 1:n_group % all the 120 trials
            for jj = 1:2 % 2 conditions: Popout vs. Conjunction
                if rand > .5 % the condition where target is present
                    CorrectAnswer{ii,jj} = 'z'; % z corresponds to present
                else
                    CorrectAnswer{ii,jj} = 'm'; % m corresponds to absent
                end % end of if loop
            end % end of loop for columns
        end % end of loop for rows
        for ii = 1:n_group % 30 trials per condition
            figure
            set(gcf,'Color','w');
            axis off
            tic % timing
            if strcmp(CorrectAnswer{ii,2},'z') == 1 % if target is present
                text(location_stim(ii,1,1),location_stim(ii,1,2),target_sym,'fontsize',Fsize,'color',t_color); % blue 'o'(target) at the first location
                text(location_stim(ii,size(location_stim,2),1),location_stim(ii,size(location_stim,2),2),...
                    non_target_sym,'fontsize',Fsize,'color',non_t_color); % green non-target at the last location
                for jj = 2:size(location_stim,2)/2 % number of green 'o's
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),target_sym,'fontsize',Fsize,'color',non_t_color); % green 'o's 
                end
                for jj = (size(location_stim,2)/2) + 1:(size(location_stim,2) - 1) % number of Conjunction, blue 'x's 
                    % i.e., nontargets with target color
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),non_target_sym,'fontsize',Fsize,'color',t_color); % blue 'x's.
                end
            else % if target is absent
                symbol_all = [repmat('x',size(location_stim,2)/2,1);repmat('o',size(location_stim,2)/2,1)]; % half 'x's, half 'o's.
                symbol_index = randperm(size(location_stim,2)); 
                symbol_all = symbol_all(symbol_index);
                for jj = 1:size(location_stim,2)
                    if strcmp(symbol_all(jj),target_sym) == 1 % if the symbol is 'o'
                        text(location_stim(ii,jj,1),location_stim(ii,jj,2),symbol_all(jj),...
                            'fontsize',Fsize,'color',non_t_color); % all 'o's in green
                    else % 'x'
                        text(location_stim(ii,jj,1),location_stim(ii,jj,2),...
                            symbol_all(jj),'fontsize',Fsize,'color',t_color); % all 'x's in blue (complete conjunction)
                    end
                end
            end
            pause(2) % pause 2000 ms for response
            RT(ii,2) = toc; % store RTs
            ActualResponse{ii} = get(gcf,'CurrentCharacter'); % collect actual response
            pause(.1); % pause for 100 ms between trials
            countTrial = countTrial + 1;
            close all
        end
        if countTrial < n_total
            interval = figure; % short break between blocks with different numbers of stimuli;
            interval.Color = 'w'; % set background color to white;
            hold on
            ttinter = text(.5,.5,'Hold on! Take a short break!','FontSize',Fsize);
            ttinter.HorizontalAlignment = 'center';
            axis off
            pause(30) % pause for 30 s between blocks
            close all
        else
            break
        end
    end
else
    for ss = 1:length(size_all)
        rs = randSetSize_popout(ss); % select corresponding set size generated randomly
        location_stim = rand(n_group,size_all(rs),2); %location of stimuli between 0 and 1 for each number of stimuli
        for ii = 1:n_group % all the 120 trials
            for jj = 1:2
                if rand > .5 % the condition where target is present
                    CorrectAnswer{ii,jj} = 'z'; % z corresponds to present
                else
                    CorrectAnswer{ii,jj} = 'm'; % m corresponds to absent
                end % end of if loop
            end % end of loop for columns
        end % end of loop for rows
        for ii = 1:n_group % 30 trials per condition
            figure
            set(gcf,'Color','w');
            axis off
            tic % timing
            if strcmp(CorrectAnswer{ii,2},'z') == 1 % if target is present
                text(location_stim(ii,1,1),location_stim(ii,1,2),target_sym,'fontsize',Fsize,'color',t_color);
                text(location_stim(ii,size(location_stim,2),1),location_stim(ii,size(location_stim,2),2),...
                    non_target_sym,'fontsize',Fsize,'color',non_t_color);
                for jj = 2:size(location_stim,2)/2
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),target_sym,'fontsize',Fsize,'color',non_t_color);
                end
                for jj = (size(location_stim,2)/2) + 1:(size(location_stim,2) - 1)
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),non_target_sym,'fontsize',Fsize,'color',t_color);
                end
            else
                symbol_all = [repmat('x',size(location_stim,2)/2,1);repmat('o',size(location_stim,2)/2,1)];
                symbol_index = randperm(size(location_stim,2));
                symbol_all = symbol_all(symbol_index);
                for jj = 1:size(location_stim,2)
                    if strcmp(symbol_all(jj),target_sym) == 1
                        text(location_stim(ii,jj,1),location_stim(ii,jj,2),...
                            symbol_all(jj),'fontsize',Fsize,'color',non_t_color);
                    else
                        text(location_stim(ii,jj,1),location_stim(ii,jj,2),...
                            symbol_all(jj),'fontsize',Fsize,'color',t_color);
                    end
                end
            end
            pause(2)
            RT(ii,2) = toc;
            ActualResponse{ii} = get(gcf,'CurrentCharacter'); % collect actual response
            pause(.1); % pause for 100 ms between trials
            countTrial = countTrial + 1;
            close all
        end
        interval = figure; % short break between blocks with different numbers of stimuli;
        interval.Color = 'w'; % set background color to white;
        hold on
        ttinter = text(.5,.5,'Hold on! Take a short break!','FontSize',Fsize);
        ttinter.HorizontalAlignment = 'center';
        axis off
        pause(30) % pause for 30 s between blocks
        close all
    end
    for ss = 1:length(size_all)
        rs = randSetSize_popout(ss); % select corresponding set size generated randomly
        location_stim = rand(n_group,size_all(rs),2); %location of stimuli between 0 and 1 for each number of stimuli
        for ii = 1:n_group % all the 120 trials
            for jj = 1:2
                if rand > .5 % the condition where target is present
                    CorrectAnswer{ii,jj} = 'z'; % z corresponds to present
                else
                    CorrectAnswer{ii,jj} = 'm'; % m corresponds to absent
                end % end of if loop
            end % end of loop for columns
        end % end of loop for rows
        for ii = 1:n_group % 30 trials per condition
            figure
            set(gcf,'Color','w');
            axis off
            tic
            if strcmp(CorrectAnswer{ii,1},'z') == 1 % if target is present
                text(location_stim(ii,1,1),location_stim(ii,1,2),target_sym,'fontsize',Fsize,'color',t_color);
                symbol_all = [repmat('x',((size(location_stim,2)/2) - 1),1);repmat('o',size(location_stim,2)/2,1)]; % store all the remaining symbols
                symbol_index = randperm(size(location_stim,2) - 1); % index for the symbols
                symbol_all = symbol_all(symbol_index);
                for jj = 2:size(location_stim,2) % 2nd to the last column
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),symbol_all(jj-1),'fontsize',Fsize,'color',non_t_color);
                    % present all the nontarget stimuli;
                end
            else
                symbol_all = [repmat('x',size(location_stim,2)/2,1);repmat('o',size(location_stim,2)/2,1)];
                symbol_index = randperm(size(location_stim,2));
                symbol_all = symbol_all(symbol_index);
                for jj = 1:size(location_stim,2)
                    text(location_stim(ii,jj,1),location_stim(ii,jj,2),symbol_all(jj),...
                        'fontsize',Fsize,'color',non_t_color); % present all stimuli;
                end
            end
            pause(2)
            RT(ii,1) = toc;
            ActualResponse{ii} = get(gcf,'CurrentCharacter'); % collect actual response
            pause(.1); % pause for 100 ms between trials;
            countTrial = countTrial + 1; % count num of trials
            close all % close all figures
        end
        if countTrial < n_total % <120
            interval = figure; % short break between blocks with different numbers of stimuli;
            interval.Color = 'w'; % set background color to white;
            hold on
            ttinter = text(.5,.5,'Hold on! Take a short break!','FontSize',Fsize);
            ttinter.HorizontalAlignment = 'center';
            axis off
            pause(30) % pause for 30 s between blocks
            close all
        else
            break % if reach 120, break the loop and terminate the experiment.
        end
    end
end

% 2c) Displaying Endwords
close all % clear stimuli
figure % open a new figure
set(gcf,'Color','w'); % set background color to white
set(gcf,'Position',ScreenSize); % full-screen size
axis off
endWord = text(.5,.5,sprintf('Now you have reached the end of the experiment.\nThank you so much for your participation.\nGoodbye and Have a nice day!'),'Color','k','FontSize',Fsize);
endWord.HorizontalAlignment = 'center'; % center the instructions

%% 3 Pruning
% 3a) Error rate
responseNumel = zeros(120,2); % store whether response is correct or not
for ii = 1:n_total
    for jj = 1:2  % 120 trials per condition
        if strcmp(CorrectAnswer{ii,jj},ActualResponse{ii,jj}) == 1
            responseNumel(ii,jj) = 1; % correct as 1
        else
            responseNumel(ii,jj) = 0; % incorrect as 0
        end
    end
end

errorRate = length(responseNumel(responseNumel == 1))./length(responseNumel); % calculate error rates
% In actual experiment, participants with error rates > 0.02 will not enter
% data analysis process

% 3b) actual response
% response other than z/m will be recoded as "nan".
responseValid = cell(120,2); % create a cell array to store all the valid responses
for ii = 1:n_total
    for jj = 1:2 % 120 trials/condition
        if strcmp(ActualResponse{ii,jj},'z') == 1
            responseValid{ii,jj} = ActualResponse{ii,jj}; % z
        elseif strcmp(ActualResponse{ii,jj},'m') == 1
            responseValid{ii,jj} = ActualResponse{ii,jj}; % m
        else
            responseValid{ii,jj} = nan; % others including other keys and no response
        end
    end
end

% 3c) response time
RT_valid = zeros(120,2); 
for ii = 1:n_total
    for jj = 1:2 % 120 trials/condition
        if RT(ii,jj) > 2 % no response + the timing function running time
            RT_valid(ii,jj) = nan;
        else
            RT_valid(ii,jj) = RT(ii,jj);
        end
    end
end

for ii = 1:n_total
    for jj = 1:2
        if isnan(responseValid{ii,jj}) == 1
            RT_valid(ii,jj) = nan; % assign RTs for trials with invalid response to nan.
        end
    end
end

% Now RT_valid should only record RTs for correct trials. All other values
% will be recoded as nans.

% 3d) participant-wise elimination
sumnans = sum(isnan(RT_valid),2); % sum up nans per row
RT_real = RT_valid(sumnans == 0); % keep rows without nans
%% 4 Data Analyzing and Visualizing
% Now that we should have stored all the valid RTs for Popout condition in the
% 1st column and all the valid RTs for Conjunction condition in the 2nd
% column.

% 4a) hypothesis a: mean difference of RTs between popout and conjunction
% conditions*
% * for space's sake, use permutation in this test only 
% Also, assume RT_xx stores mean RTs for all participants

RT_Popout = mean(RT_real(:,1)); 
RT_Conjunction = mean(RT_real(:,2)); % means for two conditions 
RT_collapsed = [RT_Popout;RT_Conjunction];

n1 = length(RT_real); % popout
n2 = 2*length(RT_real); % conjunction
numSampling = 1e3; % times of permutation
meanDiff_Condition = RT_Conjunction - RT_Popout; % empirical mean difference
meanPerm = zeros(numSampling,1); % empty matrix to store means

rng('shuffle'); % seed
for ii = 1:numSampling
    randSeq = (randperm(n2))'; % transpose the sequence into one column
    permPopout = randSeq(1:n1); % select as RTs for Popout
    permPopout = mean(permPopout); % mean RT Popout
    permConjunction = randSeq(n1+1:n2); % RTs for Conjunction
    permConjunction = mean(permConjunction); % mean RT conjunction
    meanPerm(ii) = permConjunction - permPopout; 
end 

numBins = 100;
LW = 8; % linewidth
figure
set(gcf,'Color','w');
set(gca,'TickDir','out');
box off % set features of the plot
histogram(sort(meanPerm),numBins); % plot means in the ascending way
line([meanDiff_Condition,meanDiff_Condition],[ylim],'LineWidth',LW,'Color','r'); % mark the empirical differece
meanPermsorted = sort(meanPerm); 
propNull = sum(meanPermsorted > meanDiff_Condition)./length(meanPermsorted).*100; % proportion due to chance
% Decision: reject null hypothesis or not, given probability value.

% 4b) hypothesis b: mean difference of RTs between target present and
% absent conditions
numPresent = zeros(length(RT_valid),2); 
numAbsent = numPresent; % empty matrices to store whether target is pre/ab-sent
for ii = 1:length(RT_valid)
    for jj = 1:2
        if CorrectAnswer{ii,jj} == 'z' % present
            if sumnans(ii) == 0 % valid response
            numPresent(ii,jj) = RT_valid(ii,jj); % store RT
            end
        else % absent
            if sumnans(ii) == 0 % valid response
            numAbsent(ii,jj) = RT_valid(ii,jj); 
            end
        end
    end
end
numPresent = numPresent(numPresent ~=0);
numAbsent = numAbsent(numAbsent ~=0); % keep valid responses only

% I will assume HOV, norm distribution and use t-statistic
[~,p_targetPresent,CI_targetPresent,tSTATS_targetPresent] = ttest(numPresent,numAbsent); % do repreated t-test 

% 4c) correlations between number of set size and mean RTs under popout and conjunction conditions, respectively 
meanSetSize_Popout = [mean(RT_valid(1:n_group,1),'omitnan'),mean(RT_valid(n_group+1:2*n_group,1),'omitnan'),...
    mean(RT_valid(2*n_group + 1:3*n_group,1),'omitnan'),mean(RT_valid(3*n_group + 1:4*n_group,1),'omitnan')]; 
meanSetSize_Conjunction = [mean(RT_valid(1:n_group,2),'omitnan'),mean(RT_valid(n_group+1:2*n_group,2),'omitnan'),...
    mean(RT_valid(2*n_group + 1:3*n_group,2),'omitnan'),mean(RT_valid(3*n_group + 1:4*n_group,2),'omitnan')]; % mean for each set size
[~,Index1] = sort(randSetSize_popout); 
[~,Index2] = sort(randSetSize_conjunction); % indices of set size in the ascending way
meanSetSize_Popout = meanSetSize_Popout(Index1);
meanSetSize_Conjunction = meanSetSize_Conjunction(Index2); % re-arrange means so that correspond to set size = 4,8,12,16.
% Should include mean RTs from many participants. Here, I assume the above
% two matrices of means are means of all the participants
figure
plot(Index1,meanSetSize_Popout,'LineWidth',LW,'Color','r');
hold on 
plot(Index2,meanSetSize_Conjunction,'LineWidth',LW,'Color','b'); % correlations between set size and mean RTs for Popout & Conjunction
set(gca,'TickDir','out')
legend({'Popout','Conjunction'})
xlabel('set size')
ylabel('mean RTs')
title('Relationship Between Set Size and Reaction Time for Different Conditions') % set features of the plot

indexFullPop = rep(index1,n_group); 
indexFullCon = rep(index2,n_group); % replicate indices (of order) for 30 times, each
% Will provide codes for this function
indexFullPop = indexFullPop(sumnan == 0); 
indexFullCon = indexFullCon(sumnan == 0); % keep only valid responses

corrRT_Popout = corrcoef(indexFullPop,RT_real(:,1)); % should be no diff from 0
corrRT_Conjunction = corrcoef(indexFullCon,RT_real(:,2)); % should be diff from 0
reps = length(RT_real); % will do a repeated anova, so each mean RT as a group 
[~,p_SetSize,CI_SetSize,fSTATS_SetSize] = anova2(RT_real,reps); % do repeated anvoa to assess diff between popout and conjunction conditions
