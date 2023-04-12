clear, clc, close all

rng(42)

%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\combined\combined_TCC-FreeSimilarityMatrix-workspace_230214.mat';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\211012_124119_Pollux\210422--211012_Pollux_TCC-FreeSimilarityMatrix-workspace_230222.mat';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\211108_090705_Castor\220517--211108_Castor_TCC-FreeSimilarityMatrix-workspace_230225';
filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\210609_124628_Buster\210428--210609_Buster_TCC-FreeSimilarityMatrix-workspace_230213.mat';
%filepath = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\220823_081207_Morty\220322--220823_Morty_TCC-FreeSimilarityMatrix-workspace_230213.mat';

load(filepath)

nBig = 64;

%%

sm = reshape(x, [64,64]);

% sm = zeros(nBig); % simple simulated cognitive bias
% sm(logical(eye(nBig))) = 1;
% sm(28:32,28:32) = triu(ones(5));
% sm(33:37,33:37) = tril(ones(5));

figure,
imagesc(sm)
axis square
colormap('gray')
colorbar
set(gca,'YDir','normal')

%%

%categoryCenters = [3,38]; %combined
%categoryCenters = [9,37]; %Pollux
%categoryCenters = [4,18,37,50]; %Castor
categoryCenters = [1,43]; %Buster
%categoryCenters = [37,63]; %Morty

% categoryCenters = [32]; %simple simulated cognitive bias


for categoryCenter = categoryCenters(1)

    sm_cs = circshift(sm,[nBig/2-categoryCenter,nBig/2-categoryCenter]);
    % note that that this doesn't put it full on "in the center", but in the
    % 32nd of 64 positions (there are 31 before and 32 after)

    figure,
    imagesc(sm_cs)
    axis square
    colormap('gray')
    colorbar

    xticklabels(str2num(cell2mat(xticklabels)) - (nBig/2-categoryCenter))
    yticklabels(str2num(cell2mat(yticklabels)) - (nBig/2-categoryCenter))

    xline(nBig/2,'r')
    yline(nBig/2,'r')

    figure,
    plot(sm_cs(nBig/2,:))
    xline(nBig/2)

end

%% Compute symmetry metric

range = ceil(nBig/2)+1:nBig+floor(nBig/2)+1;
t2 = NaN(nBig,nBig); % for debugging/understanding selection
for n = 1:length(range)
    t = false(nBig,nBig);
    for i = 1:nBig
        for j = 1:nBig
            if i+j == range(n) && abs(i-j) < (floor(nBig/2)+1)
                t(i,j) = 1;
                t2(i,j) = 1;
            end
        end
    end
    halfLength = floor(sum(t(:))/2);
    sm_t = sm_cs(t); % switch out for desired sm (centred on category)
    sym(n) = mean(sm_t(1:halfLength)) - mean(sm_t(end:-1:end-halfLength+1));

    % bootstrap
    for boot = 1:1000
        switch_ = logical(randi([0, 1], [1, halfLength]));
        if mod(sum(t(:)),2)
            switch_ = [switch_,false,flip(switch_)];
        else
            switch_ = [switch_,flip(switch_)];
        end
        sm_t_bs(switch_) = flip(sm_t(switch_)); % similarity matrix, temporary, bootstrap
        sm_t_bs(~switch_) = sm_t(~switch_);
        sym_bs(boot) = mean(sm_t_bs(1:halfLength)) - mean(sm_t(end:-1:end-halfLength+1)); %symmetry, bootstrap
    end
    lowerCI(n) = prctile(sym_bs,2.5);
    upperCI(n) = prctile(sym_bs,97.5);
end

figure, hold on
imagesc(sm_cs,'AlphaData',t2)
ax = gca();
ax.YDir = 'reverse';
plot([0,nBig],[0,nBig],'r')
axis square tight
colormap('gray')
colorbar

figure, hold on
plot(sym,'k','DisplayName','Symmetry')
%yline(0,'HandleVisibility','off')
axis tight
ylabel('Symmetry (one side - other side)')

plot(lowerCI,'DisplayName','lowerCI')
plot(upperCI,'DisplayName','upperCI')

legend
