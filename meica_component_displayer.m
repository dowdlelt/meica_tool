function meica_component_displayer()

%select the text file that has the timecourses of interest
ctab = spm_select(Inf, '.*ctab.txt', 'Select ctab file');
ctab_loc = strfind(ctab,'ctab.txt');

ctab_loc =ctab_loc -1;
timecourses = strcat(ctab(1:ctab_loc),'mmix.1D');


fid = fopen(ctab);
tline = fgetl(fid);

while ischar(tline)
    %disp(tline)
    if strfind(tline, '#ACC')
        accp_list = tline;
    elseif strfind(tline, '#REJ')
        rej_list = tline;
    elseif strfind(tline, '#MID')
        mid_list = tline;
    elseif strfind(tline, '#IGN')
        ign_list = tline;
    end
    tline = fgetl(fid);
end
fclose(fid);

accp_list = accp_list(6:end);
rej_list = rej_list(6:end);
mid_list = mid_list(6:end);
ign_list = ign_list(6:end);

accp_list = accp_list(1:(find(accp_list=='#')-1));
rej_list = rej_list(1:(find(rej_list=='#')-1));
mid_list = mid_list(1:(find(mid_list=='#')-1));
ign_list = ign_list(1:(find(ign_list=='#')-1));

accps = str2num(accp_list);
rejs = str2num(rej_list);
mids = str2num(mid_list);
igns = str2num(ign_list);

%Here we will turn these into useful matrices for comparison later.

%% Initialize variables.
filename = ctab;
delimiter = '\t';
startRow = 17;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
imported_ctab = cell2mat(raw);
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

%timecourses = spm_select(Inf, '.*.1D', 'Select ctab file');

%For the imported_ctab, the columns are
%Comp#  Kappa   Rho   Variance   Normed_variance

savedir = spm_select(1,'dir','Select the MEICA output folder...');
cd(savedir);

%%Making the motion plots, including framewise displacement
cfg.motionparam = 'motion.1D';
cfg.prepro_suite = 'meica';
cfg.radius = 50;

figure;
grid on; grid minor;
[fwd,~]=bramila_framewiseDisplacement(cfg); %calculate FD using script

x_axis = size(fwd,1); %Get the number of timepoints

raw_motion = load(cfg.motionparam); %get the SPM regressors

%Subplots are used here to keep everything on the same screen.
% The y axes for the 6 motion estimates are set from the min and max of
% those values.
% The FD y axes is set to top out at 3 - as that is our voxel size.
% The purpose is two fold, anything larger than that is worrisome and
% This makes it easy to jump through all the subjects and compare
% quickly.
subplot(3,1,1); plot(raw_motion(:,1:3)); axis([0 x_axis min(min(raw_motion(:,1:3))) max(max(raw_motion(:,1:3)))]);
title('translation'); grid on;
rots = (raw_motion(:,4:6)); %To convert radians to degrees.
subplot(3,1,2); plot(rots); axis([0 x_axis min(min(rots)) max(max(rots))]);
title('rotation'); grid on;
subplot(3,1,3); plot(fwd); title('Framewise Displacement'); axis([0 x_axis 0 3]);
grid on;

%Thought crosses my mind that I could then use this to create extra
%motion regressor things.





mkdir('component_plots');
savedir = [savedir, '/component_plots/'];

%plots are created here.
[~, titl, ~] = fileparts(cfg.motionparam);
print([savedir, titl], '-dpng');

timecourses_data = load(timecourses);


x_axis = size(timecourses_data,1);
figure;

for i = 1:size(timecourses_data, 2);
    
    kappa = imported_ctab(i,2);
    rho = imported_ctab(i,3);
    variance_explained = imported_ctab(i,4);
    
    %plot(timecourses_data(:,i));
    %% Prints the graphs so that the outcome can be seen
    if any(accps == (i-1))
        plot(timecourses_data(:,i), 'Color',[0 .5 0]); %Green for BOLD like
    elseif any(rejs == (i-1))
        plot(timecourses_data(:,i), 'r'); %Red for Rejected non BOLD
    elseif any(mids == (i-1))
        plot(timecourses_data(:,i), 'm'); %Magenta for R2* weighted artifacts
    elseif any(igns == (i-1))
        plot(timecourses_data(:,i), 'k'); %Black for Ignored components
    end
    %%
    axis([0 x_axis min(timecourses_data(:,i)) max(timecourses_data(:,i))]);
    
    title(strcat('Component:', num2str(i), ', on ctab: ', num2str(i-1), ', kappa: ', num2str(kappa,3), ', rho: ', num2str(rho,3), ', variance: ', num2str(variance_explained,4))); grid on;
    label = strcat('Component_', num2str(i), '_on_ctab_', num2str(i-1));
    print([savedir, label], '-dpng');
end

BOLD_var = 0;
REJ_var = 0;
MID_var = 0;
IGN_var = 0;

for i = 1:size(accps,2)
    BOLD_var = BOLD_var+ imported_ctab(accps(i)+1,5);
end

for i = 1:size(rejs,2)
    REJ_var = REJ_var+ imported_ctab(rejs(i)+1,5);
end

for i = 1:size(mids,2)
    MID_var = MID_var+ imported_ctab((mids(i)+1),5);
end

for i = 1:size(igns,2)
    IGN_var = IGN_var + imported_ctab(igns(i)+1,5);
end

y = [BOLD_var, REJ_var, MID_var, IGN_var];

bar(1, y(1),  'facecolor', [0 .5 0]);
hold on
bar(2, y(2),  'facecolor', 'r');
bar(3, y(3),  'facecolor', 'm');
bar(4, y(4),  'facecolor', 'k');
labels = {'BOLD','non-BOLD','R2* Weighted','Ignored'};
set(gca, 'XTick', 1:4, 'XTickLabel', labels);

%total_var = sum(y);
title('Total variance explained in each category');
ylabel('Variance Explained, %');

print([savedir, 'Var_exp'], '-dpng');






