
function  analysis_rat_for_vmhdm(stim_info, stim_time, spike_time, unit_id, output_name )


% stim_info = find_stim_file( stim_info );
% load(stim_info,'seq','frame_rate','repeat','template','sync_frames') 

repeat=5;
%stim_time=stim_time(1:2:end);


resp_window = [0 10]; %0.25



% stim_frame = 300;%size(stimParams.Templates{1},3);
% blank_frame = 150;

binsize = 1;%s     0.04

seq=[1;2;1;2;1;2;1;2;1;2];
template=[1;2];%length(template);
ntemplate=2;

n_unit=length(unit_id);





% if isfield(stimParams,'sync_frames') && stimParams.sync_frames > 1
%     stim1_duration = mean(diff(stim_time));
% else
%     stim1_duration = stim_time(stim_frame+1) - stim_time(1);
% end



% For FlashingSquare with blank before stim

% % if ~isempty(strfind(stim_info,'FlashingSquare'))
%     blank_duration = mean(diff(stim_time))/3;
%     stim_time = stim_time + blank_duration;
%     stim1_duration = mean(diff(stim_time))/3*2;
% % end

stim1_duration = 10;
blank_duration = 0;


% pre_stim_time = - stim1_duration;
% post_stim_time = stim1_duration +stim1_duration;

pre_stim_time = - 15; % second
post_stim_time = 120;

%% cut segment spike data relevent to this stimuli

start_time = stim_time(1)+pre_stim_time;
end_time = stim_time(end)+post_stim_time;

for i = 1:n_unit
    included_spike = spike_time{i}>start_time & spike_time{i}<end_time;
    spike_time{i} = spike_time{i}(included_spike);
end

total_spikes = zeros(n_unit,1);
for i = 1:n_unit
    total_spikes(i) = length(spike_time{i});
end

%% 

% find stim start time
nframe = length(stim_time);

nseq = length(seq);

%nstim =2;

stim_start_frame=1:nframe;


% if isfield(stimParams,'sync_frames') && stimParams.sync_frames > 1
%     stim_start_frame = 1:nframe;
%     
% else
%     
%     stim_start_frame = zeros(nstim*repeat,1); 
%     
%     k = 0;
%     last_stim_end_frame = 0;
%     
%     for i = 1:nseq
%         stim_id = seq(i);
%         k = k+1;
%         stim_start_frame(k) = last_stim_end_frame + blank_frame + 1;
%         
%         last_stim_end_frame = stim_start_frame(k) + stim_frame -1;
%         
%     end
% 
% end


%spike train
% tic
% stim_spike_time=cell(n_unit, nseq);
% for i = 1:n_unit
%     for j=1:nseq
%     included_spike = spike_time{i}>stim_time(stim_start_frame(j))-each_frame_time & spike_time{i}<stim_time(stim_end_frame(j))+each_frame_time;
%     stim_spike_time{i, j}= spike_time{i}(included_spike);
%     end
% end

%spike train
stim_spike_time=cell(n_unit, nseq);
for i = 1:n_unit
    for j=1:nseq
        spt = spike_time{i} - stim_time(stim_start_frame(j));
        included_spike = spt > pre_stim_time & spt < post_stim_time;
        stim_spike_time{i, j}= spt(included_spike);
    end
end

% organized spike time
stim_spike_time_rep=cell(n_unit, ntemplate, repeat);
for i=1:n_unit
    for j=1:ntemplate
        ind=find(seq==j);
        for jj=1:length(ind)
        stim_spike_time_rep{i,j,jj}=stim_spike_time{i, ind(jj)};
        end
    end
end

% psth
bin_edge = pre_stim_time:binsize:post_stim_time;
bin_x = bin_edge(1:end-1)+binsize/2;
nbin = length(bin_x);
psth_all = zeros(nbin, n_unit, ntemplate, repeat);

for i = 1:n_unit
    for j = 1:ntemplate
        for k = 1:repeat
            psth_all(:,i,j,k) = histcounts( stim_spike_time_rep{i,j,k}, bin_edge)/binsize;
        end
    end
end

% spike cout within specified window
spike_count = zeros(n_unit, ntemplate, repeat);
for i=1:n_unit
    for j=1:ntemplate
        for k=1:repeat
            spike_count(i,j,k) = sum( ...
                stim_spike_time_rep{i,j,k} > resp_window(1) & ...
                stim_spike_time_rep{i,j,k} < resp_window(2) ); 
        end
    end
end

firing_rate_stim = spike_count/(resp_window(2)-resp_window(1));  
%firing_rate_stim = spike_count/stim1_duration;     %  response_window

% compute anova
% p = zeros(3, 5, n_unit);
% 
% for i = 1:5
%     for j = 1:n_unit
%         ind = i*2-1;
%         ind_row1 = [stim_organization(1,ind) stim_organization(1,ind+1)];
%         ind_row2 = [stim_organization(2,ind) stim_organization(2,ind+1)];
%         y = [ squeeze(spike_count(j, ind_row1, :))'; squeeze(spike_count(j, ind_row2, :))' ];
%         p (:,i,j) = anova2(y,repeat, 'off');
%     end
% end


%% save data
rat.unit_id = unit_id;
%MS.stim_params = stimParams;
rat.total_spikes = total_spikes;
rat.p = [];
rat.p_window = resp_window;

rat.pre_stim_time = pre_stim_time;
rat.post_stim_time = post_stim_time;
rat.psth_all = psth_all;
rat.bin_x = bin_x;
rat.stim_spike_time_rep = stim_spike_time_rep;
rat.blank_duration = blank_duration;
rat.stim1_duration = stim1_duration;
rat.firing_rate_stim = firing_rate_stim;

save(output_name, 'rat')


%% figure parameters

[~, figure_name_str] = fileparts(stim_info);

figure_name_str='rat';

%% figure
if 1
    
    % stim organization
    %stim_organization = reshape(1:ntemplate, 4, 8); % 32 position
    stim_organization = reshape(1:ntemplate, 2, 1); %128 position
    [row, col] = size(stim_organization);
    
    % show stim organization
%     figure;
%     colormap gray
%     for i = 1:row
%         for j = 1:col
%             
%             stim_ind = stim_organization(i,j);
%             subplot(row, col, (i-1)*col+j)
%             imagesc(template{stim_ind}(:,:,1), [0 255]); 
%             axis off
%             axis equal
%             
%         end
%     end
    %figure_raster_psth_overlay_1( rat, stim_organization, figure_name_str )
    figure_raster_psth( rat, stim_organization, figure_name_str )

end


%% figure
if 0
    
    % stim organization
     %stim_organization = reshape(1:ntemplate, 4, 8); % 32 position
     %ntemplate =128; %%gg
     stim_organization = reshape(1: ntemplate, 4*2, 8*2); %128 position
     [row, col] = size(stim_organization);
    
    % firing_rate_stim=MS.firing_rate_stim;  %%gg
     
     
    for i = 1:n_unit    %74  %62
    
    % raster and psth
    figure('position', [1 41 700 400]);
    
    fr = mean(firing_rate_stim(i,:,:),3); % average across repeat
    
    fr = fr - mean(fr); % subtract mean
    grid = reshape(fr,[row, col]);
    imagesc(grid);
    axis off; axis equal
    
   % suptitle([ 'unit ID ' num2str(unit_id(unit_included(i)))] );
   
   %gg
    suptitle([ 'unit ID ' num2str(unit_id(i))] );
    
    %fig_filename = fullfile('fig', sprintf('ID%03d_%s.bmp', unit_id(unit_included(i)) , figure_name_str ) );
   
    %ggg
    
    fig_filename = fullfile('fig', sprintf('ID%03d_%s.bmp', unit_id(i) , figure_name_str ) );
    savefig_good(gcf, fig_filename)
    close(gcf)
    
    end
    

    
end
    
    % show stim organization
%     figure;
%     colormap gray
%     for i = 1:row
%         for j = 1:col
%             
%             stim_ind = stim_organization(i,j);
%             subplot(row, col, (i-1)*col+j)
%             imagesc(template{stim_ind}(:,:,1), [0 255]); 
%             axis off
%             axis equal
%             
%         end
%     end
% 




end 









