
function  analysis_laser_ultrasound(stim_time, stim_duration, spike_time, unit_id, output_name, figure_name_str )


binsize = 0.005;%s

n_unit=length(unit_id);

pre_stim_time = -0.02;
post_stim_time = 0.05;

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
nstim = length(stim_time);

%spike train
stim_spike_time = cell(n_unit, nstim);
for i = 1:n_unit
    for j=1:nstim
        spt = spike_time{i} - stim_time(j);
        included_spike = spt > pre_stim_time & spt < post_stim_time;
        stim_spike_time{i, j}= spt(included_spike);
    end
end

% psth
bin_edge = pre_stim_time:binsize:post_stim_time;
bin_x = bin_edge(1:end-1)+binsize/2;
nbin = length(bin_x);
psth_all = zeros(nbin, n_unit, nstim);

for i = 1:n_unit
    for j = 1:nstim
        try
            psth_all(:,i,j) = histcounts( stim_spike_time{i,j}, bin_edge)/binsize;
        catch
        end
    end
end


%% save data
O.unit_id = unit_id;
O.total_spikes = total_spikes;

O.pre_stim_time = pre_stim_time;
O.post_stim_time = post_stim_time;
O.psth_all = psth_all;
O.bin_x = bin_x;
O.stim_spike_time = stim_spike_time;
O.stim_duration = stim_duration;

save(output_name, 'O')



%% figure
if 1
    
    unit_included=find(O.total_spikes>0);
    
    n_unit = length(unit_included);
    unit_id=O.unit_id;
    repeat = size(O.stim_spike_time,2);
    stim_duration = O.stim_duration;
    bin_x = O.bin_x;
    
    %     unit_valid = false(length(O.unit_id),1);
    %     unit_valid = true()
    
    for i =1:n_unit
        
        % raster and psth
        figure('position', [1 41 800 1000]);
        
        spt_rep = O.stim_spike_time(unit_included(i), :);
        psth = squeeze( O.psth_all(:, unit_included(i), :));
        
        subplot(4, 1, 1:3)
        rectangle('Position',[0 0 stim_duration*1000 repeat], 'FaceColor',[0.7 1 0.7], 'EdgeColor', [0.7 1 0.7]);
        
        hold on
        for r=1:repeat
            rasterplot(spt_rep{r}*1000, r, 0.7, 'k');
        end
        axis([O.pre_stim_time*1000  O.post_stim_time*1000  0.5 repeat+0.5])
%         set(gca,'xticklabel',[])
%         set(gca,'yticklabel',[])

        subplot(4, 1, 4)
        sem = std(psth, 0, 2)/sqrt(repeat);
        ave = mean(psth,2);
        error_area(bin_x*1000, ave, sem, [0.7 0.7 0.7]);
        plot(bin_x*1000, ave,'k')
        xlim([O.pre_stim_time O.post_stim_time]*1000)
        %set(gca,'xticklabel',[])
        %set(gca,'yticklabel',[])
        
        suptitle([ 'unit ID ' num2str(unit_id(unit_included(i)))] );
        
        fig_filename = fullfile('fig', sprintf('ID%03d_%s_raster.bmp', unit_id(unit_included(i)) , figure_name_str ) );
        savefig_good(gcf, fig_filename)
        close(gcf)
        
    end
    
    
end





