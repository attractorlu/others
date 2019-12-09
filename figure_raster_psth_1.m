function [  ] = figure_raster_psth_1( data_info, i_unit, stim_organization )
%FIGURE_raster psth
%

[row, col] = size(stim_organization);
repeat = size(data_info.stim_spike_time_rep,3);
stim1_duration = data_info.stim1_duration;
bin_x = data_info.bin_x;
unit_id = data_info.unit_id;

%% figure raster and psth

figure('position', [ 2          96        1918         896 ]);%667   202   983   516   % 10         202        1889         516
max_ave_spikes = max(max(mean(data_info.psth_all(:, i_unit, :, :), 4))) + 1e-6;

for j = 1:row
    for k = 1:col
        
        stim_ind = stim_organization(j,k);
        if isnan(stim_ind)
            continue;
        end
        
        spt_rep = data_info.stim_spike_time_rep(i_unit, stim_ind, :);
        psth = squeeze( data_info.psth_all(:, i_unit, stim_ind, :));
        
        % raster
        subplot(row*2, col, (j-1)*col*2+k)
        %rectangle('Position',[0 0.8 stim1_duration 30], 'FaceColor',[208/255 206/255 206/255], 'EdgeColor', [208/255 206/255 206/255]);   %0.9 0.94 1   %0.61 0.84 0.95
        
        %green
      %rectangle('Position',[0 0.2 stim1_duration 30], 'FaceColor',[0.61 0.95 0.84   ], 'EdgeColor', [0.61 0.95 0.84    ]);   %0.9 0.94 1   %0.61 0.84 0.95
        % blue
        rectangle('Position',[0 0.2 stim1_duration 30], 'FaceColor',[0.61  0.84 0.95   ], 'EdgeColor', [0.61 0.84 0.95    ]);   
        alpha(.8)
        
        
        hold on
        for r=1:repeat
            rasterplot(spt_rep{r}, r, 0.7, 'k');
        end
        axis([data_info.pre_stim_time  data_info.post_stim_time  0.5 repeat+0.5])
        set(gca,'xticklabel',[])
        if k ~= 1
            set(gca,'yticklabel',[])
        end
        %ylabel('Trials')
        
        set(gca, 'LineWidth', 1.3, 'FontSize', 11, 'TickDir','out');
        
        % psth
         subplot(row*2, col, (j-1)*col*2+col+k)
           %green
      %rectangle('Position',[0 0.2 stim1_duration 1000], 'FaceColor',[0.61 0.95 0.84   ], 'EdgeColor', [0.61 0.95 0.84    ]);   %0.9 0.94 1   %0.61 0.84 0.95
        % blue
        rectangle('Position',[0 0.2 stim1_duration 1000], 'FaceColor',[0.61  0.84 0.95   ], 'EdgeColor', [0.61 0.84 0.95    ]);   
    
         
         
 
          hold on
          
        sem = std(psth, 0, 2)/sqrt(repeat);
        ave = mean(psth,2);
        error_area(bin_x, ave, sem, [0.7 0.7 0.7]);
        plot(bin_x, ave,'k')
        ylim([0 max_ave_spikes+5])
        xlim([data_info.pre_stim_time data_info.post_stim_time])
       
       
        
        if j ~= row
            set(gca,'xticklabel',[])
        end
        if k ~= 1
            set(gca,'yticklabel',[])
        end
      %  ylabel('Firing rate')
        set(gca, 'LineWidth', 1.3, 'FontSize', 11, 'TickDir','out');
    end
end

suptitle([ 'unit ID ' num2str(unit_id(i_unit))] );

end







