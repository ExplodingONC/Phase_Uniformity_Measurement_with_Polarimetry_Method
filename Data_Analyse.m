
%% start
    
    clearvars -except f_*;
    clc;
    all_fig = findall(0, 'type', 'figure');
    close(all_fig);
    clearvars all_fig;
    start = tic;
    tic;
    
%% load configuration and data

    result_from_calibration = false;
    generate_cal_profile = false;
    load('result_520nm_1.3-3.0_2x_20220611-043433.mat');
    step_cnt_X = length(step_grid_X);
    step_cnt_Y = length(step_grid_Y);
    step_size_X = ( step_grid_X(step_cnt_X) - step_grid_X(1) ) / ( step_cnt_X-1 ) * 1e-3;
    step_size_Y = ( step_grid_X(step_cnt_X) - step_grid_X(1) ) / ( step_cnt_X-1 ) * 1e-3;
    step_size = ( step_size_X + step_size_Y ) / 2;
    fprintf(' -- Data loading finished. (%6.3fs)\n',toc); tic;
    
    if result_from_calibration == true
        generate_cal_profile = false;
    end
    
%% result calculation
    
    % intensity_parallel = I0/2 * ( 1 + cos(phase_ret) )
    % intensity_perpendicular = I0/2 * ( 1 - cos(phase_ret) )
    
    % result_intensity = result_voltage / 4605;      % for 517.8nm laser,
    % 20db gain, High-Z 
    
    result_group_stitched = cell(step_cnt_Y,step_cnt_X);
    for y = 1:step_cnt_Y
        for x = 1:step_cnt_X
            % read out
            point_result_original = result_group_original{y,x};
            point_result_retarded = result_group_retarded{y,x};
%             if ~isfield(point_result_original,'pos') || ~isfield(point_result_retarded,'pos')
%                 fprintf( ' -! Point (%1d,%1d) empty\n', x, y );
%                 continue;   % skip if empty
%             end
            % phase calc
            try
                point_result_original = voltage2Phase(point_result_original);
                % save back (update)
                result_group_original{y,x} = point_result_original;
            catch
                fprintf( ' -! Failed calculation on original point (%1d,%1d) !\n', x, y );
            end
            try
                point_result_retarded = voltage2Phase(point_result_retarded);
                % save back (update)
                result_group_retarded{y,x} = point_result_retarded;
            catch
                fprintf( ' -! Failed calculation on retarded point (%1d,%1d) !\n', x, y );
            end
            % stitching
            try
                point_result_stitched = phaseStitching( point_result_original, point_result_retarded );
                result_group_stitched{y,x} = point_result_stitched;
            catch
                fprintf( ' -! Failed stitching on point (%1d,%1d) !\n', x, y );
                result_group_stitched{y,x} = result_group_original{y,x};
            end
        end
    end
    clearvars x y point_result_*;
    fprintf(' -- Phase calculation finished. (%6.3fs)\n',toc); tic;

%% phase fit (maybe change to iteration method?)

    phase_sample = zeros( step_cnt_Y, step_cnt_X, 256 );
    phase_profile = zeros( res_Y, res_X, 256 );
    
    % scan vertical sections
    for j = 1:step_cnt_Y
        % determine vertical range
        if j == 1
            border_top = 1;
        else
            border_top = fix( ( ( step_grid_Y(j-1) + step_grid_Y(j) ) / 2 ) * 1e-3 / pixel_pitch + res_Y / 2 ) + 1;
        end
        if j == step_cnt_Y
            border_bottom = res_Y;
        else
            border_bottom = fix( ( ( step_grid_Y(j) + step_grid_Y(j+1) ) / 2 ) * 1e-3 / pixel_pitch + res_Y / 2 );
        end
        % scan horizontal sections
        for i = 1:step_cnt_X
            % determine horizontal range
            if i == 1
                border_left = 1;
            else
                border_left = fix( ( ( step_grid_X(i-1) + step_grid_X(i) ) / 2 ) * 1e-3 / pixel_pitch + res_X / 2 ) + 1;
            end
            if i == step_cnt_X
                border_right = res_X;
            else
                border_right = fix( ( ( step_grid_X(i) + step_grid_X(i+1) ) / 2 ) * 1e-3 / pixel_pitch + res_X / 2 );
            end
            % fill corresponding section
            phase_sample( j, i, : ) = result_group_original{j,i}.result_avg_phase;  % maybe smooth out the curve here?
            for actual_gray_level = 0:255
                phase_profile( border_top:border_bottom, border_left:border_right, actual_gray_level+1 ) = phase_sample( j, i, actual_gray_level+1 );
            end
        end
    end
    clearvars i j border_*;
    
    % low-pass filter
    filter_size = (step_size / pixel_pitch) / (exp(1));
    phase_profile = imgaussfilt( phase_profile, filter_size );

    % compensation
    for actual_gray_level = 0:255
%         phase_profile(:,:,actual_gray_level+1) = interp1( [min(phase_sample(:,:,actual_gray_level+1),[],'all') max(phase_sample(:,:,actual_gray_level+1),[],'all')], [0 1], phase_profile(:,:,actual_gray_level+1), "linear");
%         phase_profile(:,:,actual_gray_level+1) = imadjust( phase_profile(:,:,actual_gray_level+1), [min(phase_profile(:,:,actual_gray_level+1),[],'all') max(phase_profile(:,:,actual_gray_level+1),[],'all')] );
%         phase_profile(:,:,actual_gray_level+1) = interp1( [0 1], [min(phase_sample(:,:,actual_gray_level+1),[],'all') max(phase_sample(:,:,actual_gray_level+1),[],'all')], phase_profile(:,:,actual_gray_level+1), "linear");
    end
    clearvars actual_gray_level;
    fprintf(' -- Full phase profile generated. (%6.3fs)\n',toc); tic;
    
%% calibration profile generate

    % skip if already from calibrated display
    if result_from_calibration == false && generate_cal_profile == true
        
        % data storage
        cal_profile = zeros( step_cnt_Y, step_cnt_X, 256 );
        cal_img = uint8(zeros(res_Y,res_X,256));
        dist = zeros(step_cnt_Y,step_cnt_X);
        
        % target phase determination
        phase_start = -2*pi;
        phase_end = 4*pi;
        for y = 1:step_cnt_Y
            for x = 1:step_cnt_X
                if result_group_stitched{y,x}.result_avg_phase(1) > phase_start
                    phase_start = result_group_stitched{y,x}.result_avg_phase(1);
                end
                if result_group_stitched{y,x}.result_avg_phase(256) < phase_end
                    phase_end = result_group_stitched{y,x}.result_avg_phase(256);
                end
            end
        end
        phase_start = (phase_start + phase_end) / 2 - pi;
        phase_end = phase_start + 2*pi;
        
        % calibration profile
        for y = 1:step_cnt_Y
            for x = 1:step_cnt_X
                for target_gray_level = 0:255
                    % required phase output
                    target_phase = 2*pi / 256 * target_gray_level + 0.2*pi;
                    % find required GL input
                    for actual_gray_level = 0:255
                        if result_group_stitched{y,x}.result_avg_phase(actual_gray_level+1) >= target_phase
                            cal_profile(y,x,target_gray_level+1) = actual_gray_level;
                            break;
                        end
                    end
                    if actual_gray_level == 255
                        cal_profile(y,x,target_gray_level+1) = 255;
                    end
                end
            end
        end
        clearvars x y target_phase target_gray_level actual_gray_level;
        fprintf(' -- Calibration profile finished. (%6.3fs)\n',toc); tic;
        
        % generate full-size profile
        % this is faster than interp2+fillmissing with 'nearest' method
        for y = 1:res_Y
            for x = 1:res_X
                dist = zeros( size(cal_profile,[1 2]) );
                for j = 1:step_cnt_Y
                    for i = 1:step_cnt_X
                        dist(j,i) = fix(abs( x - step_grid_X(i)*1e-3/pixel_pitch - res_X/2 )) + fix(abs( y - step_grid_Y(j)*1e-3/pixel_pitch - res_Y/2 ));
                    end
                end
                [~,index] = min( dist, [],'all','linear' );
                [j,i] = ind2sub( size(cal_profile,[1 2]) , index );
                cal_img(y,x,:) = cal_profile(j,i,:);
            end
        end
        fprintf(' -- Calibration image (LUT) generated. (%6.3fs)\n',toc); tic;
        clearvars x y i j dist index;

        % low-pass filter
        filter_size = (step_size / pixel_pitch) / (exp(1));
        %filter_h = fspecial( 'disk', filter_size );    % 'gaussian','disk','laplacian'
        %cal_img = imfilter( cal_img, filter_h, 'replicate', 'conv' );
        cal_img = imgaussfilt( cal_img, filter_size );
        fprintf(' -- Calibration image (LUT) filtered. (%6.3fs)\n',toc); tic;
        
        % compensate for contrast lost
        cal_img = single(cal_img);
        for target_gray_level = 0:255
            %current_min = min(min(cal_img( :, :, target_gray_level+1 )));
            %current_mean = mean(cal_img( :, :, target_gray_level+1 ), 'all');
            %current_max = max(max(cal_img( :, :, target_gray_level+1 )));
            %target_min = min(min(cal_profile( :, :, target_gray_level+1 )));
            %target_mean = mean(cal_profile( :, :, target_gray_level+1 ), 'all');
            %target_max = max(max(cal_profile( :, :, target_gray_level+1 )));
            %cal_img(:,:,target_gray_level+1) = interp1( [0.99*current_min current_mean 1.01*current_max], [0.95*target_min target_mean 1.05*target_max], cal_img(:,:,target_gray_level+1), 'makima' );
        end
        clearvars target_gray_level local_min local_max;
        cal_img = max(0,cal_img); cal_img = min(255,cal_img);
        
        cal_img = uint8(cal_img);
        fprintf(' -- Calibration image (LUT) compensated. (%6.3fs)\n',toc); tic;
        
    end
    
%% result display
    
    all_fig = findall(0, 'type', 'figure');
    close(all_fig);
    clearvars all_fig;
    % display enable
    display_phase_waveform = false;
    display_single_point = false;
    display_multi_point = true;
    display_phase_distribution = true;
    display_calibration = false;
    % display option
    uniformity_GL_base = 250;
    disp_point_index = sub2ind(size(result_group_stitched),ceil(size(result_group_stitched,1)/2),ceil(size(result_group_stitched,2)/2));    % center
    % disp_point_index = sub2ind(size(result_group_stitched),2,6);    % designated point
    
    % phase data comparison
    if display_phase_waveform
        f_comparison = figure('Name','Phase Waveform','NumberTitle','off','WindowState','Maximized');
        a_result_1 = gca;
        data_point = 0:10:80;
        for i = 1:9
            a_comparison(i) = subplot(3,3,i);
            plot( 0:size(result_group_stitched{disp_point_index}.result_phase_waveform,1)-1, result_group_stitched{disp_point_index}.result_phase_waveform(:,data_point(i)+1)/pi);
            xlim( [0 size(result_group_stitched{disp_point_index}.result_phase_waveform,1)-1] );
            xlabel(a_comparison(i), 'Time (us)');
            ylabel(a_comparison(i), 'Phase (\pi)');
            title(a_comparison(i), [ 'Graylevel: ' num2str(data_point(i)) ]);
        end
        clearvars i data_point a_comparison temp;
    end

    % single point data display
    if display_single_point
        f_result = figure('Name','Result','NumberTitle','off','WindowState','Maximized');
        a_result_2 = gca;
        title(a_result_2, [ num2str(1e9*wavelength) 'nm - phase retardation and flicker depth @pos [' num2str(result_group_original{disp_point_index}.pos) ']' ]);
        xlabel(a_result_2, 'Gray level');
        xlim( [0 255] );
        yyaxis left
        hold on
        plot( 0:255, (result_group_original{disp_point_index}.result_avg_phase)/pi, '-c', ...
            0:255, (result_group_retarded{disp_point_index}.result_avg_phase-result_group_stitched{disp_point_index}.phase_retardance)/pi, '-g', 'LineWidth', 0.5 );
        plot( 0:255, (result_group_stitched{disp_point_index}.result_avg_phase)/pi, '-b', 'LineWidth', 1.5 );
        plot( result_group_stitched{disp_point_index}.extremum_from_original-1, result_group_original{disp_point_index}.result_avg_phase(result_group_stitched{disp_point_index}.extremum_from_original)/pi, 'oc', 'LineWidth', 1.5 );
        plot( result_group_stitched{disp_point_index}.extremum_from_retarded-1, ( result_group_retarded{disp_point_index}.result_avg_phase(result_group_stitched{disp_point_index}.extremum_from_retarded) - result_group_stitched{disp_point_index}.phase_retardance )/pi, 'og', 'LineWidth', 1.5 );
        hold off
        ylim( [ min(result_group_stitched{disp_point_index}.result_avg_phase)/pi, 1.1*max(result_group_stitched{disp_point_index}.result_avg_phase)/pi ] );
        ylabel(a_result_2, 'Average phase retardation (\pi)');
        yyaxis right
        plot( 0:255, (max(result_group_stitched{disp_point_index}.result_phase_waveform)-min(result_group_stitched{disp_point_index}.result_phase_waveform))/pi, '--r', ...
            0:255, (std(result_group_stitched{disp_point_index}.result_phase_waveform))/pi, ':r', 'LineWidth', 1 );
        ylim( [ 0, inf ] );
        ylabel(a_result_2, 'Phase flicker (\pi)');
        grid(a_result_2, 'on');
        %legend(a_result_2, {'avg phase (original)','avg phase (QWP addon)','avg phase (stitched)','flicker (p-p)','flicker (std)'}, 'Location','northwest');
    end

    % phase uniformity
    if display_multi_point
        f_uniformity_1 = figure('Name','Phase Uniformity (compare)','NumberTitle','off','WindowState','Maximized');
        a_result_3 = gca;
        label = [];
        hold on
        for disp_point_index = 1:(step_cnt_X*step_cnt_Y)
            plot( 0:255, (result_group_original{disp_point_index}.result_avg_phase)/pi, 'LineWidth', 1 );
            label = [ label, string(mat2str(result_group_original{disp_point_index}.pos)) ];
        end
        clearvars disp_point_index;
        hold off
        title(a_result_3, [ num2str(1e9*wavelength) 'nm - phase retardation on different areas' ]);
        xlabel(a_result_3, 'Gray level'); xlim( [0 255] );
        ylabel(a_result_3, 'Average phase retardation (\pi)'); %ylim( [ 0, 3 ] );
        grid(a_result_3, 'on');
        legend(a_result_3, label, 'Location','southeast');

        f_intensity_1 = figure('Name','Intensity Uniformity (compare)','NumberTitle','off','WindowState','Maximized');
        a_result_3_2 = gca;
        label = [];
        hold on
        for disp_point_index = 1:(step_cnt_X*step_cnt_Y)
            plot( 0:255, (result_group_original{disp_point_index}.result_avg_voltage_norm), 'LineWidth', 1 );
            label = [ label, string(mat2str(result_group_original{disp_point_index}.pos)) ];
        end
        clearvars disp_point_index;
        hold off
        title(a_result_3_2, [ num2str(1e9*wavelength) 'nm - received intensity on different areas' ]);
        xlabel(a_result_3_2, 'Gray level'); xlim( [0 255] );
        ylabel(a_result_3_2, 'Average Light Intensity'); %ylim( [ 0, 3 ] );
        grid(a_result_3_2, 'on');
        legend(a_result_3_2, label, 'Location','southeast');
    end
    
    % phase distribution
    if display_phase_distribution
        % pack needed data
        packed_data.phase_profile = phase_profile;
        packed_data.phase_sample = phase_sample;
        packed_data.pixel_pitch = pixel_pitch;
        packed_data.step_grid_X = step_grid_X;
        packed_data.step_grid_Y = step_grid_Y;
        packed_data.res_X = res_X;
        packed_data.res_Y = res_Y;
        packed_data.wavelength = wavelength;
        % create UI figure
        f_uniformity_2 = uifigure('Name','Phase Uniformity (distribution)','NumberTitle','off','WindowState','Maximized');
        gl = uigridlayout(f_uniformity_2);
        gl.RowHeight = {25,'1x','fit',25};
        gl.ColumnWidth = {'1x','18x','1x'};
        a_result_4 = uiaxes(gl);
        a_result_4.Layout.Row = 2;
        a_result_4.Layout.Column = 2;
        a_result_4.View = [0 90];
        slider = uislider(gl);
        slider.Layout.Row = 3;
        slider.Layout.Column = 2;
        slider.Tooltip = 'Displayed Graylevel';
        slider.Limits = [0 255];
        slider.MajorTicks = 0:5:255;
        slider.MinorTicks = 0:1:255;
        slider.Value = 127;
        slider.ValueChangingFcn = @(slider,event) updatePhaseFigure(slider,event,a_result_4,packed_data);
        slider.ValueChangedFcn = @(slider,event) updatePhaseFigure(slider,event,a_result_4,packed_data);
        % initial uptade
        align([a_result_4 slider],'Center','none');
        updatePhaseFigure(slider,[],a_result_4,packed_data);
        clearvars packed_data;
    end
    
    % calibration profile
    if result_from_calibration == false && generate_cal_profile == true && display_calibration == true
        f_calibration = figure('Name','Calibration Profile','NumberTitle','off','WindowState','Maximized');
        a_result_5_1 = subplot(1,2,1);
        imshow( cal_profile(:,:,uniformity_GL_base+1), [0 255] );
        a_result_5_2 = subplot(1,2,2);
        imshow( cal_img(:,:,uniformity_GL_base+1), [0 255] );
    end

    fprintf(' -- Results displayed. (%6.3fs)\n',toc); tic;
    clearvars a_result_* f_*;
    
%% end
    
    fprintf( '\n' );
    toc(start);
    
%% functions

function updatePhaseFigure(slider, event, axes, packed_data)
    % unpack data
    phase_profile = packed_data.phase_profile;
    phase_sample = packed_data.phase_sample;
    pixel_pitch = packed_data.pixel_pitch;
    step_grid_X = packed_data.step_grid_X;
    step_grid_Y = packed_data.step_grid_Y;
    res_X = packed_data.res_X;
    res_Y = packed_data.res_Y;
    wavelength = packed_data.wavelength;
    % extract value (graylevel)
    if ~isempty(event)
        Value = round(event.Value);
        event.Source.Value = Value;
    else
        Value = round(slider.Value);
        slider.Value = Value;
    end
    % update figure
    point_X = fix(linspace(1,res_X,res_X/10));
    point_Y = fix(linspace(1,res_Y,res_Y/10));
    View = axes.View;
%     imagesc(axes, phase_profile(:,:,Value+1)/pi);
%     colorbar(axes);
    surfc(axes, point_X, point_Y, phase_profile(point_Y,point_X,Value+1)/pi, 'EdgeColor','none','FaceAlpha',0.9);
    hold(axes, 'on');
    stem3(axes, round(step_grid_X*1e-3/pixel_pitch+res_X/2), round(step_grid_Y*1e-3/pixel_pitch+res_Y/2), phase_sample(:,:,Value+1)/pi, 'filled', 'o', 'LineStyle','none');
    hold(axes, 'off');
%     b = bar3(axes,phase_sample(:,:,Value+1)/pi);
%     for k = 1:length(b)
%         zdata = b(k).ZData;
%         b(k).CData = zdata;
%         b(k).FaceColor = 'interp';
%     end
    axes.View = View;
    axes.DataAspectRatio = [1 1 0.005];
    axis(axes, 'ij');
    xlim(axes, [1 res_X]);
    ylim(axes, [1 res_Y]);
    zlim(axes, [ min(phase_profile,[],'all')/pi, max(phase_profile,[],'all')/pi ]);
    title(axes, [ num2str(1e9*wavelength) 'nm - phase uniformity @GL' num2str(Value) ]);
%     title(axes, [ num2str(1e9*wavelength) 'nm - phase uniformity' ]);
    xlabel(axes, 'X');
    ylabel(axes, 'Y');
    zlabel(axes, 'Average phase retardation (\pi)');
    grid(axes, 'on');
end
    