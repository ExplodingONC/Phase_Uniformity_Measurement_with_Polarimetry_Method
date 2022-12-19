function [point_result] = voltage2Phase(point_result)
%VOLTAGE2PHASE 此处显示有关此函数的摘要
%   此处显示详细说明

    point_result.result_voltage_waveform_timeavg_kalman_norm = interp1( [min(min(point_result.result_voltage_waveform_timeavg_kalman)) max(max(point_result.result_voltage_waveform_timeavg_kalman))], [0 1], point_result.result_voltage_waveform_timeavg_kalman );
    point_result.result_avg_voltage_norm = mean(point_result.result_voltage_waveform_timeavg_kalman_norm);
    temp.result_voltage_waveform = point_result.result_voltage_waveform_timeavg_kalman_norm;
    point_result.result_avg_phase = zeros(1,256);
    point_result.result_phase_waveform = acos( 1 - temp.result_voltage_waveform * 2 );
    pos = point_result.pos;   % for debug
    
    % find peaks segments from the data to unwrap phase
    point_result.result_avg_phase = mean(point_result.result_phase_waveform);
    [~,max_seg] = findpeaks( point_result.result_avg_phase,'MinPeakDistance',150,'MinPeakHeight', 0.8*pi, 'MinPeakWidth', 1 );
    [~,min_seg] = findpeaks(-point_result.result_avg_phase,'MinPeakDistance',150,'MinPeakHeight',-0.2*pi, 'MinPeakWidth', 1 );
    point_result.extremum_seg = sort( [ max_seg, min_seg ] );

    % if no peaks found
    if isempty(min_seg)
        [~,min_seg] = min(point_result.result_avg_phase,[],'all');
    end
    if isempty(max_seg)
        [~,max_seg] = max(point_result.result_avg_phase,[],'all');
    end
    % and also include starting and ending point
    extended_min_seg = min_seg;
    extended_max_seg = max_seg;
    if ( min_seg(1) > 10 )
        if ( max_seg(1) < min_seg(1) ) && ( point_result.result_avg_phase(0+1) <= point_result.result_avg_phase(max_seg(1)) )
            extended_min_seg = [0+1 extended_min_seg];
        end
    end
    if ( max_seg(1) > 10 )
        if ( min_seg(1) < max_seg(1) ) && ( point_result.result_avg_phase(0+1) >= point_result.result_avg_phase(min_seg(1)) )
            extended_max_seg = [0+1 extended_max_seg];
        end
    end
    if ( min_seg(length(min_seg)) < 247 )
        if ( max_seg(length(max_seg)) > min_seg(length(min_seg)) ) && ( point_result.result_avg_phase(255+1) <= point_result.result_avg_phase(max_seg(length(max_seg))) )
            extended_min_seg = [extended_min_seg 255+1];
        end
    end
    if ( max_seg(length(max_seg)) < 247 )
        if ( min_seg(length(min_seg)) > max_seg(length(max_seg)) ) && ( point_result.result_avg_phase(255+1) >= point_result.result_avg_phase(min_seg(length(min_seg))) )
            extended_max_seg = [extended_max_seg 255+1];
        end
    end
    min_seg = extended_min_seg;
    max_seg = extended_max_seg;
    max_count = length(max_seg);
    min_count = length(min_seg);
    point_result.min_seg = min_seg;
    point_result.max_seg = max_seg;
    
    % unwrapping
    point_result.result_phase_waveform = zeros(size(point_result.result_phase_waveform,1),256);
    if min_seg(1) < max_seg(1)
        % fisrt low peak
        for peak_count = 1 : 1
            point_result.result_phase_waveform(:,min_seg(peak_count)) = 2*pi*(peak_count-1) + acos( 1 - temp.result_voltage_waveform(:,min_seg(peak_count)) * 2 );   % unwrapping
        end
        % ascending parts
        for peak_count = 1 : max_count
            for gray_level = min_seg(peak_count)+1:max_seg(peak_count)
                point_result.result_phase_waveform(:,gray_level) = 2*pi*(peak_count-1) + acos( 1 - temp.result_voltage_waveform(:,gray_level) * 2 );   % unwrapping
            end
        end
        % descending parts
        for peak_count = 1 : min_count-1
            for gray_level = max_seg(peak_count)+1:min_seg(peak_count+1)
                point_result.result_phase_waveform(:,gray_level) = 2*pi*(peak_count) - acos( 1 - temp.result_voltage_waveform(:,gray_level) * 2 );   % unwrapping
            end
        end
    else
        % first high peak
        for peak_count = 1 : 1
            point_result.result_phase_waveform(:,max_seg(peak_count)) = 2*pi*(peak_count) - acos( 1 - temp.result_voltage_waveform(:,max_seg(peak_count)) * 2 );   % unwrapping
        end
        % descending parts
        for peak_count = 1 : min_count
            for gray_level = max_seg(peak_count)+1:min_seg(peak_count)
                point_result.result_phase_waveform(:,gray_level) = 2*pi*(peak_count) - acos( 1 - temp.result_voltage_waveform(:,gray_level) * 2 );   % unwrapping
            end
        end
        % ascending parts
        for peak_count = 1 : max_count-1
            for gray_level = min_seg(peak_count)+1:max_seg(peak_count+1)
                point_result.result_phase_waveform(:,gray_level) = 2*pi*(peak_count) + acos( 1 - temp.result_voltage_waveform(:,gray_level) * 2 );   % unwrapping
            end
        end
    end
    % update average phase
    point_result.result_avg_phase = mean(point_result.result_phase_waveform);

end

