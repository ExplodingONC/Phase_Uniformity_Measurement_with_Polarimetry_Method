function [ point_result_stitched ] = phaseStitching( point_result_original, point_result_retarded )
%PHASESTITCHING 此处显示有关此函数的摘要
%   此处显示详细说明
    
    % prepare and data copy
    if point_result_original.pos == point_result_retarded.pos
        point_result_stitched.pos = point_result_original.pos;
        pos = point_result_stitched.pos;    % for debug
    else
        fprintf( 'ERROR: Inputs position does not match!\n' );
        return
    end
    point_result_stitched.extremum_from_original = point_result_original.extremum_seg;
    point_result_stitched.extremum_from_retarded = point_result_retarded.extremum_seg;
    % vars declaration
    point_result_stitched.phase_retardance = mean(point_result_retarded.result_phase_waveform,'all') - mean(point_result_original.result_phase_waveform,'all');
    point_result_stitched.result_avg_temp = zeros(1,256);
    point_result_stitched.result_avg_phase = zeros(1,256);
    point_result_stitched.result_phase_waveform = zeros( size( point_result_original.result_phase_waveform, 1 ), 256 );
    
    % retarded phase data offset
    if point_result_stitched.phase_retardance < 0
        point_result_original.result_phase_waveform = point_result_original.result_phase_waveform - 2*pi;
        point_result_stitched.phase_retardance = mean(point_result_retarded.result_phase_waveform,'all') - mean(point_result_original.result_phase_waveform,'all');
    end
    % unwrap period
    if mean(point_result_original.result_phase_waveform(:,1)) < 0
        point_result_original.result_phase_waveform = point_result_original.result_phase_waveform + 2*pi;
        point_result_retarded.result_phase_waveform = point_result_retarded.result_phase_waveform + 2*pi;
    end
    if mean(point_result_original.result_phase_waveform(:,1)) > 2*pi
        point_result_original.result_phase_waveform = point_result_original.result_phase_waveform - 2*pi;
        point_result_retarded.result_phase_waveform = point_result_retarded.result_phase_waveform - 2*pi;
    end
    % phase shift
    %point_result_stitched.phase_retardance = 0.5*pi;
    point_result_retarded.result_phase_waveform = point_result_retarded.result_phase_waveform - point_result_stitched.phase_retardance;
    
    % stitch the phase data together
    for gray_level = 0:255
        % choose the data closer to half-intensity point
        voltage_original = abs( mean( point_result_original.result_voltage_waveform_timeavg_kalman_norm(:,gray_level+1) ) - 0.5 );
        voltage_retarded = abs( mean( point_result_retarded.result_voltage_waveform_timeavg_kalman_norm(:,gray_level+1) ) - 0.5 );
        if voltage_original < voltage_retarded
            point_result_stitched.result_phase_waveform( : , gray_level+1 ) = point_result_original.result_phase_waveform( : , gray_level+1 );
            point_result_stitched.result_avg_temp( : , gray_level+1 ) = point_result_original.result_avg_temp( : , gray_level+1 );
        else
            point_result_stitched.result_phase_waveform( : , gray_level+1 ) = point_result_retarded.result_phase_waveform( : , gray_level+1 );
            point_result_stitched.result_avg_temp( : , gray_level+1 ) = point_result_retarded.result_avg_temp( : , gray_level+1 );
        end
    end
    % update average phase
    point_result_stitched.result_avg_phase = mean(point_result_stitched.result_phase_waveform);
    temp = [flip(point_result_stitched.result_avg_phase) point_result_stitched.result_avg_phase flip(point_result_stitched.result_avg_phase)];
    temp = lowpass(temp,0.1);
    point_result_stitched.smoothed_avg_phase = temp(257:512);

end

