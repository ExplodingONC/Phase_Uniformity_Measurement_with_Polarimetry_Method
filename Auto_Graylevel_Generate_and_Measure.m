
%% start
    
    clearvars -except stage ps5000aDeviceObj;    % keep variables in case device failed to shut down properly
    clc; close all;
    start_time_stamp = tic;
    
%% load scope configuration information

    PS5000aConfig;
    
%% LCoS parameters
    
    wavelength = 520e-9;
    res_X = 1920;
    res_Y = 1080;
    pixel_pitch = 6.4e-6;
    frequency = 360;
    Vb = 1.4; Vw = 3.0;
    beam_contraction = 2;
    
%% stage parameters

    numOfAxis = 3;
    serialNo_X = '27600159';
    serialNo_Y = '27260933';
    serialNo_Z = [];
    offset_X = 9.75;
    offset_Y = 36.25;
    offset_Z = 0;
    
%% measuring parameters
    
    direction = 'perpendicular';     % 'parallel','perpendicular'
    waveform_record_length = 120000;
    center_X = 0; center_Y = 0;
    step_grid_X = -6:0.5:6; step_grid_Y = -3:0.5:3;    % full frame 0.5mm
    %step_grid_X = -6:1:6; step_grid_Y = -3:1:3;    % full frame 1.0mm
    %step_grid_X = -6:1.5:6; step_grid_Y = -3:1.5:3;    % full frame 1.5mm
    %step_grid_X = -6:2:6; step_grid_Y = -3:2:3;    % full frame 2.0mm
    %step_grid_X = [-3 0 3]; step_grid_Y = [-1.5 0 1.5];   % quater frame test
    step_cnt_X = length(step_grid_X); step_cnt_Y = length(step_grid_Y);
    
    use_calibration = false;        % calibration file
    result_from_calibration = use_calibration;
    if result_from_calibration
        load('cal_profile_1.3-2.7_2x.mat');
    end
    
%% connecting

    if exist('stage', 'var')
        try
            ShutDown(stage);
        catch
        end
    end
    stage = Thorlabs_Translation_Stage( numOfAxis, serialNo_X,serialNo_Y,serialNo_Z );
    Connect(stage, 10);
    
%% scope device connection

    % Check if an Instrument session using the device object |ps5000aDeviceObj|
    % is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
    if (exist('ps5000aDeviceObj', 'var') && ps5000aDeviceObj.isvalid && strcmp(ps5000aDeviceObj.status, 'open'))
        openDevice = questionDialog(['Device object ps5000aDeviceObj has an open connection. ' ...
            'Do you wish to close the connection and continue?'], ...
            'Device Object Connection Open');
        if (openDevice == PicoConstants.TRUE)
            % Close connection to device.
            disconnect(ps5000aDeviceObj);
            delete(ps5000aDeviceObj);
        else
            % Exit script if User selects 'No'.
            ShutDown(stage);
            return;
        end
    end

    % Create a device object. 
    ps5000aDeviceObj = icdevice('picotech_ps5000a_generic', ''); 
    % Connect device object to hardware.
    connect(ps5000aDeviceObj);

%% set channels
    % Channels       : 0 - 1 (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A & PS5000A_CHANNEL_B)
    % Enabled        : 1 (PicoConstants.TRUE)
    % Type           : 1 (ps5000aEnuminfo.enPS5000ACoupling.PS5000A_DC)
    % Range          : 8 (ps5000aEnuminfo.enPS5000ARange.PS5000A_5V)
    % Analog Offset  : 0.0 V

    % Find current power source
    [status.currentPowerSource] = invoke(ps5000aDeviceObj, 'ps5000aCurrentPowerSource');
    if (ps5000aDeviceObj.channelCount == PicoConstants.QUAD_SCOPE && status.currentPowerSource == PicoStatus.PICO_POWER_SUPPLY_CONNECTED)
        [status.setChC] = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', 2, 0, 1, 8, 0.0);
        [status.setChD] = invoke(ps5000aDeviceObj, 'ps5000aSetChannel', 3, 0, 1, 8, 0.0);
    end

%% set device resolution

    % Max. resolution with 2 channels enabled is 15 bits.
    [status.setResolution, data_resolution] = invoke(ps5000aDeviceObj, 'ps5000aSetDeviceResolution', 15);

%% Verify timebase index and maximum number of samples

    % for 15bit resolution:
    % sample_interval = (timebase–2) / 125,000,000Hz = (timebase–2) * 8ns
    status.getTimebase2 = PicoStatus.PICO_INVALID_TIMEBASE;
    timebaseIndex = 127;    % this equals to 1us interval

    while (status.getTimebase2 == PicoStatus.PICO_INVALID_TIMEBASE)
        [status.getTimebase2, timeIntervalNanoseconds, maxSamples] = invoke(ps5000aDeviceObj, 'ps5000aGetTimebase2', timebaseIndex, 0);
        if (status.getTimebase2 == PicoStatus.PICO_OK)
            break;
        else
            timebaseIndex = timebaseIndex + 1;
        end
    end
    timeInterval = timeIntervalNanoseconds / 1e9;
    fprintf('Timebase index: %d, sampling interval: %d ns\n', timebaseIndex, timeIntervalNanoseconds);
    clearvars timeIntervalNanoseconds;

    % Configure the device object's |timebase| property value.
    set(ps5000aDeviceObj, 'timebase', timebaseIndex);

%% set trigger
    % Set a trigger on channel A, with an auto timeout - the default value for
    % delay is used.

    % Trigger properties and functions are located in the Instrument
    % Driver's Trigger group.
    triggerGroupObj = get(ps5000aDeviceObj, 'Trigger');
    triggerGroupObj = triggerGroupObj(1);

    % Set the |autoTriggerMs| property in order to automatically trigger the
    % oscilloscope after 0.01 second if a trigger event has not occurred. Set to 0
    % to wait indefinitely for a trigger event.
    set(triggerGroupObj, 'autoTriggerMs', 10);

    % Channel     : 0 (ps5000aEnuminfo.enPS5000AChannel.PS5000A_CHANNEL_A)
    % Threshold   : 10000 mV
    % Direction   : 2 (ps5000aEnuminfo.enPS5000AThresholdDirection.PS5000A_RISING)
    [status.setSimpleTrigger] = invoke(triggerGroupObj, 'setSimpleTrigger', 0, 10000, 2);

%% set block parameters
    % Capture a block of data and retrieve data values for channels A and B.

    % Block data acquisition properties and functions are located in the 
    % Instrument Driver's Block group.
    blockGroupObj = get(ps5000aDeviceObj, 'Block');
    blockGroupObj = blockGroupObj(1);

    % Set pre-trigger and post-trigger samples as required - the total of this
    % should not exceed the value of |maxSamples| returned from the call to
    % |ps5000aGetTimebase2()|. The number of pre-trigger samples is set in this
    % example but default of 10000 post-trigger samples is used.

    % Set pre-trigger samples.
    set(ps5000aDeviceObj, 'numPreTriggerSamples', 0);
    set(ps5000aDeviceObj, 'numPostTriggerSamples', waveform_record_length);

%% settings for retrieving data

% --capture function--:
% [status.runBlock] = invoke(blockGroupObj, 'runBlock', 0);

% Retrieve data values:
startIndex              = 0;
segmentIndex            = 0;
downsamplingRatio       = 1;
downsamplingRatioMode   = ps5000aEnuminfo.enPS5000ARatioMode.PS5000A_RATIO_MODE_NONE;

% --retrieve function--:
% [numSamples, overflow, chA, chB] = invoke(blockGroupObj, 'getBlockData', startIndex, segmentIndex, downsamplingRatio, downsamplingRatioMode);

%% stage homing

    if IsHomed(stage) ~= 1
        Home(stage);
    else
        fprintf('Stage already homed.\n\n');
    end
    SetSoftwareHome(stage, [ offset_X, offset_Y, offset_Z ]);
    Return(stage);
    
%% display area prepare
    
    close all;
    f_ctrl = figure('Name','Direct Display','NumberTitle','on','WindowState','maximized','MenuBar','none','ToolBar','none','Position',[3840 400 1920 1080]);
    img_ctrl = imshow( 0*ones(res_Y,res_X), [0 255] );
    set(gca,'DataAspectRatioMode','auto');
    set(gca,'Position',[0 0 1 1]);
    f_ctrl.Position = [3840 360 1920 1080];
    pause(1);   % wait for the display window to stabilize
    
%% auto stepping and measuring
    
    result_waveform_period = 1 / ( timeInterval * frequency );    % 360Hz LC voltage drive
    point_result.result_avg_voltage = zeros(1,256);
    point_result.result_avg_temp = zeros(1,256);
    point_result.result_voltage_waveform_timeavg = zeros( round(result_waveform_period), 256 );
    point_result.result_voltage_waveform_timeavg_kalman = zeros( round(result_waveform_period), 256 );
    
    % guarantee positive integer number
    step_cnt_X = max(1,floor(step_cnt_X));
    step_cnt_Y = max(1,floor(step_cnt_Y));
    result_group_original = cell(step_cnt_Y,step_cnt_X);
    result_group_retarded = cell(step_cnt_Y,step_cnt_X);
    % stepping for 2 times
    for measure_index = 1:1
        % hint for QWP
        if measure_index == 1
            questdlg('Please filp DOWN the QWP', 'QWP State', 'Yes','是','Yes');
        elseif measure_index == 2
            questdlg('Please filp UP the QWP', 'QWP State', 'Yes','是','Yes');
        end
        measuring_time_stamp = tic;
        % measuring
        for step_index_Y = 1 : step_cnt_Y
            pos_Y = center_Y + step_grid_Y(step_index_Y);
            for step_index_X = 1 : step_cnt_X
                pos_X = center_X + step_grid_X(step_index_X);
                % stage movement
                point_result.pos = [ pos_X, pos_Y, 0 ];
                Move(stage, [ offset_X+pos_X, offset_Y-pos_Y, offset_Z ], 2);   % Y axis is invert-mounted
                % data measuring
                for gray_level = 0:255
                    % display
                    if result_from_calibration
                	    display_bitmap = cal_img(:,:,gray_level+1);            % calibrated
                    else
                        display_bitmap = gray_level * ones( res_Y, res_X );    % uncalibrated
                    end
                    img_ctrl.CData = display_bitmap;
                    img_ctrl.CDataMapping = 'direct';
                    % display stabilize
                    if gray_level == 0 
                        pause(1); 
                    else
                        pause(0.2); 
                    end
                    % capture
                    [status.runBlock] = invoke(blockGroupObj, 'runBlock', 0);
                    % retrieve
                    [~, ~, chA, chB] = invoke(blockGroupObj, 'getBlockData', startIndex, segmentIndex, downsamplingRatio, downsamplingRatioMode);
                    point_result.result_avg_temp( gray_level+1 ) = mean( (chB/50) );
                    % data post-process
                    temp_waveform_timeavg = zeros( round(result_waveform_period), 1 );
                    for index = 0 : floor( waveform_record_length/result_waveform_period ) - 1
                        temp_startpoint = round( index * result_waveform_period ) + 1;
                        temp_waveform_timeavg = temp_waveform_timeavg + chA( temp_startpoint : temp_startpoint + round(result_waveform_period) - 1 );
                    end
                    temp_waveform_timeavg = temp_waveform_timeavg / (index+1);
                    temp_waveform_timeavg_kalman = kalmanFilter( [temp_waveform_timeavg;temp_waveform_timeavg] );
                    temp_waveform_timeavg_kalman = temp_waveform_timeavg_kalman( round(result_waveform_period)+1:2*round(result_waveform_period) );
                    % data save
                    point_result.result_voltage_waveform_timeavg( : , gray_level+1 ) = temp_waveform_timeavg;
                    point_result.result_voltage_waveform_timeavg_kalman( : , gray_level+1 ) = temp_waveform_timeavg_kalman;
                    point_result.result_avg_voltage( gray_level+1 ) = mean( temp_waveform_timeavg_kalman );
                    fprintf( ' -- @[%4.2f %4.2f %4.2f] Sample:%d/%d, lvl %3d done, result is %6.2fmv @%4.2fC\n', ...
                        point_result.pos, (step_index_Y-1)*step_cnt_X+step_index_X, step_cnt_X*step_cnt_Y, gray_level, ...
                        point_result.result_avg_voltage( gray_level+1 ), point_result.result_avg_temp( gray_level+1 ) );
                    fprintf( ' -- Measurement No.%d, time elapsed: %ds ...\n\n', measure_index, fix(toc(measuring_time_stamp)) );
                end
                clearvars gray_level index chA chB temp_*;
                % data save
                if measure_index == 1
                    result_group_original{ step_index_Y, step_index_X } = point_result;
                elseif measure_index == 2
                    result_group_retarded{ step_index_Y, step_index_X } = point_result;
                end
            end
        end
    end
    clearvars step_index_* pos_* gray_level;
    
    Return(stage);
    close(f_ctrl);
    
%% stop the device

    % Stop ps5000a
    [status.stop] = invoke(ps5000aDeviceObj, 'ps5000aStop');
    
%% disconnect the device

    % Disconnect device object from hardware.
    ShutDown(stage);
    disconnect(ps5000aDeviceObj);
    delete(ps5000aDeviceObj);
    
%% result processing
    
    % timetable for easier waveform analysis
    waveform_timeline = seconds( timeInterval : timeInterval : round(result_waveform_period) * timeInterval ).';
    waveform_timetable = timetable( waveform_timeline, point_result.result_voltage_waveform_timeavg_kalman );
    
    % intensity_parallel = I0/2 * ( 1 + cos(phase_ret) )
    % intensity_perpendicular = I0/2 * ( 1 - cos(phase_ret) )

%% data save
    
    % format
    time = now;
    time_str = datestr( time, 'yyyymmdd-HHMMSS' );
    if result_from_calibration
        result_file_name = sprintf( 'result_%dnm_%.1f-%.1f_%dx_postcal_%s.mat', wavelength*1e9,Vb,Vw,beam_contraction,time_str );
    else
        result_file_name = sprintf( 'result_%dnm_%.1f-%.1f_%dx_%s.mat', wavelength*1e9,Vb,Vw,beam_contraction,time_str );
    end
    saved_vars = ["data_resolution","direction","frequency","pixel_pitch","res_*","result_from_calibration",...
                "result_group_*","step_grid_*","timeInterval","Vb","Vw","waveform_record_length","wavelength"];
    % and save
    save( result_file_name, "time", '-v7.3' );   % time stamp, also create the file
    for index = 1 : length(saved_vars)
        try
            save( result_file_name, saved_vars(index), '-append' );
        catch
            fprintf( '  -!  Variable %s not found.  !-  \n', saved_vars(index) );
        end
    end
    fprintf( '\n\n  --  Results saved.  --  \n' );
    clearvars time time_str;

%% result display
    
    close all;
    disp_point_index = sub2ind(size(result_group_original),ceil(size(result_group_original,1)/2),ceil(size(result_group_original,2)/2));    % center
    point_result = result_group_original{disp_point_index};
    
    f_result = figure('Name','Result','NumberTitle','off','WindowState','Maximized');
    % raw
    a_result_1 = subplot(2,1,1); 
    plot( 0:255, point_result.result_avg_voltage, 'b', ...
          0:255, max(point_result.result_voltage_waveform_timeavg), 'c', ...
          0:255, min(point_result.result_voltage_waveform_timeavg), 'c' );
    axis( [ 0, 255,  0, max(max(point_result.result_voltage_waveform_timeavg)) ] );
    xlabel(a_result_1, 'Gray level');
    ylabel(a_result_1, 'Voltage / Intensity (mV)');
    grid(a_result_1, 'on');
    title(a_result_1, [ num2str(1e9*wavelength) 'nm - intensity curve (with flicker) timeavg' ]);
    % time averaged
    a_result_2 = subplot(2,1,2); 
    plot( 0:255, point_result.result_avg_voltage, 'b', ...
          0:255, max(point_result.result_voltage_waveform_timeavg_kalman), 'c', ...
          0:255, min(point_result.result_voltage_waveform_timeavg_kalman), 'c' );
    axis( [ 0, 255,  0, max(max(point_result.result_voltage_waveform_timeavg_kalman)) ] );
    xlabel(a_result_2, 'Gray level');
    ylabel(a_result_2, 'Voltage / Intensity (mV)');
    grid(a_result_2, 'on');
    title(a_result_2, [ num2str(1e9*wavelength) 'nm - intensity curve (with flicker) timeavg+kalman' ]);
    
    % temp
    f_comparison = figure('Name','Comparison (timeavg+kalman)','NumberTitle','off','WindowState','Maximized');
    data_point = 70:10:150;
    for i = 1:9
        a_comparison(i) = subplot(3,3,i);
        plot(point_result.result_voltage_waveform_timeavg_kalman( 1:round(result_waveform_period) , data_point(i)+1 ));
        title(a_comparison(i), data_point(i));
    end
    clearvars i data_point a_comparison temp;
       
%% end
    
    fprintf( '\n' );
    toc(start_time_stamp);
    
    