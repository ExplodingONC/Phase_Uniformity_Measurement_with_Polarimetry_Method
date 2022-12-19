function [output_waveform] = kalmanFilter(input_waveform)
%KALMAN_FILTER 此处显示有关此函数的摘要
%   此处显示详细说明
    
    % data length
    length = size(input_waveform,1);
    % size of array. n_iter行，1列
    sz = [length, 1];
    % 过程方差，反应连续两个时刻温度方差。更改查看效果
    Q = 36e-4;
    % 测量方差，反应温度计的测量精度。更改查看效果
    R = 0.25;
    % 对数组进行初始化
    % 对温度的后验估计。即在k时刻，结合温度计当前测量值与k-1时刻先验估计，得到的最终估计值
    xhat = zeros(sz); 
    % 后验估计的方差
    P = zeros(sz); 
    % 温度的先验估计。即在k-1时刻，对k时刻温度做出的估计
    xhatminus = zeros(sz);
    % 先验估计的方差
    Pminus = zeros(sz);
    % 卡尔曼增益，反应了温度计测量结果与过程模型（即当前时刻与下一时刻温度相同这一模型）的可信程度
    K = zeros(sz); 
    % intial guesses
    xhat(1) = mean(input_waveform); %温度初始估计值为23.5度
    P(1) =1; % 误差方差为1
   
    for k = 2:length
        % 时间更新（预测）
        % 用上一时刻的最优估计值来作为对当前时刻的温度的预测
        xhatminus(k) = xhat(k-1);
        % 预测的方差为上一时刻温度最优估计值的方差与过程方差之和
        Pminus(k) = P(k-1)+Q;
        % 测量更新（校正）
        % 计算卡尔曼增益
        K(k) = Pminus(k)/( Pminus(k)+R );
        % 结合当前时刻温度计的测量值，对上一时刻的预测进行校正，得到校正后的最优估计。该估计具有最小均方差
        xhat(k) = xhatminus(k)+K(k)*(input_waveform(k)-xhatminus(k));
        % 计算最终估计值的方差
        P(k) = (1-K(k))*Pminus(k);
    end

    output_waveform = xhat;

end

