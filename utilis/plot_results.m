function plot_results(azimuth_stds, rmse_ekf, rmse_lg)
    % plot_results: Generates a comparison plot for the Monte Carlo results
    
    figure('Color', 'w', 'Position', [150, 150, 700, 500]);
    
    % Plot EKF results
    plot(azimuth_stds, rmse_ekf, 'r-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'r'); 
    hold on;
    
    % Plot Lie Group EKF results
    plot(azimuth_stds, rmse_lg, 'b-s', ...
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'b');
    
    % Formatting
    grid on;
    set(gca, 'FontSize', 11);
    
    % Labels and Title
    xlabel('Azimuth Standard Deviation (^\circ)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('RMSE of Position (m)', 'FontSize', 13, 'FontWeight', 'bold');
    title('Monte Carlo Comparison: D-EKF vs. D-LG-EKF', 'FontSize', 15);
    
    % Legend
    legend({'Standard Discrete EKF', 'Discrete Lie Group EKF'}, ...
        'FontSize', 12, ...
        'Location', 'northwest', ...
        'Box', 'on');
    
    % Adjust limits slightly for better visibility
    xlim([min(azimuth_stds)-0.5, max(azimuth_stds)+0.5]);
    ylim([0, max([rmse_ekf; rmse_lg]) * 1.2]);
    
    fprintf('\nPlotting complete. Check the figure window for results.\n');
end