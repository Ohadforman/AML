function CSV_to_plot_all_in_dir(dir_path)
% Plot all CSV files in the given directory

% Get a list of all the CSV files in the directory
csv_files = dir(fullfile(dir_path, '*.csv'));

% Loop over each CSV file and generate a plot
for i = 1:numel(csv_files)
    % Load CSV file and extract X and Y data columns
    csv_path = fullfile(csv_files(i).folder, csv_files(i).name);
    [~, file_name, ~] = fileparts(csv_path);
    data = readtable(csv_path);
    x_col_name = data.Properties.VariableNames{1};
    y_col_name = data.Properties.VariableNames{2};
    x_data = table2array(data(:, 1));
    y_data = table2array(data(:, 2));

    % Plot the data
    figure;
    plot(x_data, y_data, 'Color', [i/numel(csv_files), (i/numel(csv_files))^2, (i/numel(csv_files))^3]);
    xlabel(x_col_name);
    ylabel(y_col_name);
    title(file_name);
end
