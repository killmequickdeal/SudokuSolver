import image_functions
img_funcs = image_functions();

close all;
global savedir;
global testname;
savedir = 'C:\Users\Riley\Documents\MATLAB\Images';
testname = 'sudo-2';

sudo = imread(strcat(testname,'.jpg'));
sudo = rgb2gray(sudo);
% sudo = imresize(sudo, .1);
% sudo = imrotate(sudo, -90);

figure; imshow(sudo); title('Original Image');
filename = strcat(testname, '_original.png');
saveas(gcf, fullfile(savedir, filename));

original_image = sudo;
[rows,cols,channel] = size(sudo);


sudo_blur = img_funcs.gaussian_avg(sudo, .5, 5);
sudo = sudo_blur/max(sudo_blur(:));
figure; imshow(sudo); title('Gaussian blur; sigma = .5; radius=5');
filename = strcat(testname, '_blurred.png');
saveas(gcf, fullfile(savedir, filename));

sudo = img_funcs.adaptive_thresh(sudo, .9, 11);
figure; imshow(sudo); title('Adaptive Thresh; Sensitivity = .9; Neighborhood size = 11');
filename = strcat(testname, '_adaptive_thresh.png');
saveas(gcf, fullfile(savedir, filename));

sudo = img_funcs.dilate_image(sudo);
figure; imshow(sudo); title('Dilated Image');
filename = strcat(testname, '_dilated.png');
saveas(gcf, fullfile(savedir, filename));

sudo = logical(img_funcs.find_largest_component(sudo));
figure; imshow(sudo); title('Largest Component');
filename = strcat(testname, '_largestcomponent.png');
saveas(gcf, fullfile(savedir, filename));


potential_corners = img_funcs.harris_corner_detector(sudo);

[sudo_X, sudo_Y, image_X, image_Y] = img_funcs.detect_corners(potential_corners);


sudo = img_funcs.projective_transform(original_image, sudo_X, sudo_Y, image_X, image_Y);
figure; imshow(sudo); title('Projective transform');
filename = strcat(testname, '_transform.png');
saveas(gcf, fullfile(savedir, filename));

transformed_sudo = sudo;

estimate_lines(sudo, rows, cols, img_funcs);
hough_lines(sudo, img_funcs);




% solve sudoku
global sudoku_matrix;
global unsolved_sudo;
% solution for test 1
% sudoku_matrix = [0,0,0,6,0,4,7,0,0;
%                  7,0,6,0,0,0,0,0,9;
%                  0,0,0,0,0,5,0,8,0;
%                  0,7,0,0,2,0,0,9,3;
%                  8,0,0,0,0,0,0,0,5;
%                  4,3,0,0,1,0,0,7,0;
%                  0,5,0,2,0,0,0,0,0;
%                  3,0,0,0,0,0,2,0,8;
%                  0,0,2,3,0,1,0,0,0;];
% solution for test 2
sudoku_matrix = [0,3,9,1,0,0,0,0,0;
                 4,0,8,0,6,0,0,0,2;
                 2,0,0,5,8,0,7,0,0;
                 8,0,0,0,0,0,0,0,0;
                 0,2,0,0,0,9,0,0,0;
                 3,0,6,0,0,0,0,4,9;
                 0,0,0,0,1,0,0,3,0;
                 0,4,0,3,0,0,0,0,8;
                 7,0,0,0,0,0,4,0,0;];

unsolved_sudo = sudoku_matrix;

solved_sudo = solve();

% place numbers
estimate_num_locs(transformed_sudo, rows, cols);

% solving code adapted from https://www.youtube.com/watch?v=G_UYXzGuqvM
function result = possible(x,y,n)
    global sudoku_matrix
    % if any value in x range is equal to the test number, it doesnt work
    for i = 1:9
        if sudoku_matrix(x,i) == n
            result = false;
            return;
        end
    end
    % if any value in y range is equal to test number, it doesnt work
    for i = 1:9
        if sudoku_matrix(i,y) == n
            result = false;
            return;
        end
    end
    
    x0 = floor((x-1)/3)*3;
    y0 = floor((y-1)/3)*3;
    % if any value in shared 3x3 block is equal to test number, it doesnt
    % work
    for i = 1:3
        for j = 1:3
            if sudoku_matrix(x0+j, y0+i) == n
                result = false;
                return;
            end
        end
    end
    
    % if we make it past these conditions, it can work!
    result = true;
end

function result = solve()
    global sudoku_matrix
    result = sudoku_matrix;
    % for each cell in sudoku
    for x = 1:9
        for y = 1:9
            % if the cell is empty
            if sudoku_matrix(x,y) == 0
                % for all possible values in a sudoku
                for i = 1:9
                    % if the value can fit in the location
                    if possible(x,y,i)
                        % try it, and solve for the next cell recursively
                        % until nothing works, this means somewhere a cell
                        % was incorrect so go back up the recursive tree
                        sudoku_matrix(x,y) = i;
                        result = solve();
                        % if every value in the sudoku has nonzero values,
                        % it is solved so we can return
                        if all(sudoku_matrix, 'all')
                            return
                        end
                        sudoku_matrix(x,y) = 0;
                    end
                end
                return;
            end
        end
    end
    result = sudoku_matrix;
end

function estimate_num_locs(sudo, rows, cols)
    global sudoku_matrix;
    global unsolved_sudo;
    global testname;
    global savedir;
    row_gap = round(rows/9);
    col_gap = round(cols/9);
    figure; imshow(sudo); title("Solved Sudoku");
    % for each cell in matrix, check if the matrix was missing a number
    % and if it was - try to place the solution in the middle of the cell
    for j = 1:9
        for i = 1:9
            if unsolved_sudo(j,i) == 0
                text((col_gap*i) - col_gap/2,(row_gap*j) - row_gap/2, string(sudoku_matrix(j,i)), 'FontSize', 20, 'Color', 'cyan');
            end
        end
    end
    filename = strcat(testname, '_solved.png');
    saveas(gcf, fullfile(savedir, filename));
end

function estimate_lines(sudo, rows, cols, img_funcs)
    row_gap = round(rows/9);
    col_gap = round(cols/9);
    % make 9 rows and 9 cols estimated from image size and sudoku structure
    figure; imshow(sudo); title('Estimated Lines');
    for i = 1:9
        line([1,cols],[row_gap*i,row_gap*i], 'Color', 'red');
        line([col_gap*i, col_gap*i],[1, rows], 'Color', 'green');
    end
    global testname;
    global savedir;
    filename = strcat(testname, '_estimatedlines.png');
    saveas(gcf, fullfile(savedir, filename));

    % for each cell in sudoku, crop out the cell. To be used in image
    % recognition
    for j = 1:9
        for i = 1:9
            b1 = max(1,i*row_gap-row_gap);
            b2 = max(1,j*col_gap-col_gap);
            width = row_gap;
            height = col_gap;
            cropped_image = imcrop(sudo,[b1 b2 width height]);
%             figure; imshow(cropped_image);
%             cropped_image = img_funcs.adaptive_thresh(cropped_image, .9, 11);
%             figure; imshow(cropped_image);
%             res = logical(img_funcs.find_largest_component(cropped_image));
%             figure; imshow(res);
        end
    end
end

function hough_lines(sudo, img_funcs)
    transformed_sudo = sudo;
    sudo_blur = img_funcs.gaussian_avg(sudo, .5);
    sudo = sudo_blur/max(sudo_blur(:));

    [sudo, gradient_direction] = img_funcs.canny(sudo, 1, 1);
    figure; imshow(sudo); title('Canny Edge Detection');
    global testname;
    global savedir;
    filename = strcat(testname, '_canny.png');
    saveas(gcf, fullfile(savedir, filename));


    [hough_results, T, R] = img_funcs.hough_transform(sudo, gradient_direction, 1, -90:89);
    P = img_funcs.hough_peaks(hough_results, 200);
    figure; imshow(imadjust(rescale(hough_results)),'XData',T,'YData',R,...
        'InitialMagnification','fit'); 
    xlabel('\theta'), ylabel('\rho');
    axis on;
    axis normal;
    hold on;
    plot(P(:,2),P(:,1),'r*');
    title('Hough peaks');
    colormap(gcf,hot);
    filename = strcat(testname, '_hough.png');
    saveas(gcf, fullfile(savedir, filename));
    
    img_funcs.hough_lines(transformed_sudo, hough_results, P);
end