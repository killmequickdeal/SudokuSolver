classdef image_functions

    properties
        
    end
    
    methods
        
        % used for equalizing contrast by redistribution of gray levels
        function [new_image, hist] = histogram_equalization(~, image)
            % create a histogram of elements values from the image to
            % determine frequencies
            hist = imhist(image);
            [x,y,z] = size(image);
            N = numel(hist);
            c = zeros(256,1);
            
            % cumulative distribution of area under histogram
            c(1,1) = hist(1,1)/N;
            for i = 2:256
                c(i,1) = c(i-1,1) + hist(i,1)/N;
            end
            
            % normalize c
            for i = 1:256
               c(i,1) = c(i,1)*(255/(c(256,1) - c(1,1))); 
            end
            
            % find image(i,j) value in c and set the new value in image to
            % the value in c
            for i = 1:x
                for j = 1:y
                    image(i,j) = c(image(i,j)+1, 1);
                end
            end
            
            hist = c;
            new_image = image;
        end

        
        % goal is to compress the range of image thus bringing out details
        % in dark regions
        function new_image = log_mapping(~, image)
            [x,y,z] = size(image);
            M = max(image(:));
            
            % calculate the scaling constant from the maximum image value
            c = 255/log10(1.0*double(M));
            % for each value in image take the log and scale it
            for i = 1:x
                for j = 1:y
                    image(i,j) = c*log10(1 + double(image(i,j)));
                end
            end
            
            new_image = image;
        end

        % rotate an image
        function new_image = rotate(~, image, theta)
            [x,y,z] = size(image);
            % determine the new rows and columns after rotation. +2 for
            % added padding, as the rounding isnt perfect.
            rows = ceil(x*abs(cosd(theta)) + y*abs(sind(theta))) + 2;
            cols = ceil(x*abs(sind(theta)) + y*abs(cosd(theta))) + 2;
            
            % calculate midpoints of old and new sizes
            mid_original_x = round((x)/2);
            mid_original_y = round((y)/2);
            
            mid_final_x = ceil(rows/2);
            mid_final_y = ceil(cols/2);
            
            new_image = zeros(rows, cols);
         
            for i = 1:rows
                for j = 1:cols
                    % determine the place in original image which the i,j
                    % in the new image needs to come from
                 
                    inputRow = round((i-mid_final_x)*cosd(theta) + (j-mid_final_y)*sind(theta)) + mid_original_x;
                    inputCol = round(-(i-mid_final_x)*sind(theta) + (j-mid_final_y)*cosd(theta)) + mid_original_y;
                    
                    % if input location is outside the image, leave it 0
                    if inputRow <= x && inputCol <= y && inputRow > 0 && inputCol > 0
                        new_image(i, j) = image(inputRow, inputCol);
                    end
                end
            end
            
            new_image = uint8(new_image);
            
        end

        % goal is to blur images, and remove detail/noise
        function new_image = gaussian_avg(~, image, sigma, varargin)
            [x,y,z] = size(image);
            % vary radius based on sigma
            if length(varargin) == 0
                if sigma < 1
                    radius = 2;
                else
                    radius = round(2*sigma);
                end
            else
                radius = varargin{1};
            end
                
            
            kernel_size = 2*radius + 1;
            kernel = zeros(kernel_size,kernel_size);

            % build the gaussian kernel based on the gaussian 2D function
            for i = 1:kernel_size
                for j = 1:kernel_size
                    kernel(i,j) = (exp(-1*((i-(radius+1))^2 + (j-(radius+1))^2)/(2*sigma*sigma)))/(2*pi*sigma*sigma);
                end
            end
            % normalize filter so sum == 1
            normalization_constant = sum(sum(kernel));
            kernel = kernel / normalization_constant;
            
            
            new_image = image;
            % pad image with black, I am using Vignetting for dealing with
            % the borders
            image = padarray(image, [radius, radius]);
            % for each pixel of output image, apply the gaussian 5x5 kernel
            % to a 5x5 chunk of the original image. sum the results to get
            % the new value.
            for i = 1+radius:x+radius
                for j = 1+radius:y+radius
                    gaus_section = double(image(i-radius:i+radius, j-radius:j+radius));
                    new_image(i-radius,j-radius) = sum(sum(gaus_section.*kernel));
                end
            end
            new_image = double(new_image);
        end

        % removes isolated pixels/lines, preserves step edges. Doesnt blur
        % edges like smoothing filters(gaussian). Costly.
        function new_image = median_filtering(~, image)
            [x,y,z] = size(image);
            new_image = image;
            % for each value in image (ignoring outside border)
            for i = 2:x-1
               for j = 2:y-1
                  % take a 3x3 chunk, and find the median value
                  % this value is the new image value
                  chunk = sort(reshape(image(i-1:i+1,j-1:j+1),[1,9]));
                  new_image(i,j) = chunk(1,5);
               end
            end
        end
        
        function response = harris_corner_detector(obj, image)
            image = gaussian_avg(obj, image, .5);
%             figure; imshow(image/max(image(:)));
%             title('Gaussian masked building.pnm');
            % Spatial derivative calculation
            [x_gradient, y_gradient] = sobel(obj, image);
            
%             figure; imshow(x_gradient/max(x_gradient(:)));
%             title('X gradient building.pnm');
%             figure; imshow(y_gradient/max(y_gradient(:)));
%             title('Y gradient building.pnm');

            % structure tensor setup
            x_gradient_sq = x_gradient.*x_gradient;
            y_gradient_sq = y_gradient.*y_gradient;
            xy_gradient_sq = x_gradient.*y_gradient;
%             figure; imshow(x_gradient_sq/max(x_gradient_sq(:)));
%             title('X gradient squared building.pnm');
%             figure; imshow(y_gradient_sq/max(y_gradient_sq(:)));
%             title('Y gradient squared building.pnm');
%             figure; imshow(xy_gradient_sq/max(xy_gradient_sq(:)));
%             title('XY gradient building.pnm');
            
          
            x_gradient_sq = gaussian_avg(obj, x_gradient_sq, 1);
            y_gradient_sq = gaussian_avg(obj, y_gradient_sq, 1);
            xy_gradient_sq = gaussian_avg(obj, xy_gradient_sq, 1);
%             figure; imshow(x_gradient_sq/max(x_gradient_sq(:)));
%             title('X gradient squared, Smoothed sigma=1 building.pnm');
%             figure; imshow(y_gradient_sq/max(y_gradient_sq(:)));
%             title('Y gradient, Smoothed sigma=1 building.pnm');
%             figure; imshow(xy_gradient_sq/max(xy_gradient_sq(:)));
%             title('XY gradient, Smoothed sigma=1 building.pnm');
            % Harris response calculation
            alpha = .04;
            R = x_gradient_sq.*y_gradient_sq - xy_gradient_sq.^2 - alpha .* (x_gradient_sq + y_gradient_sq).^2;
            maxVal = max(R, [], 'all');
%             figure; imshow(R/max(R(:)));
%             title('Harris Response building.pnm');

            % Non maximum suppression
            response = window_suppression(obj, R, maxVal);
%             figure; imshow(response);
%             title('Non-maximum suppression building.pnm');

        end
        
        % Perform non-maximum suppression with a 3x3 window keeping only
        % the highest value
        function results = window_suppression(~, R, maxVal)
            [rows,cols] = size(R);
            results = zeros(rows,cols);
            for i = 1:rows
                for j = 1:cols
                   if R(i,j) > .01*maxVal % threshold to only allow above a certain value compared to the max
                       MatSection = R(max(i-1, 1):min(i+1, rows), max(j-1, 1):min(j+1, cols));
                       MaxSectionVal = max(MatSection(:));
                       if MaxSectionVal == R(i,j)
                           results(i,j) = 1;
                       end
                   end
                end
            end
        end
        
        function [x_gradient, y_gradient] = sobel(~, image)
            sobel_ymask = [-1,-2,-1;0,0,0;1,2,1];
            sobel_xmask = transpose(sobel_ymask);

            % tried both conv2 and imfilter, but imfilter gave better
            % results for some reason. Matlab documentation says they
            % should be the same
            
%             x_gradient = conv2(image, sobel_xmask, 'same');
%             y_gradient = conv2(image, sobel_ymask, 'same');
            x_gradient = imfilter(image, sobel_xmask, 'conv', 'replicate');
            y_gradient = imfilter(image, sobel_ymask, 'conv', 'replicate');


        end
       
        function [results, gradient_direction] = canny(obj, image, sigma, suppress)
            [rows,cols,z] = size(image);
            
            % gaussian filter
            filtered_image = gaussian_avg(obj, image, sigma);
            
%             figure; imshow(filtered_image/max(filtered_image(:)));
%             title('Gaussian masked hinge.pnm');

            % find intensity gradients
            [x_gradient, y_gradient] = sobel(obj, filtered_image);
            
%             figure; imshow(x_gradient);
%             title('X gradient hinge.pnm');
%             figure; imshow(y_gradient);
%             title('Y gradient hinge.pnm');
            
            x_gradient = gaussian_avg(obj, x_gradient, sigma);
            y_gradient = gaussian_avg(obj, y_gradient, sigma);
            
%             figure; imshow(x_gradient);
%             title('X gradient, smoothed sigma=1 hinge.pnm');
%             figure; imshow(x_gradient);
%             title('X gradient, smoothed sigma=1 hinge.pnm');
            
            x_gradient_sq = x_gradient.*x_gradient;
            y_gradient_sq = y_gradient.*y_gradient;
            
%             figure; imshow(x_gradient_sq/max(x_gradient_sq(:)));
%             title('X gradient squared hinge.pnm');
%             figure; imshow(y_gradient_sq/max(y_gradient_sq(:)));
%             title('Y gradient squared hinge.pnm');
            
            gradient_magnitude = sqrt(x_gradient_sq + y_gradient_sq);
            gradient_direction = atan2d(y_gradient, x_gradient);
            gradient_direction = gradient_direction + 180;
%             gradient_direction = wrapTo360(gradient_direction);
            global testname;
            global savedir;
            
            figure; imshow(gradient_magnitude/max(gradient_magnitude(:)));
            title('Gradient magnitude');
            filename = strcat(testname, '_magnitude.png');
            saveas(gcf, fullfile(savedir, filename));
            figure; imshow(gradient_direction/max(gradient_direction(:)));
            title('Gradient direction');
            filename = strcat(testname, '_direction.png');
            saveas(gcf, fullfile(savedir, filename));
            
            if suppress == 1
                % apply non-maximum suppression
                suppressed_values = zeros(rows, cols);
                for i = 2:rows-1
                    for j = 2:cols-1
                        suppressed = true;
                        cur_direction = gradient_direction(i,j);
                        cur_magnitude = gradient_magnitude(i,j);
                        % x = 0 (0)
                        if (cur_direction >= 157.5 && cur_direction < 202.5) ||... 
                            (cur_direction >= 337.5 || cur_direction < 22.5)

                            if cur_magnitude >= gradient_magnitude(i, j-1) && cur_magnitude > gradient_magnitude(i, j+1)
                                % keep local maximum values
                                suppressed = false;
                            end
                        % y = 0 (90)
                        elseif (cur_direction >= 67.5 && cur_direction < 112.5) ||...
                                (cur_direction >= 247.5 && cur_direction < 292.5)

                            if cur_magnitude >= gradient_magnitude(i-1, j) && cur_magnitude > gradient_magnitude(i+1, j)
                                suppressed = false;
                            end
                        % y = x (45)        
                        elseif (cur_direction >= 22.5 && cur_direction < 67.5) ||...
                                (cur_direction >= 202.5 && cur_direction < 247.5)

                            if cur_magnitude >= gradient_magnitude(i-1, j-1) && cur_magnitude > gradient_magnitude(i+1, j+1)
                                suppressed = false;
                            end
                        % y = x (135)  
                        elseif (cur_direction >= 112.5 && cur_direction < 157.5) ||...
                                (cur_direction >= 292.5 && cur_direction < 337.5)

                            if cur_magnitude >= gradient_magnitude(i+1, j-1) && cur_magnitude > gradient_magnitude(i-1, j+1)
                                suppressed = false;
                            end
                        end

                        if suppressed == false
                            suppressed_values(i,j) = cur_magnitude;
                        end
                    end
                end

                % normalize the gradient
                gradient_magnitude = suppressed_values/max(suppressed_values(:));
                figure; imshow(gradient_magnitude);
                title('non-maximum suppression');
                filename = strcat(testname, '_nonmaximumsuppression.png');
                saveas(gcf, fullfile(savedir, filename));
                % calculate the thresholds for double thresholding based on the
                % "weighted mean" of the gradient. I'm not sure what this is
                % actually called but i just did the mean using only the
                % nonzero elements. before all the zero elements were pulling
                % my mean down making it kind of useless. This seems to work
                % well.
                weighted_mean = sum(sum(gradient_magnitude))/nnz(gradient_magnitude);
                % between 3:1 and 2:1 recommended by canny
                % most talk about 1.25 : .75 or 1.33 : .66, I just experimented
                % to find what worked for me 

                strong_thresh = weighted_mean*.8; % .8

                weak_thresh = weighted_mean*.26; %.26


                edge_strengths = zeros(rows,cols); % 0 is bad, 1 is weak, 2 is strong

                % build a stack of all strong edges for hysteresis
                import java.util.Stack;
                strong_stack = Stack();           

                for i = 1:rows
                    for j = 1:cols
                        current_magnitude = gradient_magnitude(i,j);
                        if current_magnitude <= weak_thresh
                            gradient_magnitude(i,j) = 0;
                        elseif current_magnitude > weak_thresh && current_magnitude < strong_thresh
                            edge_strengths(i,j) = 1;
                        else
                            edge_strengths(i,j) = 2;
                            strong_stack.push([i,j]);
                        end
                    end
                end

                figure; imshow(gradient_magnitude);
                title('Post Double Thresholding');
                filename = strcat(testname, '_thresholded.png');
                saveas(gcf, fullfile(savedir, filename));

                % hysteresis until out of strong edges
                while strong_stack.size() > 0
                    strong_edge = strong_stack.pop();
                    i = strong_edge(1);
                    j = strong_edge(2);
                    % for each strong edge check if any nearby weak edges
                    for k = -1:1:1
                        new_i = min(max(0, i+k), rows);
                        for m = -1:1:1
                            new_j = min(max(0, j+m), cols);

                            % make the nearby weak edges strong and add to
                            % stack
                            if edge_strengths(new_i, new_j) == 1
                                edge_strengths(new_i, new_j) = 2;
                                strong_stack.push([new_i, new_j]);
                            end
                        end
                    end
                end

                figure; imshow(gradient_magnitude);
                title('Post Hysteresis');
                filename = strcat(testname, '_hysteresis.png');
                saveas(gcf, fullfile(savedir, filename));

                % set the remaining weak edges to zero since they are not
                % connected to a strong edge
                for i = 1:rows
                    for j = 1:cols
                        if edge_strengths(i,j) == 1
                            gradient_magnitude(i,j) = 0;
                        end
                    end
                end

                figure; imshow(gradient_magnitude);
                title('Post Weak edge elimination');
                filename = strcat(testname, '_weakelimination.png');
                saveas(gcf, fullfile(savedir, filename));
                
                % max the intensity of any edges we kept
                for i = 1:rows
                    for j = 1:cols 
                        if gradient_magnitude(i,j) > 0
                            gradient_magnitude(i,j) = 255;
                        end
                    end
                end

                results = logical(gradient_magnitude);
            else
                gradient_magnitude = gradient_magnitude/max(gradient_magnitude(:));
                results = gradient_magnitude;
            end
        end
        
        function [results, T, R] = hough_transform(~, image, gradient_direction, RhoResolution, theta)
            [rows,cols,z] = size(image);
            gradient_direction = wrapTo180(gradient_direction);
            T = theta;
            
            % convert 0, 360 gradient direction into -90,90
            for i = 1:rows
                for j = 1:cols
                    if gradient_direction(i,j) > 90
                        gradient_direction(i,j) = round(gradient_direction(i,j) - 180);
                    elseif gradient_direction(i,j) < -90
                        gradient_direction(i,j) = round(gradient_direction(i,j) + 180);
                    else
                        gradient_direction(i,j) = round(gradient_direction(i,j));
                    end
                end
            end
            
            % https://www.mathworks.com/help/images/ref/hough.html for
            % formulas
            % Basically calculate the number a rho values, default rhoresolution 1 = #
            % pixels. 
            D = sqrt((rows-1)^2 + (cols-1)^2);
            nrho = 2*(ceil(D/RhoResolution)) + 1;
            diagonal = RhoResolution*ceil(D/RhoResolution);
            ntheta = length(theta);
            % cache sin and cos values for a decent speedup
            cos_computations = cosd(theta);
            sin_computations = sind(theta);
            R = zeros(1, nrho);
            R = -diagonal:diagonal;
            accumulator = zeros(nrho, ntheta);
            results = zeros(nrho, ntheta);
            
            % for each pixel in image
            for i=1:rows
                for j=1:cols
                    % if it is an edge
                    if image(i,j) == 1
                        % vote from +-22.5 degrees of the gradient
                        % direction. This should reduce the complexity of
                        % the hough space.
                        low_theta = round(max(gradient_direction(i,j) - 22.5, -89)) + 90;
                        high_theta = round(min(gradient_direction(i,j) + 22.5, 90)) + 90;
                        for cur_theta=low_theta:high_theta
                            r = round(j * cos_computations(cur_theta) + i * sin_computations(cur_theta)) + diagonal;
                            accumulator(r, cur_theta) = accumulator(r, cur_theta) + 1;
                        end
                    end
                end
            end
            
            results = accumulator;
        end
        
        function results = hough_peaks(~, H, k)
            [rows, cols, channels] = size(H);
            thresh = ceil(max(H(:))*.3); % between .1*max and .5*max seems common
            vals_above_thresh = zeros(rows, cols);
            
            % for all hough results
            for i=1:rows
               for j=1:cols
                   % check if the result is above a certain threshold
                   if H(i,j) > thresh
                       % check that the point is a local maximum (3x3 area)
                       MatSection = H(max(i-2, 1):min(i+2, rows), max(j-2, 1):min(j+2, cols));
                       MaxSectionVal = max(MatSection(:));
                       if MaxSectionVal == H(i,j)
                           vals_above_thresh(i,j) = H(i,j);
                           H(max(i-2, 1):min(i+2, rows), max(j-2, 1):min(j+2, cols)) = 0;
                       end
                   end
               end
            end
            % get the maximum K values from those which are valid
            [B, I] = maxk(vals_above_thresh(:), k);
            [row, col] = ind2sub([rows, cols], I);
            % normalize values
            results = [minus(row, rows/2), minus(col, cols/2)];
        end
        
        function hough_lines(~, image, H, P)
            [rows, cols] = size(image);
            [nrho, ntheta] = size(H);

            figure; imshow(image);
            axis on; title('Hough Lines'); hold on;
            
            % for each peak given
            for peaks = 1:size(P)
                rho = P(peaks,1);% - nrho/2;
                theta = P(peaks,2);% - ntheta/2;
                
                % calculate two endpoints based on closeness to horiz/vert
                % line
                if theta <= 45 && theta >= -45
                    x1 = 0;
                    x2 = cols;
                    y1 = (-(cosd(theta)/sind(theta)))*x1 + rho/sind(theta);
                    y2 = (-(cosd(theta)/sind(theta)))*x2 + rho/sind(theta);
                else
                    y1 = 0;
                    y2 = rows;
                    x1 = (-(sind(theta)/cosd(theta)))*y1 + rho/cosd(theta);
                    x2 = (-(sind(theta)/cosd(theta)))*y2 + rho/cosd(theta);
                end
                
                slope = (y2-y1)/(x2-x1);
                y_int = y1-(slope*x1);
                x_int = -1*y_int / slope;
                line([x1,x2],[y1,y2], 'Color', 'green');
            end
            global testname;
            global savedir;
            filename = strcat(testname, '_houghlines.png');
            saveas(gcf, fullfile(savedir, filename));
        end
        
        function RANSAC(~, image)
            import ransactest;

            % treat edge points as points in an xy plane
            figure; imshow(image); 
            [row,col] = find(image);
            data = [col,row];
           
            % for i = 1: Trials
            %   select 2 data points randomly
            %   estimate features
            %   if sufficient number of nearby points
            %       success!
            %   end
            % end
            %
            % fail
            ransactest(data,2,100,10,0.1);
%             plot(data(:,1),data(:,2),'o')
%             set(gca,'YDir','reverse');
        end
        
        function results = adaptive_thresh(~, image, sensitivity, neighborhood_size)
            [rows, cols, channels] = size(image);
            
            % determine size of neighborhoods using the matlab defaults
            % from here https://www.mathworks.com/help/images/ref/adaptthresh.html
            % this method will use a set number of neighborhood divisions
            % for simplicity. For more accurate thresholding perform it
            % around every pixel.
            neighborhood = [neighborhood_size,neighborhood_size];
            x_neigh = ceil(rows / neighborhood(1,1));
            y_neigh = ceil(cols / neighborhood(1,2));
            
            % for each neighborhood
            for i = 1:x_neigh
                for j = 1:y_neigh
                    
                    % create the start and end pixel locations
                    x_start = max(1, neighborhood(1,1)*(i-1)+1);
                    x_end = min(rows, neighborhood(1,1)*i);
                    y_start = max(1, neighborhood(1,2)*(j-1)+1);
                    y_end = min(cols, neighborhood(1,2)*j);
                    
                    % calculate the local mean
                    MatSection = image(x_start:x_end, y_start:y_end);
                    matsize = size(MatSection);
                    N = matsize(1,1) * matsize(1,2);
                    total = sum(sum(MatSection));
                    average = total/N;
                    
                    % for all pixels in neighborhood if they are below 
                    % defined percent of the mean, set to black,
                    % otherwise white
                    for k = x_start:x_end
                        for m = y_start:y_end
                            if image(k,m) > average*sensitivity
                                image(k,m) = 0;
                            else
                                image(k,m) = 255;
                            end
                        end
                    end
                end
            end
            results = logical(image);
        end
        
        function results = dilate_image(~, image)
            [rows, cols, channels] = size(image);

            dilated_image=zeros(rows,cols);
            % perform a 3x3 'diamond' dilation on the image
            % 0 1 0
            % 1 1 1 
            % 0 1 0
            for i = 1:rows
                for j = 1:cols
                    % if value is 1, stretch in 4 major directions
                    if image(i,j) == 1
                        dilated_image(i,j)=1;
                        dilated_image(max(i-1, 1), j) = 1;
                        dilated_image(min(i+1, rows), j) = 1;
                        dilated_image(i, max(j-1,1)) = 1;
                        dilated_image(i, min(j+1, cols)) = 1;
                    end
                end
            end
            results = logical(dilated_image);
        end
        
        function results = negate_image(~, image)
            [rows, cols, channels] = size(image);
            results = zeros(rows, cols);
            for i = 1:rows
                for j = 1:cols
                   if image(i,j) == 1
                      results(i,j) = -1; 
                   end
                end
            end
        end
        
        function LB = recursive_conn_comp(obj, B)
            [rows, cols, channels] = size(B);
            % set foreground to -1 for labeling prep
            LB = negate_image(obj, B);
            label = 0;
            LB = find_components(obj, LB, label, rows, cols);
        end
        
        function LB = find_components(obj, LB, label, rows, cols)
            % for each index check if it is part of the foreground
            for i = 1:rows
                for j = 1:cols
                    % if part of foreground then label it and check the 8x8
                    % neighborhood
                    if LB(i,j) == -1
                       label = label + 1; 
                       LB = search(obj, LB, label, i, j, rows, cols);
                    end
                end
            end
        end
        
        function LB = search(obj, LB, label, i, j, rows, cols)
           % set the index to be part of the current component
           LB(i,j) = label;
           % check if any neighbors are also a member, if so then recurse
           for k=max(1,i-1):min(rows, i+1)
               for L = max(1, j-1):min(cols, j+1)
                   if LB(k,L) == -1
                      LB = search(obj, LB, label, k, L, rows, cols); 
                   end
               end
           end
        end
        
        function result = find_largest_component(obj, image)
            [rows, cols, channels] = size(image);
            
            % find all connected components
            result = recursive_conn_comp(obj, image);
            coloredLabels = label2rgb (result, 'hsv', 'k', 'shuffle');
            figure; imshow(coloredLabels); title("Colored connected components");
             global testname;
            global savedir;
            filename = strcat(testname, '_coloredcomponents.png');
            saveas(gcf, fullfile(savedir, filename));
            % find the number of connected components
            num_components = max(result(:));
            max_comp = -1;
            max_comp_area = -1;
            % for each connected component
            for i = 1:num_components
                % find the number of occurences
                occurence_matrix = result == i;
                num_cur_comp = sum(occurence_matrix(:));
                % if this is the new largest blob, set it
                if num_cur_comp > max_comp_area
                    max_comp = i;
                    max_comp_area = num_cur_comp;
                end
            end
            
            % set all non largest blobs to 0
            for i = 1:rows
                for j = 1:cols
                    if result(i,j) ~= max_comp
                        result(i,j) = 0;
                    end
                end
            end
        end
        
        % theory & code based on 
        % https://wp.optics.arizona.edu/visualopticslab/wp-content/uploads/sites/52/2016/08/Lectures6_7.pdf
        function result = projective_transform(~, image, sudo_X, sudo_Y, image_X, image_Y)
            [rows,cols,channel] = size(image);
            
            % create A (8x8), similar to our project 1 but with 4 points
            A = [sudo_X, sudo_Y, ones(size(sudo_X)), zeros(4,3), -sudo_X.*image_X, -sudo_Y.*image_X...
                zeros(4,3), sudo_X, sudo_Y, ones(size(sudo_X)), -sudo_X.*image_Y, -sudo_Y.*image_Y];
            A = reshape(A', 8, 8)';

            % B is column vector of coords we are mapping to
            B = [image_X, image_Y];
            B = reshape(B', 8, 1);
            
            % A*lambda = B
            % lambda = (A^T * A)^-1 * A^T * B
            lambda = inv(A' * A) * A' * B;
            
            % reshape the result into 'H' from the lecture slides
            H = transpose(reshape([lambda(1:8)' 1], 3, 3)');
            inv_H = inv(H);
            
            result = zeros(rows,cols);
            
            for j = 1:rows
                for i = 1:cols
                    final_loc = [i,j,1];
                    
                    % if X' = HX then 
                    % X = H^-1 * X'
                    % calculate which original pixel maps to the new image
                    % pixel
                    starting_loc_x = round((inv_H(1,1)*i + inv_H(2,1)*j + inv_H(3,1)) / (inv_H(1,3)*i + inv_H(2,3)*j + inv_H(3,3)));
                    starting_loc_y = round((inv_H(1,2)*i + inv_H(2,2)*j + inv_H(3,2)) / (inv_H(1,3)*i + inv_H(2,3)*j + inv_H(3,3)));
                    result(j,i) = image(min(rows,starting_loc_y), min(cols,starting_loc_x));
%                     result(j,i) = image(min(cols,starting_loc_x), min(rows,starting_loc_y));

                end
            end
            
            result = uint8(result);
        end
        
        function [sudo_X, sudo_Y, image_X, image_Y] = detect_corners(~, image)
            [rows,cols,channel] = size(image);
            top_left = [0,0];
            top_left_sum = 99999;

            bottom_left = [0,0];
            bottom_left_sum = 99999;

            top_right = [0,0];
            top_right_sum = 0;

            bottom_right = [0,0];
            bottom_right_sum = 0;
            % for every cell in image
            for i = 1:rows
                for j = 1:cols
                    % if it is a corner
                    if image(i,j) == 1
                        % check if it is one of the four extreme corners
                        % in the sudoku and keep track
                        if i+j < top_left_sum
                           top_left = [j,i]; 
                           top_left_sum = i+j;
                        end

                        if i-j > top_right_sum
                            top_right = [j,i];
                            top_right_sum = i-j;
                        end

                        if i+j > bottom_right_sum
                           bottom_right = [j,i];
                           bottom_right_sum = i+j;
                        end

                        if i-j < bottom_left_sum
                           bottom_left = [j,i];
                           bottom_left_sum = i-j;
                        end
                    end
                end
            end
            figure; imshow(image); hold on; title('Sudoku Corners');

            plot(bottom_left(1,1),bottom_left(1,2), 'g*');
            plot(top_right(1,1),top_right(1,2), 'g*');

            plot(top_left(1,1),top_left(1,2), 'g*');
            plot(bottom_right(1,1),bottom_right(1,2), 'g*');
            
            global savedir
            global testname
            filename = strcat(testname, '_corners.png');
            saveas(gcf, fullfile(savedir, filename));
            % store the sudoku xy and image xy in a format to be used by
            % the projective transform
            sudo_X = [top_left(1,1);bottom_left(1,1);top_right(1,1);bottom_right(1,1)];
            sudo_Y = [top_left(1,2);bottom_left(1,2);top_right(1,2);bottom_right(1,2)];
            image_X = [1; cols; 1; cols];
            image_Y = [1; 1; rows; rows];
        end
    end
end
