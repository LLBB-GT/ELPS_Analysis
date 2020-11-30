% Read the video file
vid_name = uigetfile();
% Create a videoreader object
cont_vid = VideoReader(vid_name);

% Read and display the first frame
I_init = readFrame(cont_vid);
I_init = rgb2gray(I_init);
figure
imshow(I_init);
% Select the vessel from the first frame
init_rect = getrect;
close

[~, col] = size(I_init);
box_width = 10;   % Width of the window for binarizing the frame
n_pts = floor(col/box_width);   % Total number of binarizing windows
% Initialize the variables for calculating and storing the diameter data
dia_up_mat = [];
dia_down_mat = [];
dia_mat = [];

% Binarize and calculate the diameter for each frame
while hasFrame(cont_vid)
    
    I = readFrame(cont_vid);
    I = rgb2gray(I);   % Convert the frame to grayscale
    windows = cell(1,n_pts);
    bw_windows = cell(1,n_pts);
    bw_img = [];
    for i = 0:n_pts-1
        wind_rect = [i*box_width+1,init_rect(2),box_width-1,init_rect(4)];
        windows{1,i+1} = imcrop(I,wind_rect);   % Divide the frame into windows
        thresh = graythresh(windows{1,i+1});   % Calculate the hreshold
        bw_windows{1,i+1} = im2bw(windows{1,i+1},thresh);   % Binarize the frame
        bw_img = [bw_img bw_windows{1,i+1}];   % Stich the binarized images together
    end

    [~, s_col] = size(bw_img);
    dia_up = zeros(1,s_col);
    dia_down = zeros(1,s_col);
    
    % Calculate the upper and lower wall of the binarized image
    for line = 1:s_col
        dia_up(1,line) = find(~bw_img(:,line),1,'first');
        dia_down(1,line) = find(~bw_img(:,line),1,'last');
    end
    
    dia_up_mat = [dia_up_mat; dia_up];
    dia_down_mat = [dia_down_mat; dia_down];
    
end

% Calculate the diameter
dia_mat = dia_down_mat - dia_up_mat;

% Save the diameter file
save(strrep(vid_name,'.avi','_dia_mat'),...
    'dia_up_mat','dia_down_mat','dia_mat','-v7.3')
