close all;
clear all;

% paramters
imThresh = 128; % value to threshold images
centFactor = 100; % twice this value is the centeral region selected

% calibration parameters from K-readings spreadsheet
slope = 14.604; % slope from calibration
intercept = 11.467;

% debug params
plot_fig = 1; % 0 = do not do debug plots

% factor - what portion of the central image to use
factor = 64;

% set working directory for image files
wd = '/media/ruckerlab/tiny/KERA/Falk_Matt_Kera'; %% edit here %%

% set the data dir
datadir = '/media/ruckerlab/tiny/KERA/Falk_Matt_Kera';

% edit with name of .csv logfile -- HERE

% create a raw datafile and header
big_filename = 'falk_kera_all.csv';  %% edit here %%

light_d = {'white' , 'yellow', 'yb', 'yk'}; %% edit here -- condition1 %%
freq_d = {'scgx'}; %% edit here -- condition2 %%
pre_post_str = {'PRE', 'POST'};
eye = {'OD','OS'}


kera_file = fopen(big_filename,'w+'); 
% write header
% bird# pre/post, light, freq
fprintf(kera_file,'bird,condition,light,meas,eye,m,j0,j45\n');

% loop through all dirs with image files
for brCtr = 1:length(light_d)
    for freqCtr = 1:length(freq_d)
        for pre_post = 1:length(pre_post_str)
            % change to the wd
            cd(sprintf('%s/%s/%s', wd, freq_d{freqCtr}, light_d{brCtr}));
            fprintf('dir: %s\n', pwd);
            if pre_post == 1    
                % get pre files
                im_dir= dir('*PRE*');
            elseif pre_post == 2
                % get post files
                im_dir = dir('*POST*');
            elseif pre_post == 3
                % get rec files
                im_dir = dir('*REC*');
            end
            % loop through pre/post, data mat, bird #, light, freq, pre/post, M, J0, J45

            birdList = im_dir;

            % pre-allocate wh var
            wh = []; tmp =[];

            counter = 1;
            pcounter = 1;
              
            for imCtr = 3:size(im_dir)

                fprintf('bird: %s\n',birdList(imCtr).name(1:4));
                eye_str = birdList(imCtr).name(6:7);    
                % read in image file
                im = imread(birdList(imCtr).name);
                im = double(rgb2gray(im));

                % convert to thresholded logical matrix
                imlv = (im >= imThresh);

                % get image size    
                sz = size(imlv);
                cnt_x = sz(1)/2;
                cnt_y = sz(2)/2;

                % get the center of the image
                cnt_imlv = imlv(cnt_x-factor:cnt_x+factor,cnt_y-factor:cnt_y+factor);

                % get the region properties using centroid method
                s  = regionprops(cnt_imlv, 'centroid');
                centroids = cat(1, s.Centroid);

                if ~isempty(centroids)
                    wh = [max(centroids(:,1)) - min(centroids(:,1)) ...
                    max(centroids(:,2)) - min(centroids(:,2))];
                end
                %centroids = centroids -  (wh'/2 * ones(1,size(centroids,1)));
                ocentroids = centroids;
                for tt = 1:size(centroids,1)
                    centroids(tt,:) = centroids(tt,:) - wh/2;
                end

                [A, c] = MinVolEllipse(centroids', 0.2);
                [Ao, co] = MinVolEllipse(ocentroids', 0.2);

                % test code to plot region finding - debug flag at top
                if plot_fig == 1
                    %figure(counter+100);
                    %imagesc(im);
                    %colormap('gray')
                    h = figure(pcounter);
                    %imshow(cnt_imlv);        
                    imagesc(im(cnt_x-factor:cnt_x+factor,cnt_y-factor:cnt_y+factor));
                    colormap('gray')
                    hold on;        
                    plot(ocentroids(:,1), ocentroids(:,2), 'g*');
                    Ellipse_plot(Ao,co);
                    %saveas(h,sprintf('marked_%s',birdList(imCtr).name));
                    %cd('..')
                    pcounter = pcounter + 1;
                end;

                % "singular value decomposition" to extract the orientation and the
                [U D V] = svd(A);
                
                % get the major and minor axes lengths
                a(counter) = 1/sqrt(D(1,1)); 
                b(counter) = 1/sqrt(D(2,2));

                % conver to mm;
                a(counter) = (a(counter) - intercept)./slope; 
                b(counter) = (b(counter) - intercept)./slope;

                % which axis is major / minor, compute D
                if a(counter) > b(counter)
                    S(counter) = 1/a(counter) / 1.3379;
                    J(counter) = 1/(b(counter) - a(counter)/2) / 1.3379;
                else
                    S(counter) = 1/b(counter) / 1.3379;    
                    J(counter) = 1/(a(counter) - b(counter)/2) / 1.3379;
                end    

                % get angle of the ellipse    
                Ax(counter) =  atan(V(2,1)/V(1,1))*180/pi;

                % convert to power vector notiation
                M(counter) = S(counter) + J(counter);
                J0(counter) = J(counter) .* cos(2.*(deg2rad(Ax(counter) - 90)));
                J45(counter) = J(counter) .* sin(2.*(deg2rad(Ax(counter) - 90)));
                bName(counter) = str2num(birdList(imCtr).name(1:4));
                
                % get current analysis dir
                oldDir = pwd;
                cd(wd);
                % write bird#, light, condition, pre/post/rec, M, J0, J45
                % strings
                fprintf(kera_file,sprintf('%s,%s,%s,%s,%s,%f,%f,%f\n', num2str(bName(counter)),...
                    freq_d{freqCtr}, light_d{brCtr},pre_post_str{pre_post},...
                    eye_str, M(counter), J0(counter), J45(counter)));    
                % increment counter
                counter = counter + 1;
                % go back to analysis dir
                cd(oldDir);
                which_eye = 2;
            end
            % save M, J0, J4 in .csv files, pre/post, light, freq
            data = [bName' M' J0' J45'];
            oldDir = pwd;
            cd(datadir);
            eval(sprintf('csvwrite(''%s_%s_%s_kera.csv'',data);', pre_post_str{pre_post},light_d{brCtr}, freq_d{freqCtr}));
            cd(oldDir);
        end    
    end
end
fclose(kera_file);

% %% correlate pre-measures and save
% zgo('autoKera');
% ctr = 1;
% for brCtr = 1:length(light_d)
%     for freqCtr = 1:length(freq_d)
%         eval(sprintf('dPRE = csvread(''PRE_%s_%s_kera.csv'');', light_d{brCtr}, freq_d{freqCtr}));
%         for avgCTR = 2:2:size(dPRE,1)
%             dPRE_N(ctr,:) = dPRE(avgCTR,2:4);
%             dPRE_X(ctr,:) = dPRE(avgCTR-1,2:4);
%             ctr = ctr + 1;
%         end
%     end
% end
% figure;
% plot(dPRE_N(:,1),dPRE_X(:,1),'ks');
% title('M');
% figure;
% plot(dPRE_N(:,2),dPRE_X(:,2),'ks');
% title('J0');
% figure;
% plot(dPRE_N(:,3),dPRE_X(:,3),'ks');
% title('J45');

%% post processing
% loop through B/R, all frequencies
cd(datadir);
for brCtr = 1:length(light_d)
    for freqCtr = 1:length(freq_d)
        clear('dPRE','dPOST','Bird','dTMP','dPREavg','dPOSTavg');
        eval(sprintf('dPRE = csvread(''PRE_%s_%s_kera.csv'');', light_d{brCtr}, freq_d{freqCtr}));
        eval(sprintf('dPOST = csvread(''POST_%s_%s_kera.csv'');', light_d{brCtr}, freq_d{freqCtr}));        
        ctr = 1;
        for avgCTR = 2:2:size(dPRE,1)
            dPREavg(ctr,:) = (dPRE(avgCTR,2:4) + dPRE(avgCTR-1,2:4))/2;
            dPOSTavg(ctr,:) = (dPOST(avgCTR,2:4) + dPOST(avgCTR-1,2:4))/2;
            Bird(ctr) = dPRE(avgCTR);
            Light(ctr) = brCtr;
            Freq(ctr) = freqCtr;
            ctr = ctr + 1;
        end
        dDiff = dPOSTavg - dPREavg;
        [r,c] = find(dDiff == 0);
        if ~isempty(r)
            dDiff(r(1)-1:end,:) = [];
            Bird(r(1)-1:end) = [];
            Light(r(1)-1:end) = [];
            Freq(r(1)-1:end) = [];
        end
        dTMP = [Bird' Light' Freq' dDiff]
        eval(sprintf('csvwrite(''%s_%s_kera.csv'',dTMP);', light_d{brCtr}, freq_d{freqCtr}));
    end
end