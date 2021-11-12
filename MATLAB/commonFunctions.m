%% Common scripts used throughout modules
classdef commonFunctions
	methods (Static)
        % Loading data - output will always be x-y-frame-intensity
        % Load data from CSV file (ThunderSTORM)
        function output = readCSVThunderSTORM(filelocation)
            %Read the input file as a table (so we keep header information)
            output = readtable(filelocation,'VariableNamingRule','Preserve');
            %Select which headers to keep
            columnTitles = {'frame','x [nm]','y [nm]','intensity [photon]'};
            %Find the position of these headers
            columnIndeces = zeros(1,size(columnTitles,2));
            for i = size(columnTitles,2):-1:1
                try
                    columnIndeces(i) = find(strcmp(output.Properties.VariableNames,columnTitles{i}));
                catch
                    %Catch if the column index doesn't exist
                    columnIndeces(i) = 0;
                end
            end
            if min(columnIndeces) > 0
                %Select only the headers we want, and output as an array
                output = table2array(output(:,columnIndeces));
            else %If not all column indeces existed
                %Create a new output of correct size
                outputnew = array2table(zeros(size(output,1),size(columnIndeces,2)));
                %Loop over the columnIndeces, and fill if it's non-zero, or
                %leave as zeros if it isn't included
                for i = 1:size(columnIndeces,2)
                    if columnIndeces(i) > 0
                        outputnew(:,i) = output(:,columnIndeces(i));
                    end
                end
                %output as an array
                output = table2array(outputnew);
            end
        end
        
        % Loading data - output will always be x-y-frame-intensity
        % Load data from TXT file (RapidSTORM)
        function data = readTXTRapidSTORM(filelocation)
            dataunordered = readmatrix(filelocation);

            fid = fopen(filelocation);
            header = convertCharsToStrings(fgets(fid));
            fclose(fid);
            clear fid

            %Extract data from the header
            fieldstartindeces = strfind(header,'<field');
            nrcols = size(fieldstartindeces,2);
            %find all spaces
            allspaces = strfind(header,' ');

            columnorder = zeros(1,4);
            %Get the info from the fields
            for i = 1:nrcols
                %find first index of a space
                spacelist = (allspaces(allspaces>fieldstartindeces(i)));
                nextspace = spacelist(2);
                %find info
                offset = 19;
                infostring = extractBetween(header,fieldstartindeces(i)+offset,nextspace-2);
                %Based on the information in this infostring, we change the columnorder
                switch infostring
                    case 'Position-0-0' %x
                        columnorder(2) = i;
                    case 'Position-1-0' %y
                        columnorder(3) = i;
                    case 'ImageNumber-0-0' %frame
                        columnorder(1) = i;
                    case 'Amplitude-0-0' %intensity
                        columnorder(4) = i;
                end  
            end

            %And we order the data based on the obtained columnorder
            data = dataunordered(:,columnorder);
            data(:,1) = data(:,1)+1; %start frame at 1 instead of 0
        end
        
%         function [tracked_data] = nearestNeighbourFinding(data)
%             % This function is used in modules 4a, 7, 8, and 9a. It adds
%             % the nearest neighbour ID (i.e. the row index) and distance to
%             % each datapoint.
%             try
%             % The list of frames is put into uniqueFrame variable in parallel with the number of events per each frame.
%             % The result of this is simply an array of all frames that are encountered in the localization list, as well as an array with how many localizations are in all these frames. These arrays are used later.
%             [uniqueFrame, cumFrameLocCounts] = unique(data(:,1));
%             % Two empty columns are added to the dataset for later use
%             twoColumns = zeros(size(data,1),2);
%             tracked_data = [zeros(size(data)) twoColumns];
% 
%             %We use a counter to fill tracked_data
%             counter = 1;
% 
%             % Dataset is processed in a for loop, frame-by-frame to improve the computation speed.
%             for frame = min(uniqueFrame):max(uniqueFrame)-1
%                 % First, the data from the current and the next frame are stored in separate variables.
%                 current_frame = data(data(:,1)==frame,:);
%                 next_frame = data(data(:,1)==frame+1,:);
%                 
%                 % Only continue if both frames have at least 1 localization
%                 if size(current_frame,1)>0 && size(next_frame,1)>0
%                     % Now we can find all the nearest neighbours between all localizations on
%                     % this frame and on the next frame
%                     [neighbourID, neighbourDist] = knnsearch(next_frame(:,2:3),current_frame(:,2:3));
%                     % The index of the nearest next-rame neighbor is saved in the current frame data variable together with the distance between both datapoints. Pay attention that we are adding the cumulative count of events so the index will be correct once the dataset will be fully processed.
%                     current_frame(:,5) = neighbourID + cumFrameLocCounts(find(uniqueFrame == frame)+1) - 1;
%                     current_frame(:,6) = neighbourDist;
% 
%                     %Now we have to remove duplicate links in case a single
%                     %emitter is linked to >1 emitter in the next frame
%                     %We loop over the unique track ids in this frame
%                     uniqueList = unique(current_frame(:,5));
%                     for k = 1:size(uniqueList,1)
%                         t = uniqueList(k);
%                         %We check how many there are, and if there are more
%                         %than 1, we continue
%                         current_linkages = find(current_frame(:,5)==t);
%                         if size(current_linkages,1)>1
%                             %We find the index of the minimum track length,
%                             %because we do not want to remove this one
%                             index_min_trackLength = find(current_frame(current_linkages,6)==min(current_frame(current_linkages,6)));
%                             %Next, we find which linkage indexes need to be
%                             %removed (i.e. the others)
%                             linkages_to_be_removed = current_linkages(setdiff(1:end,index_min_trackLength));
%                             %And we set those to zero
%                             current_frame(linkages_to_be_removed,[5,6])=zeros(size(linkages_to_be_removed,1),2);
%                         end
%                     end
%                     % Finally, the processed frame is stored into a new variable which contains six columns: X-coordinate, Y-coordinate, time coordinate, intensity, the index of the nearest next-frame neighbor and the distance to it.
%                     tracked_data(counter:counter+size(current_frame,1)-1,:) = current_frame;
%                     counter = counter+size(current_frame,1);
%                 else %If one of the frames does not have localizations, still add the current_frame to the tracked_data if it has localisation
%                     if size(current_frame,1)>0
%                         tracked_data(counter:counter+size(current_frame,1),:) = current_frame;
%                     end
%                 end
%             end
%             % The last frame processed data are directly appended at the end of the array as they are not processed.
%             tracked_data(counter:end,1:4) = next_frame;
%             catch
%                 keyboard
%             end
%         end
function [tracked_data,tracksCounter] = creating_trajectories(tracked_data,frame1,frame2,maxDistance,tracksCounter)
            %Make sub-matrix of this frame and the next frame
            framematrix = tracked_data(tracked_data(:,1)==frame1,:);
            nextframematrix = tracked_data(tracked_data(:,1)==frame2,:);
            
            %Now we check that there are localizations in both of these frames - if one
            %of the frames does not have localizations, we cannot try to track the
            %localizations. Note that we already know that in this case, frame 0 and
            %frame 1 have localizations, but this will not always be the case.
            if (size(framematrix,1) > 0 && size(nextframematrix,1) > 0)
                %Now we can find all the nearest neighbours between all localizations
                % on this frame and on the next frame
                [neighbourID, neighbourDist] = knnsearch(nextframematrix(:,2:3),framematrix(:,2:3));
                % We store this data in a 2-column array
                foundnn = [neighbourID neighbourDist];
            
                % Now we have to select only those entries that are < maxDistance. We want
                % to extract the IDs of the localizations in frame n+1 that belong to them.
                NeighbourIDs = foundnn(foundnn(:,2) < maxDistance,1);
                % We also want to extract the ID of the original localizations by
                % finding the row index of the nearest neighbours.
                OriginIDs = find(foundnn(:,2) < maxDistance);
            
                % For every found neighbour, we make both neighbours the same track-id.
                % First, we check that the neighbour on frame n has or doesn't have an
                % track-id yet, then we set the neighbour on frame n+1
                % to the same value, or to a new track-id if none is assigned yet
                % We loop over all found IDs
                for i = 1:size(NeighbourIDs,1)
                    %We get the localization-ID of the neighbour in frame n+1
                    neighbourID = NeighbourIDs(i);
                    %We also get the localization-ID of the original localization in frame n
                    originID = OriginIDs(i);
                    %Prevent linkage if the neighbour in frame n+1 is already linked
                    if nextframematrix(neighbourID,5) == 0
                        %We check that the localization that will be included in a trajectory
                        %is not yet part of an existing trajectory - if it is part of
                        %an existing track, the index in column 5 will be higher than 0
                        if framematrix(originID,5) == 0
                            %If it's not linked yet, set it and the neighbour to a new track-id value
                            framematrix(originID,5) = tracksCounter;
                            nextframematrix(neighbourID,5) = tracksCounter;
                            tracksCounter = tracksCounter+1;
                        else
                            %If it is linked, set it to the track-id value of the neighbour in frame n
                            nextframematrix(neighbourID,5) = framematrix(originID,5);
                        end
                        %Finally, we  provide the distance to the next emitter in the
                        %track on the next frame
                        framematrix(originID,6) = foundnn(originID,2);
                    end
                end
            end
            %Now we have fully filled frame matrix and nextframematrix variables, but
            %we need to fill these back in into the original tracked_data matrix. We do
            %this by looking up the values via the idCounter value
            if size(framematrix,1) > 0
                tracked_data(tracked_data(:,1)==frame1,:) = framematrix;
            end
            if size(nextframematrix,1) > 0
                tracked_data(tracked_data(:,1)==frame2,:) = nextframematrix;
            end
        end
        
        function binim = TwoDimBinning(localizations,pxsize,histSize)
            %Check if histSize is used
            if exist('histSize','var')
                numberBins = histSize;
            else
                %We specify the number of bins based on the px size
                numberBins = [ceil(max(localizations(:,2))/pxsize), ceil(max(localizations(:,3))/pxsize)];
            end
            %And the binning itself is done via the histcounts2 function.
            binim = histcounts2(localizations(:,2),localizations(:,3),numberBins,'XBinLimits',[0, numberBins(1)*pxsize],'YBinLimits',[0, numberBins(2)*pxsize]);
        end
	end
end
%%
