classdef DAIO
    %DAIO Summary of this class goes here
    %   DAIO stands for data analysis input output
    % this class contains core functions that are called by DA class. It is
    % seldom necessary to call DAIO directly from your script
    
    methods (Static=true)
        
        function [dat,wf,t,wl] = readdat(path, fname)
            % this function read the .dat file containing the 2D intensity
            % plot
            
            % read the .wf files to get the parameters of the waveform
            % filename should be imported as something .wf (handled by DA
            % class)
            wf = DAIO.readwf(path, fname);
            
            % the 3rd parameter from .wf file is the # of points in a
            % single trace, same as the length of the time array
            num_points = wf(3);
            
            % now keep the same file name but swap the extension. The new
            % filename will look for a .dat file with the same name as .wf
            fname = regexprep(fname, '.wf', '.dat');
            
            % the full file name is path + name
            file = [path, fname];
            
            % open a file as read-only and return an ID to the file handle
            fid = fopen(file);
            
            % read from the file with the ID = fid. Inf specifies the
            % number of bytes to read, in this case everything. 'uint16'
            % specifies the format of the binary data which is unsigned 16
            % bit integer (the number of points on the y axes of the scope
            % starting from zero. 'b' stands for bigendian. Labview uses
            % IEEE standard to record numbers as binary, we tell matlab how
            % to read that number as an integer
            rawdat = fread(fid, inf, 'uint16', 'b');
            
            % close the connection to the file (a good practice to always
            % do that)
            fclose(fid);
            
            % rawdata will be a long array of integers. The processdat
            % function in this class will reshape and scale the data to
            % give voltages in mV. I have eliminated the need to offset
            % data in the core function.
            dat=DAIO.processdat(rawdat,num_points,wf);
            
            % build time array based on the number of points, division per
            % point and the offset
            t=DAIO.buildt(wf,num_points);
            
            % this function will read the wl file associated with other
            % files that have the same name but ends in .wl
            wl=DAIO.wlfileuni(path,fname);
            
            % make sure all the connection to files are closed
            fclose('all');
        end
        function t=buildt(wf,num_points)
            t0=wf(6);   % initial scope time
            tinc=wf(5); % time increment
            tf=tinc*num_points+t0; % final time
            
            % generate an array with specified start, end and number of
            % points
            t=1e6*linspace(t0,tf,num_points);
        end
        function dat=processdat(rawdat,num_points,wf)
            % reshape the rawdata based on the number of point in each
            % single trace
            dat = reshape(rawdat, num_points, max(size(rawdat))/...
                num_points); 
            %just reshape it
            
            % convert integers to double
            dat = double(dat'); %convert to double
            
            % convert raw values to mV
            dat=-dat*wf(8); %convert to mV
        end
        
        function wf=readwf(path,fname)
            % this function will read the waveform file that containst some
            % information about the
            
            % build full path to the file
            fname = [path fname];
            
            % open file as read-only
            fID=fopen(fname);
            
            % read one unit of binary data in the format specified as 'b' =
            % bigendian. The unit will depend on the format. int16 wil read
            % 2 bytes, while double will read 8 bytes
            wf(1)=double(fread(fID,1,'int16=>int16','b'));
            wf(2)=fread(fID,1,'int16=>int16','b');
            wf(3)=fread(fID,1,'int32=>int32','b');
            wf(4)=fread(fID,1,'int32=>int32','b');
            wf(5)=fread(fID,1,'double=>double','b');
            wf(6)=fread(fID,1,'double=>double','b');
            wf(7)=fread(fID,1,'int32=>int32','b');
            wf(8)=fread(fID,1,'double=>double','b');
            wf(9)=fread(fID,1,'double=>double','b');
            wf(10)=fread(fID,1,'int32=>int32','b');
            
            % close connection to all files
            fclose('all');
        end
        function [t,d,wf]=readwl(path,fname,lines)
            fname2=regexprep(fname,'.wf','.wl');
            d=importdata([path fname2],'\t',lines);
            wf=DAIO.readwf(path,fname)';
            d=d.data;
            t=(d(:,1)'+wf(6))*1e6;
            d=d(:,2)'*wf(8);
            fclose('all');
        end
        function [t,d,wf]=readwlraw(path,fname,lines)
            fname2=regexprep(fname,'.wf','.wl');
            d=importdata([path fname2],'\t',lines);
            wf=DAIO.readwf(path,fname)';
            d=d.data;
            t=(d(:,1)'+wf(6))*1e6;
            d=d(:,2)';
            fclose('all');
        end
        function wl=wlfileuni(path,fname)
            % this function does the same thing as
            % readwl. However, it will not require you to specify the
            % number of lines of commnets
            
            % prepare the name of the file with the same name as the .dat
            % file but with .wl extension
            fname2=regexprep(fname,'.dat','.wl');
            %             fname2=[fname(1:end-3) '.wl'];
            
            % open the file with read-only permission
            fid=fopen([path fname2]);
            %             [path fname2]
            % get the first line
            tline = fgetl(fid);
            
            % while the first line or the next one are not #START continue
            % reading lines
            while ~strcmp(tline,'#START')
                tline=fgetl(fid);
                if tline==-1
                    break;
                end
            end
            
            q=1;
            
            % the curser to the fid must be right after #START, where
            % comments end and data begins
            while ~feof(fid)
                % read a new line
                tline=fgetl(fid);
                % split the line by \t = tab between the numbers
                tline=strsplit(tline,'\t');
                % convert the array of text cells to an array of doubles
                wl(q,:)=str2double(tline);
                
                % advance the index
                q=q+1;
            end
            % close connection to all files
            fclose(fid);
            fclose('all');
        end
    end
    
end

