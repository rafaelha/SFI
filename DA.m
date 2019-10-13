% DA class stands for data analysis done for NO experiments at Prof.
% Grant's laboratory. This class call essential functions that loads data
% from a selected directory. There are some functions that facilitate 
% plotting and some other tasks.

classdef DA
    methods (Static)
        function d=readspec(path)
            % readspec take path as input. The path should indicate the
            % location of the folder that contains pfi, w2, or w1 spectrum
            
            % find files that end with .wf
            f=dir(fullfile(path,'*.wf'));
            
            % initiate a data structure
            d=struct;
            
            % loop over all .wf files
            for i=1:length(f)
                % print the file name in the command window
                fprintf('%s\n',f(i).name);
                
                % read dataset from a set of files
                % the inputs are the main path to the folder and the name
                % of each file with .wf in the folder path
                % the outputs are 2D data, specifications of the waveform,
                % time array, and the scanning variable (delay, wavelength,
                % ...)
                [d(i).data,d(i).wf,d(i).t,d(i).wl]=...
                    DAIO.readdat(path,f(i).name);
                
                % add a new field called name which is the same as the file
                % name except the date at the beginning and .wf at the end
                d(i).name=f(i).name(13:end-3);
            end
        end
        function d=readPS(path, lines)
            % readPS take path and lines as input and returns pulse shape
            % data. Path is the folder location of the pulse shapes and
            % lines is the number of comment lines on top of the data. The
            % code will eliminate the top part of the file as specified by
            % number of lines and reads the rest as data
            
            % find files that end with .wf
            f=dir(fullfile(path,'*.wf'));
            
            % initiate a data structure
            d=struct;
            for i=1:length(f)
                % readwl is the core function that read a pulse shape file.
                % it take the original path and the name of each pulseshape
                % in the folder and the number of lines of comments and
                % processes the data
                [d(i).t,d(i).v,d(i).wf]=DAIO.readwl(path,f(i).name, lines);
                
                % make sure that there is no offset to the voltage
                d(i).v=d(i).v-d(i).v(1);
                
                % add a new field name from the file name
                d(i).name=f(i).name;
                
                % EF (electric field) used to be voltage divided by
                % distance. Right now EF is redundant as it is the same as
                % v. Let's keep it that way for now.
                d(i).EF=d(i).v;
                
                % fit a smooth function to the shape of the electric field
                % without worrying much about its functional form
                d(i).cf=fit(d(i).t(:),sort(d(i).EF(:)),...
                    fittype('smoothingspline'));
            end
        end
        function x=fitns(wl,g)
            % this function take a set of wavelength at maxiumum absorbtion
            % and an initial guess of g for the n0 = quantum number and
            % delta = quantum defect and attempts to fit IP, delta, and
            % quantum number
            % wavelengths must be consequtive. It can be ascending or
            % descending. Wavelength is of the doubled frequency and vacuum
            % corrected
            
            % no matter what or wl is, make is ascending
            wl=sort(wl);
            
            % compute wavenumber from double frequency wavelength in vacuum
            wavenum=1e7./wl; % input is in nm. output is in cm^-1
            
            % Rydberg constant cm^-1
            Ryd=109737.316;
            
            % an objective function that returns the error of fit
            function o=obj(x)
                % quantum defect
                delta=x(3);
                
                % build a range of numbers that when added to n0 = g(1)
                % will give the range of quantum numbers.
                n=length(wl)-1:-1:0;
                
                % convert the range of numbers to range of quantum numbers
                n=n+fix(x(2));
                
                % compute difference of experimental wavenumbers and the
                % computed one
                o=wavenum'-(x(1)-Ryd./(n-delta).^2);
            end
            
            % minimize error of the objective function to find the best fit
            %                          IP             n0    delta
            % x is 3 by 1 array: IP, n0, and delta
            x=lsqnonlin(@(x)obj(x),[30544.2333070219,g(1),g(2)]);
            
        end
        function x=fitnsqdfixed(wl,g)
            % this function take a set of wavelength at maxiumum absorbtion
            % and an initial guess of g for the n0 = quantum number and
            % delta = quantum defect and attempts to fit IP, delta, and
            % quantum number
            % wavelengths must be consequtive. It can be ascending or
            % descending. Wavelength is of the doubled frequency and vacuum
            % corrected
            
            % no matter what or wl is, make is ascending
            wl=sort(wl);
            
            % compute wavenumber from double frequency wavelength in vacuum
            wavenum=1e7./wl; % input is in nm. output is in cm^-1
            
            % Rydberg constant cm^-1
            Ryd=109737.316;
            
            % quantum defect
            delta=g(2);
            
            
            % an objective function that returns the error of fit
            function o=obj(x)
                
                % build a range of numbers that when added to n0 = g(1)
                % will give the range of quantum numbers.
                n=length(wl)-1:-1:0;
                
                % convert the range of numbers to range of quantum numbers
                n=n+fix(x(2));
                
                % compute difference of experimental wavenumbers and the
                % computed one
                o=wavenum'-(x(1)-Ryd./(n-delta).^2);
            end
            
            % minimize error of the objective function to find the best fit
            %                          IP             n0    delta
            % x is 3 by 1 array: IP, n0, and delta
            x=lsqnonlin(@(x)obj(x),[30544.2333070219,g(1)]);
            
        end
        function wl = ntowl(ns,x)
            % compute the wavelength of the double frequency light in
            % vacuume (make sure multiply by 2 if you're working with the
            % fundamental wavelength out of ND6000
            wl = 1e7./(x(1)-cons.Ryd/100./(ns-x(2)).^2);
        end
    end
end