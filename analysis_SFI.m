%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Code for analysis of SFI spectra  %%%%%%%%%%
%%%%%%%%%                                    %%%%%%%%%%
%%%%%%%%%      05/10/2017 Rafael Haenel      %%%%%%%%%%
%%%%%%%%%                                    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc; %clear Workspace and console
close all; %close all windows
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  !The following values must be updated!      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

originaldate='10-May-2017';

%specify working directory
folder_path='C:\Users\rafae\Documents\Thesis UCP\Data\2017\03mar\27\'; 
%the folder must include 
%   (1) a folder named 'ps' that includes the pulse-shape files
%   (2) a second folder named 'sfi'
%       (A) The 'sfi' folder includes a number of folders for different ramp
%       delays. They must be named in the following way: 'RampDelay2.0u',
%       'RampDelay100.0n' (1 decimal place!)
%           (a) Each ramp delay folder must include folders for all w1w2
%           dealys. They must be named in the following way: 'w1w2Delay7n',
%           'w1w2Delay1u'

%set the list of dw1w2 dealays. The experiment should be set up to return
%equally spaced data points. The following values have been calculated to achieve that: [0 5 7 15 23 31 40 49 59 69 79 90 102 114 127 141 156 172 190 209 230 254 281 312 348 391 447 524 647 800 1000]
dw1w2=[0 5 7 15 23 31 40 49 59 69 79 90 102 114 127 141 156 172 190 209 230 254 281 312 348 391 447 524 647 800 1000];

lifetime = 206; %lifetime of first excited state in ns

PQN=49; %principle quantum number

nt=100; %number of traces per w1w2-delay-measurement (as set in the LabView program)

grid_position=110; %gridposition in mm

x1=1730; %x1 and x2 specify the interval (in pixels) of the trace that should be considered. 
x2=6800-1392; %All pixels from x1 to x2 will be read into the matrix 'traces'

bg1=1; % set background noise range (with respect to x1)
bg2=1000; %e.g. bg1=1 and bg2=1000 means that the first 1000 pixels will be interpreted as background-noise

n_bins=35; %number of horizontal bins (pixels) of the SFI spectrum

max_den=2.0; %maximum rel. density (upper y-limit of the SFI-plot)

maxcolor=3.0; %color limit for SFI false color scale

separation =350; %the location x that separates plasma signal and Rydberg signal. All pixels left of 
%separation are interpreted as plasma signal, all pixels to the right are
%Rydberg signal


%load ps file and evaluate ramp
dist=(155.5-grid_position)*1e-1;
shift=0.042; %the electrons need 42ns to reach the detector
ps=DA.readPS([folder_path, 'ps\'], 19); %read ps (skip the first 19 lines)
EF=ps.cf(ps.t(x1:x2)-shift)/dist; %calculated EM field

%find all ramp delay folders
folder=dir([folder_path,'sfi combined\']);
G=zeros(1,length(folder)-2);
G_folder=cell(1,length(folder)-2);
for s=1:(length(folder)-2)
    f_name=folder(s+2).name;
    f_num=replace(f_name,'RampDelay','');
    f_num=replace(f_num,'n','');
    f_num=replace(f_num,'u','*1000');
    G(s)=str2num(f_num);
    G_folder{s}=f_name;
end
[G, f_order]=sort(G);


den=1*exp(-dw1w2/lifetime); %rel. density values

nG=length(G); %number of gridpositions
nD=length(dw1w2); %number of w1w2-delays

%Preallocate a bunch of arrays
traces=zeros(nD*nt,x2-x1+1); %preallocate the SFI-matrix
av_sig=zeros(1,nD);
std_err_sig=zeros(1,nD);

peaks=[];

ev=zeros(n_bins+1,x2-x1+1,length(G)); %tensor of 3rd order that contains time evolution of SFI spectra

plasma_signal=zeros(n_bins+1,nG); %average plasma signal for all ramp delays
plasma_error=zeros(n_bins+1,nG);
rydberg_signal=zeros(n_bins+1,nG); %average Rydberg signal for all ramp delays
rydberg_error=zeros(n_bins+1,nG);

for i=1:nG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%          Load all data                %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:nD
        path=strcat([folder_path, 'sfi combined\',G_folder{f_order(i)},'\w1w2Delay',num2str(dw1w2(j)),'n\']);
        path = replace(path, 'Delay0n','Delay0');
        path = replace(path, '000n','u');
        data=DA.readspec(path); %read trace set
        ts=data(1).data; %ts contains the traces
        ts=ts-mean(ts(:,bg1:bg2),2)*ones(1,length(data(1).t)); %background subtraction
        traces(((j-1)*nt+1):(j*(nt)),:)=ts(1:nt,x1:x2); %add traceset into traces-matrix
        av_sig(j)=mean(sum(ts(:,x1:x2),2)); %average signal for this w1w2 delay (needed for calibration)
        std_err_sig(j)=std(sum(ts(:,x1:x2),2))/sqrt(nt);
    end
    [sortedtotal,order]=sort(sum(traces,2),1,'descend');
    sortedtraces=traces(order,:);
    

    figure('Position',[0 0 900 500]) 
    subplot(1,2,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%       Density Calibration             %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    errorbar(av_sig, den, 0*std_err_sig, 0*std_err_sig,std_err_sig,std_err_sig,'.') %plot total signal over rel. density
    xlim([min(av_sig)*1.05 max(av_sig)*1.05])
    ylim([min(den)*1.05 max(den)*1.05])
    f = polyfit(av_sig',den',2); %fit with polynomial of specified order
    x = linspace(min(av_sig)-300,max(av_sig)+300);
    yfit = polyval(f,x);
    hold on;
    plot(x,yfit,'r-.'); %plot the fit in same graph
    xlabel('count','Interpreter','latex')
    ylabel('rel. density','Interpreter','latex')
    title('Density Calibration','Interpreter','latex');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%              Binning                  %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bin=zeros(n_bins+1,x2-x1+1);
    
    plasma=cell(n_bins+1,1);
    rydberg=cell(n_bins+1,1);
    
    bincount=zeros(1,n_bins+1); %number of traces in each bin
    densities=zeros(1,n_bins+1); %density of the bin

    d=polyval(f,sortedtotal); %initial densities as calculated by calibration function
    w_bins=max_den/(n_bins); %width of bins
    max_d=0;
    min_d=0;
    
    for k=1:nD*nt
        v=sortedtotal(k);
        if( d(k)>= 0 && d(k) <= max_den)%if(v>=min(av_sig) && v<= max(av_sig) && d(k)>= 0 && d(k) <= max_den) % must be positive signal
            n=round(d(k)/w_bins)+1;
            bin(n,:)=bin(n,:)+sortedtraces(k,:);
            bincount(n)=bincount(n)+1;
            densities(n)=densities(n)+d(k);
            plasma{n}=[plasma{n} sum(sortedtraces(k,1:separation))];
            rydberg{n}=[rydberg{n} sum(sortedtraces(k,separation:end))];
        end
    end
    mat=bin./(bincount');
    densities=densities./(bincount);
    
    for q=1:n_bins+1
        plasma_signal(q,i)=mean(plasma{q});
        plasma_error(q,i)=std(plasma{q})/sqrt(bincount(q));
        rydberg_signal(q,i)=mean(rydberg{q});
        rydberg_error(q,i)=std(rydberg{q})/sqrt(bincount(q));
    end
    
    ev(:,:,i)=mat; %save current SFI spectrum in array

    ylimits=[min(densities),max(densities)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%            Plot SFI                   %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2)
    clim=[0 maxcolor];
    imagesc(EF,ylimits,mat,clim)
    set(gca,'YDir','normal');
    colormap jet
    xlabel('Field in V/cm','interpreter','latex');
    ylabel('rel. density','Interpreter','latex');
    title(['t=',num2str(G(i)),'ns'],'Interpreter','latex');
    
    R=109737.30; % Rydberg constant in (cm^-1)
    B=1.98; % Rotational constant for excited NO.
    BE_0_N=R/(PQN^2); % Binding Energy of N+=0
    BE_2_N=(R/(PQN^2))-(B*2*(2+1)); %Binding Energy of N+=2
    F_0_N=(BE_0_N/4.12)^2; % Field treshold for N+=0
    F_2_N=(BE_2_N/4.12)^2; % Field treshold for N+=0
    
    
    ll=get(gca,'ylim');
    for nn=10:70
        if(nn==PQN)
            sl=0.1;
        else
            sl=0.05;
        end
        BE_0_N=R/(nn^2); % Binding Energy of N+=0
        BE_2_N=(R/(nn^2))-(B*2*(2+1)); %Binding Energy of N+=2
        F_0_N=(BE_0_N/4.12)^2; % Field treshold for N+=0
        F_2_N=(BE_2_N/4.12)^2;
        hold on
        plot([F_2_N F_2_N],[0 sl],'w-');
        plot([F_0_N F_0_N],[ll(2)*(1-sl*0.7) ll(2)],'w-');
    end

%     subplot(2,2,3)
%     mat_normalized=mat./(sum(mat,2));
%     mat_normalized(isnan(mat_normalized))=0;
% %     r=round(size(mat_normalized,1)*0.5):size(mat_normalized,1);
%     imagesc(t,ylimits,mat_normalized)%,[0 max(max(mat_normalized(r,3000:end)))])
%     set(gca,'YDir','normal');
%     colormap jet
%     %xlim(xlimits);
%     xlabel('time in $\mu$s','interpreter','latex');
%     ylabel('rel. density','Interpreter','latex');
%     title(['t=',num2str(G(i)),'ns'],'Interpreter','latex');
%    
%     subplot(2,2,4)
%     mat_den=mat./(densities');
%     mat_den(isnan(mat_den))=0;
% %     r=round(size(mat_den,1)*0.5):size(mat_den,1);
%     imagesc(t,ylimits,mat_den)%,[0 max(max(mat_den(r,3000:end)))])
%     set(gca,'YDir','normal');
%     colormap jet
%     %xlim(xlimits);
%     xlabel('time in $\mu$s','interpreter','latex');
%     ylabel('rel. density','Interpreter','latex');
%     title('Normalized by Density','Interpreter','latex');
% 
%     %save as pdf
      filename=['sfi_n',num2str(PQN), '_', num2str(n_bins)];
      export_fig([filename,'_spectra'], '-pdf', '-append')
      close
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Plot plasma and Rydber signal curves %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=flip(1:(n_bins+1))
    figure('Position',[0 0 1000 500])
    subplot(1,2,1)
    errorbar(G/1000,plasma_signal(k,:),plasma_error(k,:),'.')
    xlabel('time in $\mu$s','Interpreter','latex')
    ylabel('plasma signal (arbitr.)','Interpreter','latex')
    title(['initial density $\rho_0=$',num2str(round(densities(k),2)),'$\rho_{peak}$'],'Interpreter','latex')
    xlim([-1,21])
    subplot(1,2,2)
    errorbar(G/1000,rydberg_signal(k,:),rydberg_error(k,:),'.')
    xlabel('time in $\mu$s','Interpreter','latex')
    ylabel('Rydberg signal (arbitr.)','Interpreter','latex')
    title(['initial density $\rho_0=$',num2str(round(densities(k),2)),'$\rho_{peak}$'],'Interpreter','latex')
    xlim([-1,21])
    
%     subplot(1,3,3)
%     errorbar(G/1000,rydberg_signal(k,:)+plasma_signal(k,:),rydberg_error(k,:)+plasma_error(k,:),'.')
%     xlabel('time in $\mu$s','Interpreter','latex')
%     ylabel('Total signal (arbitr.)','Interpreter','latex')
%     title(['initial density $\rho_0=$',num2str(round(densities(k),2)),'$\rho_{peak}$'],'Interpreter','latex')
%     xlim([-1,21])
     

    export_fig([filename,'_plasma_rydber_evolution'], '-pdf', '-append')
    close
end
toc


