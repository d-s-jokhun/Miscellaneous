%%% Don't forget to select the appropriate channel
%%% written by D.S.JOKHUN on 23/11/2016


clear all


sample = 'L';
filenames = dir ([sample,'*.tif']);
actin_channel=1;

cell_num=0;
acute_angles_frm_mean=[];
percent_explained=[];
preS=[];
acute_angles_frm_mean_combined=[];
percent_explained_combined=[];
S_combined=[];

for file_count=1:size(filenames,1);
    
    filename = filenames(file_count).name

    Reader = bfGetReader (filename);
    OmeMeta = Reader.getMetadataStore();

    MetaData.SeriesCount = Reader.getSeriesCount();
    MetaData.TimePoints = OmeMeta.getPixelsSizeT(0).getValue();
    MetaData.Num_of_Ch = OmeMeta.getPixelsSizeC(0).getValue();
    MetaData.Num_of_Pixels_Z = OmeMeta.getPixelsSizeZ(0).getValue();
    MetaData.Num_of_Pixels_X = OmeMeta.getPixelsSizeX(0).getValue();
    MetaData.Num_of_Pixels_Y = OmeMeta.getPixelsSizeY(0).getValue();
    MetaData.Voxel_Size_X = double(OmeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in µm
    MetaData.Voxel_Size_Y = double(OmeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in µm
    MetaData.Voxel_Size_Z = double(OmeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM)); % in µm
%     MetaData.Plane_Origin_X = double(OmeMeta.getPlanePositionX(0,0).value);
%     MetaData.Plane_Origin_Y = double(OmeMeta.getPlanePositionY(0,0).value);
%     MetaData.Plane_Origin_Z = double(OmeMeta.getPlanePositionZ(0,0).value);

    MetaData.ChannelID = [];
    for ch_count = 1:MetaData.Num_of_Ch ;
        chID_temp = ['   ' char(num2str(ch_count-1)) '.'] ;
        chNAME_temp= [char(OmeMeta.getChannelName(0,ch_count-1))];
        MetaData.ChannelID = [MetaData.ChannelID  chID_temp chNAME_temp];  % (series 0, channel ch_count)
    end
    MetaData


    
    %%% loading particular 3D frames
    for iSeries = 1:MetaData.SeriesCount;  %%%choosing a specific XY 3D point from the multipoint image
        iSeries 
        XYZ=[];
        Reader.setSeries(iSeries - 1);
        for iT=1:MetaData.TimePoints;
            iT
            for iCh=actin_channel%1:MetaData.Num_of_Ch;
                XYZ_temp =uint16([]);
                for iZ=1:MetaData.Num_of_Pixels_Z;
                    iPlane = Reader.getIndex(iZ-1, iCh-1, iT-1) + 1;     %%% The last '1-1' is for timepoint 0 (the 1st timepoint)
                	XYZ_temp(:,:,iZ)= bfGetPlane(Reader, iPlane);
                end
                XYZ{iCh}=XYZ_temp;   %%% 1st element of XYZ will be a 3D matrix of series 1 (XY1) in Channel 1
                                 %%% 2st element of XYZ will be a 3D matrix of series 1 (XY1) in Channel 2                    
            end 
        
        
       %%% {XYZ} currently has the 3D intensity matrices of all the channels at iT(timepoint iT) from iSeries(multipoint i) found in file f 
       %% PERFORM ANALYSIS BELOW!!!
       cell_num=cell_num+1;
       angle_distri_cell_n=[];
       percent_explained_cell_n=[];
       
       ch=uint16(sum(XYZ{1,iCh},3)); 
       
       bckgrnd_multiple =1;
       max_back= [max(max(ch(4:18,4:18))) , max(max(ch(4:18,end-17:end-3))) , max(max(ch(end-17:end-3,4:18))) , max(max(ch(end-17:end-3,end-17:end-3)))];   % finding the intensities at the corners.
       max_back = double(max(max_back));
       cell_thresh_lvl= (bckgrnd_multiple*max_back)/65535;  %adjust the multiple to adjust the threshold level
       cell_bw = im2bw(ch,cell_thresh_lvl);
       cell_bw = imfill(cell_bw,'holes');
       erode_radius=2 %in micron
       erode_radius = round(erode_radius/MetaData.Voxel_Size_X) ; %in pixels
       cell_bw = imerode (cell_bw,strel('disk',erode_radius));
       
       
       actin=uint16(cell_bw).*ch;
       imtool(actin,[])
       
       box_size=2; % in micron
       box_size=round(box_size/MetaData.Voxel_Size_X); % in pixel
       if box_size<6
           box_size=6;
       end
       ['box size=',num2str(box_size),' pixels']
       ['box size=',num2str(box_size*MetaData.Voxel_Size_X),' microns']
       num_box_x=size(actin,2)-(box_size-1);
       num_box_y=size(actin,1)-(box_size-1);
       
       for box_num_x = 1:box_size:num_box_x  %1:num_box_x 
           cell_num
           ['box num in x=',num2str(box_num_x),' of ',num2str(num_box_x)]
           for box_num_y = 1:box_size:num_box_y  %1:num_box_y 
               box_n=actin(box_num_y:box_num_y+(box_size-1),box_num_x:box_num_x+(box_size-1));
               if sum(sum(box_n>0))>=(box_size^2)/4%==box_size^2 %making sure that only full boxes are considered
                    freq_distri_box_n=[];
                    cumu_freq_box_n=0;
                    for county=1:size(box_n,1)
                        for countx=1:size(box_n,2)
                            if box_n(county,countx)>0
                                freq_distri_box_n(cumu_freq_box_n+1:cumu_freq_box_n+box_n(county,countx),1)=countx;
                                freq_distri_box_n(cumu_freq_box_n+1:cumu_freq_box_n+box_n(county,countx),2)=size(box_n,1)+1-county;
                                cumu_freq_box_n=cumu_freq_box_n+box_n(county,countx);
                            end
                        end
                    end
                    [coeff,~,~,~,explained] = pca(freq_distri_box_n);
                    percent_variance_box_n=explained(1,1);
                    angle_box_n=round(atand (coeff(2,1)/coeff(1,1)));
       
                    if percent_variance_box_n>50
                        if percent_variance_box_n<100
                            angle_distri_cell_n(end+1,1)=angle_box_n;
                            percent_explained_cell_n(end+1,1)=percent_variance_box_n;
                        end
                    end

               end
           end
       end
       
       
       
       angles_frm_mean_cell_n=angle_distri_cell_n-round(mean(angle_distri_cell_n));
       %getting only the acute angles
       a=((angles_frm_mean_cell_n<-90)*180)+((angles_frm_mean_cell_n<-90).*angles_frm_mean_cell_n);
       b=((angles_frm_mean_cell_n>90).*angles_frm_mean_cell_n)-((angles_frm_mean_cell_n>90)*180);
       c=((-90<=angles_frm_mean_cell_n).*(angles_frm_mean_cell_n<=90)).*angles_frm_mean_cell_n;
       acute_angles_frm_mean_cell_n=a+b+c;
       
%        preS_cell_n=cosd(2*acute_angles_frm_mean_cell_n);
       preS_cell_n=((-1/45)*abs(acute_angles_frm_mean_cell_n))+1;
       acute_angles_frm_mean(1:size(acute_angles_frm_mean_cell_n,1),cell_num)=acute_angles_frm_mean_cell_n;
       percent_explained(1:size(percent_explained_cell_n,1),cell_num)=percent_explained_cell_n; 
       preS(1:size(preS_cell_n,1),cell_num)=preS_cell_n; 
       acute_angles_frm_mean_combined=vertcat(acute_angles_frm_mean_combined,acute_angles_frm_mean_cell_n);
       percent_explained_combined=vertcat(percent_explained_combined,percent_explained_cell_n);
       S_combined=vertcat(S_combined,mean(preS_cell_n));
        
        
        
        
       %% PERFORM ANALYSIS ABOVE!!!            
        end
    end
     
end


