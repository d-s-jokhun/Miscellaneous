
%%% written by D.S.JOKHUN on 05/12/2016
%%% Change line 10 according to the file/s you want to study
%%% results are saved in the variable 'result'

clear all
clc


sample = 'transfer';
filenames = dir ([sample,'*1800cir*','.nd2']);
rslt=[];
result_Pos=[];
result_pro_nuc_A=[];
result_AR=[];
result_Eccentricity=[];
result_Orient=[];
result_Corr=[];

cell_num=0;
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
%     MetaData.Voxel_Size_Z = double(OmeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM)); % in µm
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
    for iSeries =1:MetaData.SeriesCount;  %%%choosing a specific XY 3D point from the multipoint image
        iSeries 
        cell_num=cell_num+1;
        XYZ=[];
        Reader.setSeries(iSeries - 1);
        for iT=1:MetaData.TimePoints;
            iT
            for iCh=1%1:MetaData.Num_of_Ch;
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

       % Finding the nucleus
       
       Ch1_pro=uint16(sum(XYZ{1,1},3));
       Ch1_filt_pro=uint16(sum(medfilt3(XYZ{1,1}),3));

       
%           %% peak method%%
%           Ch1_pro_dist=imhist(Ch1_pro);
%           Ch1_pro_dist(1,1)=0;
%           [peak_vals,peak_ind]=findpeaks(Ch1_pro_dist);
%           peak = horzcat (peak_vals,peak_ind);
%           for p=1:size(peak,1)
%               if peak(p,1)<0.9*max(peak(:,1))
%                  peak(p,1)=0;
%                  peak(p,2)=0;
%               end
%           end
%           std_dv_Ch1_pro=std(nonzeros(double(reshape(Ch1_pro,[],1))));
%           
%           Ch1_pro_ThreshLevel = (mean(nonzeros(peak(:,2)))-(1*std_dv_Ch1_pro*(255/65535)))/(size(Ch1_pro_dist,1)-1);  %%threshold is mode + 3 x standard deviation
%           %%%%
          
          %% mean method%%
          mean_Ch1_filt_pro=mean(nonzeros(double(reshape(Ch1_filt_pro,[],1))));
          std_dv_Ch1_filt_pro=std(nonzeros(double(reshape(Ch1_filt_pro,[],1))));
          Ch1_filt_pro_ThreshLevel = (mean_Ch1_filt_pro+(1.5*std_dv_Ch1_filt_pro))/65535;
          %%%%%
       
          
          
          if cell_num==1
              thresh_multiplier=0.8
          end
          if cell_num==2
              thresh_multiplier=0.8
          end
          if cell_num==3
              thresh_multiplier=0.8
          end
          
          
          Ch1_filt_pro_bw=[];
          Ch1_filt_pro_bw=im2bw(Ch1_filt_pro,thresh_multiplier*Ch1_filt_pro_ThreshLevel);
          
          Ch1_filt_pro_bw = imfill (Ch1_filt_pro_bw,'holes');
          
          eroded_Ch1_filt_pro_bw=imerode(Ch1_filt_pro_bw,strel('disk',1));
       
          eroded_Ch1_filt_pro_bw_CC = bwconncomp(eroded_Ch1_filt_pro_bw,6);
          real_nuc_pro_bw=[];
          real_nuc_pro_bw(1:MetaData.Num_of_Pixels_Y,1:MetaData.Num_of_Pixels_X) = 0;
          for count=1:size(eroded_Ch1_filt_pro_bw_CC.PixelIdxList,2);
           if size(eroded_Ch1_filt_pro_bw_CC.PixelIdxList{count},1)>5000    %5000 pix is the minimum size required to consider an object as the nucleus
               real_nuc_pro_bw(eroded_Ch1_filt_pro_bw_CC.PixelIdxList{count})=1;
           end
          end
                    
          
          nuc_stats = regionprops(real_nuc_pro_bw,'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation', 'BoundingBox');
          rslt(size(rslt,1)+1,1)=file_count;
          rslt(size(rslt,1),2)=iSeries;
          rslt(size(rslt,1),3)=iT;
          rslt(size(rslt,1),4:5)=(nuc_stats.Centroid) .* MetaData.Voxel_Size_X;
            result_Pos (iT,(cell_num*2)-1:(cell_num*2))=(nuc_stats.Centroid) .* MetaData.Voxel_Size_X; % XY XY XY
          rslt(size(rslt,1),6)=nuc_stats.Area*MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y;
            result_pro_nuc_A (iT,(cell_num*2)-1)=iT;
            result_pro_nuc_A (iT,(cell_num*2))=nuc_stats.Area*MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y; % tA tA tA tA
          rslt(size(rslt,1),7)=nuc_stats.MajorAxisLength/nuc_stats.MinorAxisLength;
            result_AR (iT,(cell_num*2)-1)=iT;
            result_AR (iT,(cell_num*2))=nuc_stats.MajorAxisLength/nuc_stats.MinorAxisLength; % tAR tAR tAR tAR
          rslt(size(rslt,1),8)=nuc_stats.Eccentricity;
            result_Eccentricity (iT,(cell_num*2)-1)=iT;
            result_Eccentricity (iT,(cell_num*2))=nuc_stats.Eccentricity; % tE tE tE tE
          rslt(size(rslt,1),9)=nuc_stats.Orientation;
            result_Orient (iT,(cell_num*2)-1)=iT;
            result_Orient (iT,(cell_num*2))=nuc_stats.Orientation; % tO tO tO tO
            
          
         if iT==1
             ref_orient = nuc_stats.Orientation;
         end
          rotated_translated_Ch1_pro = imtranslate (Ch1_pro,[-(nuc_stats.Centroid(1,1)-(MetaData.Num_of_Pixels_X/2)) -(nuc_stats.Centroid(1,2)-(MetaData.Num_of_Pixels_Y/2))]);
          rotated_translated_Ch1_pro = imrotate (rotated_translated_Ch1_pro, -(nuc_stats.Orientation-ref_orient), 'crop');
          

          
         if iT==1
             factor=1.5;
             cropping_mask=[(MetaData.Num_of_Pixels_X/2)-((factor*nuc_stats.BoundingBox(3))/2) (MetaData.Num_of_Pixels_Y/2)-((factor*nuc_stats.BoundingBox(4))/2) factor*nuc_stats.BoundingBox(3) factor*nuc_stats.BoundingBox(4)];
             ref_img=imcrop(rotated_translated_Ch1_pro,cropping_mask);
         end
         current_img=imcrop(rotated_translated_Ch1_pro,cropping_mask);
         
         pre_corr_mask= ((ref_img>0)+(current_img>0))==2;
         pre_corr_ref_img= uint16(pre_corr_mask).*ref_img;
         pre_corr_current_img= uint16(pre_corr_mask).*current_img;
         
         img_corr2 = corr2(pre_corr_current_img,pre_corr_ref_img);
          
          
          rslt(size(rslt,1),10)=img_corr2;
            result_Corr (iT,(cell_num*2)-1)=iT;
            result_Corr (iT,(cell_num*2))=img_corr2; % tC tC tC tC
          
          
          
          
       %% PERFORM ANALYSIS ABOVE!!!            
        end
    end
     
end

% plot(rslt(:,3),rslt(:,10))
rslt_headers = {'File' 'MultiPoint' 'TimePoint' 'X pos /um' 'Y pos /um' 'Nuc pro Area /um2' 'Aspect Ratio' 'Eccentricity' 'Orientation /deg' 'correlation with img 0' };
result=vertcat(rslt_headers,num2cell(rslt));