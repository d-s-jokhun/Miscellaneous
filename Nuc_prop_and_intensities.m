
%%% written by D.S.JOKHUN on 06/03/2017
%%% results are saved in the variable 'result'



clear all
clc


sample = 'cir_1hr';
filenames = dir (['*',sample,'*','.nd2']);
rslt=[];
result=[];



cell_num=0;

for file_count=1:size(filenames,1);
    
    
    background_box =[];
    
    
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
    for iSeries =1:MetaData.SeriesCount;  %%%choosing a specific XY 3D point from the multipoint image
        iSeries 
        XYZ=[];
        Reader.setSeries(iSeries - 1);
        for iT=1:MetaData.TimePoints;
%             iT
            for iCh=1:MetaData.Num_of_Ch;
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

       
       % Mask for background intensities
       stripe_thickness=10; %in pixels
       background_box(1:MetaData.Num_of_Pixels_Y,1:MetaData.Num_of_Pixels_X,1:MetaData.Num_of_Pixels_Z) =1;
       background_box(stripe_thickness+1:MetaData.Num_of_Pixels_Y-stripe_thickness,stripe_thickness+1:MetaData.Num_of_Pixels_X-stripe_thickness,1:MetaData.Num_of_Pixels_Z)=0;
       num_of_bckgrnd_pix=sum(sum(sum(background_box)));
       
       
       % Finding the nucleus
       
       Ch1_background_box = XYZ{1,1}.*uint16(background_box);
       Ch1_bckgrnd_per_pix = mean(nonzeros(Ch1_background_box));
       Ch1_bckgrnd_per_pix_SD = std(double(nonzeros(Ch1_background_box)));
       Nuc_Thresh = (5*(Ch1_bckgrnd_per_pix+(5*Ch1_bckgrnd_per_pix_SD)))/((2^16)-1);   %%% thresholding the nucleus
       nuc_bw=reshape(imbinarize(reshape(XYZ{1,1},[],1),Nuc_Thresh),[MetaData.Num_of_Pixels_Y,MetaData.Num_of_Pixels_X,MetaData.Num_of_Pixels_Z]);

       
       
       
       
       nuc_bw_CC=bwconncomp(nuc_bw,6);
       real_nuc_bw=[];
       for count=1:size(nuc_bw_CC.PixelIdxList,2);
           if size(nuc_bw_CC.PixelIdxList{count},1)>10000    %5000 pix is the minimum size required to consider an object as the nucleus
               real_nuc_bw(1:MetaData.Num_of_Pixels_Y,1:MetaData.Num_of_Pixels_X,1:MetaData.Num_of_Pixels_Z) = 0;
               real_nuc_bw(nuc_bw_CC.PixelIdxList{count})=1;
               
               
               % nuc prop
               real_nuc_bw_filled = imfill (real_nuc_bw,'holes');
               nuc_2D_stats = regionprops(sum(real_nuc_bw_filled,3)>0,'Area','Eccentricity','MajorAxisLength','MinorAxisLength');
               
               if nuc_2D_stats.Area*MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y>70
                   
                   cell_num=cell_num+1
                   
                   real_nuc_edge=imdilate(edge(sum(real_nuc_bw,3)>0),strel('disk',1));
imtool(sum(XYZ{1,1},3)+(real_nuc_edge*max(max(sum(XYZ{1,1},3))))+((sum(background_box,3)>0)*max(max(sum(XYZ{1,1},3)))),[]);
                   
                   
                   rslt(size(rslt,1)+1,1)=file_count;
                   rslt(size(rslt,1),2)=iSeries;
                   rslt(size(rslt,1),3)=cell_num;
                   rslt(size(rslt,1),4)=nuc_2D_stats.Area*MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y;
                   rslt(size(rslt,1),5)=sum(sum(sum(real_nuc_bw_filled)))*MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y*MetaData.Voxel_Size_Z;
                   rslt(size(rslt,1),6)=nuc_2D_stats.MajorAxisLength/nuc_2D_stats.MinorAxisLength;
                   rslt(size(rslt,1),7)=nuc_2D_stats.Eccentricity;
                   
                   
                   
                   % Ch1 analyses
                   Ch1_background_box=XYZ{1,1}.*uint16(background_box);
                   Ch1_bckgrnd_per_pix = mean(nonzeros(Ch1_background_box));
                   
                   Ch1_in_nuc = XYZ{1,1}.*uint16(real_nuc_bw);
                   net_nuc_Ch1_int = sum(sum(sum(Ch1_in_nuc))) - (Ch1_bckgrnd_per_pix*sum(sum(sum(real_nuc_bw))));
                   net_nuc_Ch1_int_per_pix=net_nuc_Ch1_int/sum(sum(sum(real_nuc_bw)));
                   
                   rslt(size(rslt,1),8)= net_nuc_Ch1_int;
                   rslt(size(rslt,1),9)= net_nuc_Ch1_int_per_pix;
                   
                   
                   
                   % Ch2 analyses
                   Ch2_background_box=XYZ{1,2}.*uint16(background_box);
                   Ch2_bckgrnd_per_pix = mean(nonzeros(Ch2_background_box));
                   
                   Ch2_in_nuc = XYZ{1,2}.*uint16(real_nuc_bw);
                   net_nuc_Ch2_int = sum(sum(sum(Ch2_in_nuc))) - (Ch2_bckgrnd_per_pix*sum(sum(sum(real_nuc_bw))));
                   net_nuc_Ch2_int_per_pix=net_nuc_Ch2_int/sum(sum(sum(real_nuc_bw)));
                   
                   rslt(size(rslt,1),10)= net_nuc_Ch2_int;
                   rslt(size(rslt,1),11)= net_nuc_Ch2_int_per_pix;
                   
                   rslt(size(rslt,1),12)=net_nuc_Ch2_int/net_nuc_Ch1_int;
                   
                   
                   % Ch3 analyses
                   Ch3_background_box=XYZ{1,3}.*uint16(background_box);
                   Ch3_bckgrnd_per_pix = mean(nonzeros(Ch3_background_box));
                   
                   Ch3_in_nuc = XYZ{1,3}.*uint16(real_nuc_bw);
                   net_nuc_Ch3_int = sum(sum(sum(Ch3_in_nuc))) - (Ch3_bckgrnd_per_pix*sum(sum(sum(real_nuc_bw))));
                   net_nuc_Ch3_int_per_pix=net_nuc_Ch3_int/sum(sum(sum(real_nuc_bw)));
                   
                   rslt(size(rslt,1),13)= net_nuc_Ch3_int;
                   rslt(size(rslt,1),14)= net_nuc_Ch3_int_per_pix;
                   
                   rslt(size(rslt,1),15)=net_nuc_Ch3_int/net_nuc_Ch1_int;
                   
                   rslt(size(rslt,1),16)=net_nuc_Ch2_int/net_nuc_Ch3_int;   %gH2AX to total H2AX ratio
                   
                   
%                    
%                    % In absence of Ch3 (H2AX)
%                    rslt(size(rslt,1),13)= 0;
%                    rslt(size(rslt,1),14)= 0;
%                    rslt(size(rslt,1),15)= 0;
%                    rslt(size(rslt,1),16)= 0;
%                    
                   


%                    %% Finding gH2AX foci
%                    
%                    Ch2_in_nuc_pro=uint16(sum(Ch2_in_nuc,3));
%                    
%                    
%                    if mod(floor(5/MetaData.Voxel_Size_X),2)==1
%                        adaptthresh_size=floor(5/MetaData.Voxel_Size_X);
%                    else
%                        adaptthresh_size=ceil(5/MetaData.Voxel_Size_X);
%                    end
%                    Ch2_in_nuc_bw=imbinarize((Ch2_in_nuc_pro),1.5*adaptthresh(uint16(sum(XYZ{1,2},3)),0.5,'NeighborhoodSize',[adaptthresh_size adaptthresh_size]));
%                    Ch2_nodes_in_nuc=Ch2_in_nuc_pro.*uint16(Ch2_in_nuc_bw);
%                    Ch2_nodes_in_nuc_gauss=imgaussfilt3(Ch2_nodes_in_nuc,0.5);
%                    Ch2_nodes_2_gauss_med=medfilt3(Ch2_nodes_in_nuc_gauss);
%                    Ch2_nodes_2_gauss_med_eroded=imerode(Ch2_nodes_2_gauss_med,strel('disk',1));
% %                    Ch2_nodes_2_gauss_med_eroded_dilated=imdilate(Ch2_nodes_2_gauss_med_eroded,strel('disk',1));
%                    
%                    
%                    imtool(Ch2_in_nuc_pro,[])
%                    imtool(sum(Ch2_nodes_2_gauss_med_eroded,3))
% 
%                    
%                    Ch2_nodes_CC=bwconncomp(Ch2_nodes_2_gauss_med_eroded>0,4);
%                    rslt(size(rslt,1),13)=Ch2_nodes_CC.NumObjects;
                   





               end
           end
       end
       
       
          
       %% PERFORM ANALYSIS ABOVE!!!            
        end
    end
     
end


rslt_headers = {'File' 'MultiPoint' 'Cell Num' 'Nuc pro Area /um2' 'Nuc Volume /um3' 'Aspect Ratio' 'Eccentricity' 'Net DAPI in nuc' 'Net DAPI per pix in nuc' 'Net gH2AX in Nuc' 'Net gH2AX per pix in Nuc' 'Net gH2AX/Net DAPI in nuc' 'Net H2AX in Nuc' 'Net H2AX per pix in Nuc' 'Net H2AX/Net DAPI in nuc' 'Net gH2AX/Net H2AX in Nuc' };
result=vertcat(rslt_headers,num2cell(rslt));

