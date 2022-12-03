function [rgb,numframes]  =    yuv2rgb(pj_dir,data_dir,test_file)
yuv_path          =      fullfile(pj_dir,data_dir,strcat(test_file,'.yuv'));
fid               =      fopen(yuv_path,'rb');
while (fid == -1)
      fid         =      fopen(yuv_path,'rb');
end
if test_file(end-3) == 'q'
    row           =      176 ;
    col           =      144;
else
    row           =      176*2;
    col           =      144*2;
end
numframes         =      50; % total=300
Y                 =      zeros(row,col,numframes);
U                 =      zeros(row/2,col/2,numframes);
V                 =      zeros(row/2,col/2,numframes);
UU                =      zeros(row,col,numframes);%
VV                =      zeros(row,col,numframes);
yuv               =      zeros(col,row,3,numframes);


for frame=1:numframes
    [Y(:,:,frame),count]              =     fread(fid,[row,col],'uchar');
    [U(:,:,frame),count1]             =     fread(fid,[row/2,col/2],'uchar');
    [V(:,:,frame),count2]             =     fread(fid,[row/2,col/2],'uchar');
    
    UU(1:2:row-1,1:2:col-1,frame)     =     U(:,:,frame);%
    UU(1:2:row-1,2:2:col,frame)       =     U(:,:,frame);%
    UU(2:2:row,1:2:col-1,frame)       =     U(:,:,frame);%
    UU(2:2:row,2:2:col,frame)         =     U(:,:,frame);%
    
    VV(1:2:row-1,1:2:col-1,frame)     =     V(:,:,frame);
    VV(1:2:row-1,2:2:col,frame)       =     V(:,:,frame);
    VV(2:2:row,1:2:col-1,frame)       =     V(:,:,frame);
    VV(2:2:row,2:2:col,frame)         =     V(:,:,frame);
    
%     yuv(:,:,1,frame)                  =     Y(:,:,frame)';
%     yuv(:,:,2,frame)                  =     UU(:,:,frame)';
%     yuv(:,:,3,frame)                  =     VV(:,:,frame)';
    rgb(:,:,1,frame) = (Y(:,:,frame) + 1.140 * (VV(:,:,frame)-128 ))';
    rgb(:,:,2,frame) = (Y(:,:,frame) - 0.395 * (UU(:,:,frame)-128 ) - 0.581 *(VV(:,:,frame)-128))';
    rgb(:,:,3,frame) = (Y(:,:,frame) + 2.032 *(UU(:,:,frame)-128))';
    
   
    
end
rgb(rgb<0)                           =     0;
rgb(rgb>255)                         =     255;
fclose(fid);