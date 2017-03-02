% Copyright (c) 2016.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Pakorn Kanchanawong
% Contact: biekp@nus.edu.sg
function AllFragments = IterGenFragment(OriginImg, ...
    AllFragments, ...
    iteration, ...
    R, ...
    ROI_Mask, ...
    NofOrientations_FT, ...
    Iter_RemoveR, ...
    Thresh, ...
    R_Junc, ...
    MIN_FragmentLength);


MasterMask = ROI_Mask;

NewImg = zeros(size(AllFragments));
NewImg((R+1):(size(OriginImg,1)+R),(R+1):(size(OriginImg,2)+R)) = OriginImg;

[temp1 num0] = bwlabel(AllFragments);
AllNumbers = num0;
for ii = 1:iteration
    display(['Information for Iteration ',num2str(ii),'.  Extraction in Progress ...']);
    [temp1 num0] = bwlabel(AllFragments);
    display(['     Previous Number of Fragment =  ',num2str(num0)]);
    
    OldAllFragments = AllFragments;
    NewImg = zeros(size(AllFragments));
    NewImg((R+1):(size(OriginImg,1)+R),(R+1):(size(OriginImg,2)+R)) = OriginImg;
    NewImg = mat2gray(NewImg);
    [x y] = find(AllFragments==1);
    
    mask2 = zeros(size(NewImg));
    for i = (-Iter_RemoveR):Iter_RemoveR
        for j = (-Iter_RemoveR):Iter_RemoveR
            mask2(sub2ind(size(mask2),x+i,y+j)) = 1;
        end
    end
    MasterMask = ROI_Mask - mask2;
    [OFT_Img, LFT_Img, LFT_Orientations] = LFT_OFT_mex(double(NewImg),double(R),double(NofOrientations_FT),double(MasterMask));
    OFT_Img = mat2gray(OFT_Img);
    
    BW = im2bw(OFT_Img,Thresh);
    RawSke = bwmorph(BW,'thin',Inf);
    [x y] = find(RawSke==1);
    
    % --- Remove Junctions
    RawSke(1:R,:) = 0;
    RawSke((end-R):end,:) = 0;
    RawSke(:,1:R) = 0;
    RawSke(:,(end-R):end) = 0;
    [x y] = find(RawSke==1);
    Allpts = [x y];
    N1 = RawSke(sub2ind(size(RawSke),Allpts(:,1)-1,Allpts(:,2)-1));
    N2 = RawSke(sub2ind(size(RawSke),Allpts(:,1)-1,Allpts(:,2)));
    N3 = RawSke(sub2ind(size(RawSke),Allpts(:,1)-1,Allpts(:,2)+1));
    N4 = RawSke(sub2ind(size(RawSke),Allpts(:,1),Allpts(:,2)-1));
    N5 = RawSke(sub2ind(size(RawSke),Allpts(:,1),Allpts(:,2)+1));
    N6 = RawSke(sub2ind(size(RawSke),Allpts(:,1)+1,Allpts(:,2)-1));
    N7 = RawSke(sub2ind(size(RawSke),Allpts(:,1)+1,Allpts(:,2)));
    N8 = RawSke(sub2ind(size(RawSke),Allpts(:,1)+1,Allpts(:,2)+1));
    N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;
    CrPts = Allpts(find(N>=3),1:2);
    for i = 1:size(CrPts,1)
        RawSke(CrPts(i,1)-R_Junc:CrPts(i,1)+R_Junc,CrPts(i,2)-R_Junc:CrPts(i,2)+R_Junc) = 0;
    end
    [x y] = find(RawSke==1);
    Allpts = [x y];
    N1 = RawSke(sub2ind(size(RawSke),x-1,y-1));
    N2 = RawSke(sub2ind(size(RawSke),x-1,y));
    N3 = RawSke(sub2ind(size(RawSke),x-1,y+1));
    N4 = RawSke(sub2ind(size(RawSke),x,y-1));
    N5 = RawSke(sub2ind(size(RawSke),x,y+1));
    N6 = RawSke(sub2ind(size(RawSke),x+1,y-1));
    N7 = RawSke(sub2ind(size(RawSke),x+1,y));
    N8 = RawSke(sub2ind(size(RawSke),x+1,y+1));
    N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;
    SinglePts = Allpts(find(N==0),1:2);
    for i = 1:size(SinglePts,1)
        RawSke(SinglePts(i,1),SinglePts(i,2)) = 0;
    end
    [x y] = find(RawSke==1);
    Allpts = [x y];
    
    % --- Remove Short Fragment
    RawSke = bwareaopen(RawSke, MIN_FragmentLength);
    [RawSke num] = bwlabel(RawSke);
    h = waitbar(0,'Checking the Curvature of Newly Detected Fragments...');
    for g = 1:num
        waitbar(g/num,h);
        [x y] = find(RawSke==g);
        flag =  LargeCurFlag(x, y);
        if flag == 1
            RawSke(sub2ind(size(RawSke),x,y)) = 0;
        end
    end
    close(h);
    RawSke(find(RawSke~=0))=1;
    [RawSke num] = bwlabel(RawSke);
    [x y] = find(RawSke~=0);
    
    AllFragments(sub2ind(size(AllFragments),x,y)) = 1;
    
    % --- Remove Junctions
    AllFragments(1:R,:) = 0;
    AllFragments((end-R):end,:) = 0;
    AllFragments(:,1:R) = 0;
    AllFragments(:,(end-R):end) = 0;
    [x y] = find(AllFragments==1);
    Allpts = [x y];
    N1 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)-1));
    N2 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)));
    N3 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)+1));
    N4 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)-1));
    N5 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)+1));
    N6 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)-1));
    N7 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)));
    N8 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)+1));
    N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;
    CrPts = Allpts(find(N>=3),1:2);
    for i = 1:size(CrPts,1)
        AllFragments(CrPts(i,1)-R_Junc:CrPts(i,1)+R_Junc,CrPts(i,2)-R_Junc:CrPts(i,2)+R_Junc) = 0;
    end
    [x y] = find(AllFragments==1);
    Allpts = [x y];
    N1 = AllFragments(sub2ind(size(AllFragments),x-1,y-1));
    N2 = AllFragments(sub2ind(size(AllFragments),x-1,y));
    N3 = AllFragments(sub2ind(size(AllFragments),x-1,y+1));
    N4 = AllFragments(sub2ind(size(AllFragments),x,y-1));
    N5 = AllFragments(sub2ind(size(AllFragments),x,y+1));
    N6 = AllFragments(sub2ind(size(AllFragments),x+1,y-1));
    N7 = AllFragments(sub2ind(size(AllFragments),x+1,y));
    N8 = AllFragments(sub2ind(size(AllFragments),x+1,y+1));
    N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;
    SinglePts = Allpts(find(N==0),1:2);
    for i = 1:size(SinglePts,1)
        AllFragments(SinglePts(i,1),SinglePts(i,2)) = 0;
    end
    
    % --- Remove Short Fragment
    AllFragments = bwareaopen(AllFragments, MIN_FragmentLength);
    [L num] = bwlabel(AllFragments,8);
    display(['     Current Number of Fragment  =  ',num2str(num)]);
end







