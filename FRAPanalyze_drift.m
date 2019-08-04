function out = FRAPanalyze_drift(fname,bleachFrame,filterWidth,channelNum,outfname)
% FRAPanalyze
% 
% Automatically detects the region in a FRAP movie that was photobleached,
% and sums the total intensity inside and outside of this region for each
% frame in the movie.
% 
% usage:
% out = FRAPanalyze_drift(fname,bleachFrame,filterWidth,channelNum,outfname)
% 
% fname - name of input file
% bleachFrame - frame number after which bleaching occurred
% filterWidth - width of Gaussian filter to use for smoothing the
% difference image
% channelNum - index of channel to use for bleaching detection/ROI
% definition. Numbering in BioFormats starts with 0.
% outfname - base name of output file
% 
% Requires BioFormats package for MATLAB
% (https://www.openmicroscopy.org/bio-formats/downloads/)
% Add BioFormats to the MATLAB path using the MATLAB addpath command prior
% to running this script
% 
% Thomas Graham, Tjian-Darzacq lab, 20190727
% 
% CHANGELOG:
% 20190803 - new version FRAPanalyze_drift corrects for drift using the
% function getMyDrift_nonlin

inparams.fname = fname;
inparams.bleachFrame = bleachFrame;
inparams.filterWidth = filterWidth;
inparams.channelNum = channelNum;
inparams.outfname = outfname;

out.inparams = inparams;


r = bfGetReader(fname);

[offx, offy] = getMyDrift_nonlin(r,bleachFrame+1,'channel',0:r.getSizeC-1);

imBefore = double(bfGetPlane(r,r.getIndex(0,channelNum,0)+1));
if bleachFrame > 1
    for j=1:bleachFrame-1
        imBefore = imBefore + double(bfGetPlane(r,r.getIndex(0,channelNum,j)+1));
    end
end

imAfter = double(bfGetPlane(r,r.getIndex(0,channelNum,bleachFrame)+1));
for j=bleachFrame:bleachFrame+bleachFrame-1
    imAfter = imAfter + double(bfGetPlane(r,r.getIndex(0,channelNum,j)+1));
end

imBefore = imBefore./mean(imBefore(:));
imAfter = imAfter./mean(imAfter(:));
imDiff = imBefore-imAfter;

gaussKern = gKern(size(imDiff),filterWidth);
imFilt = fftshift(ifft2(fft2(imDiff).*fft2(gaussKern)));

mask = imFilt > max(imFilt(:))/2;
invmask = (1-mask);

out.imBefore = imBefore;
out.imAfter = imAfter;
out.imDiff = imDiff;
out.imFilt = imFilt;
out.mask = mask;
out.offx = offx;
out.offy = offy;

figure;
subplot(1,3,1)
imagesc(imDiff); axis off; colormap gray;
subplot(1,3,2)
imagesc(imFilt); axis off;
subplot(1,3,3)
imagesc(mask); axis off;


nframes = r.getSizeT;
nchan = r.getSizeC;

fin = zeros(nframes,nchan);
fout = zeros(nframes,nchan);


w = waitbar(0,['Analyzing frame 1 of ', num2str(nframes)]);

figure;
for j=1:nframes
    currMask = circshift(mask,[-offx(max(1,j-bleachFrame)),-offy(max(1,j-bleachFrame))]);
    currInvMask = 1-currMask;
    waitbar(j/nframes,w,['Analyzing frame ' num2str(j) ' of ', num2str(nframes)]);
	for k=1:nchan
        im = double(bfGetPlane(r,r.getIndex(0,k-1,j-1)+1));
        fin(j,k) = sum(sum(im.*currMask));
        fout(j,k) = sum(sum(im.*(currInvMask)));
    end
    
	imshow(cat(3,imnorm(im),mask,currMask))
    drawnow
end

fnorm = fin./(fin+fout);
fnorm = fnorm./mean(fnorm(1:bleachFrame,:));

out.fin = fin;
out.fout = fout;
out.fnorm = fnorm;

figure;
for j=1:nchan
    plot(fnorm(:,j)); hold on;
end
legend(num2str(transpose(1:nchan)))

close(w)

save(outfname,'out')

end


function k = gKern(imageSize, width)
% k = gKern(imageSize, width)
% return normalized Gaussian kernel of specified width for 
% 
% inputs: 
% imageSize - size of 2-dimensional image
% width = width of kernel
% 
% output:
% k - gaussian kernel the same size as the image, centered in the middle
nrows = imageSize(1);
ncols = imageSize(2); 
[R,C] = meshgrid(linspace(-nrows/2,nrows/2,nrows),linspace(-ncols/2,ncols/2,ncols));
rsq = R.^2 + C.^2;
k = exp(-rsq/(2*width^2));
k = k/sum(k(:));

end

function [offx, offy] = getMyDrift_nonlin(reader,framestart,varargin) %frameend,channel,smoothingKernelWidth)%,fname)
% [offx, offy] = getMyDrift_nonlin(refim,reader,framestart,frameend,period,border,fname)
%
% Use image correlation to calculate the optimal linear translation in x and y.
%
% getMyDrift_nonlin(reader,framestart,frameend,channel,smoothingKernelWidth)
% 
% reader - BioFormats reader object (output of bfGetReader function)
% 
% framestart - first frame to start with
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Optional parameters:
% 
% frameend - frame to end with; if not specified, all frames will be
% analyzed
% 
% channel - channel number(s) to use for drift correction. Accepts either a
% single number or an array. If multiple channels are selected, it sums the
% images in all of those channels. default = 0
% 
% smoothingKernelWidth - width in pixels of Gaussian kernel to use for
% image smoothing. default = 10
% 
% movingAverageWidth - width of moving average box (to prevent jumping back
% and forth by 1 px). default = 10

ip = inputParser;
ip.addRequired('reader');
ip.addRequired('framestart');
ip.addParameter('frameend',inf,@isnumeric);
ip.addParameter('channel',0,@isnumeric);
ip.addParameter('smoothingKernelWidth',10,@isnumeric);
ip.addParameter('movingAverageWidth',20,@isnumeric);
ip.parse(reader,framestart,varargin{:});
frameend = ip.Results.frameend;
channel = ip.Results.channel;
smoothingKernelWidth = ip.Results.smoothingKernelWidth;
maw = ip.Results.movingAverageWidth;

if isinf(frameend)
    frameend = reader.getSizeT;
end

nframes = frameend - framestart + 1;
w = waitbar(0,['Determining offsets: Frame 0 of ' num2str(nframes)]);

if numel(channel) > 1
    refim = double(bfGetPlane(reader,reader.getIndex(0,channel(1),framestart-1)+1));
    for j=2:numel(channel)
        refim = refim + double(bfGetPlane(reader,reader.getIndex(0,channel(j),framestart-1)+1));
    end
else
    refim = double(bfGetPlane(reader,reader.getIndex(0,channel,framestart-1)+1));
end
refim = imnorm(refim);

% figure(3); imagesc(refim); colormap gray;

% Gaussian filtering kernel to use for imcorr2
imsize = size(refim);
[X,Y] = meshgrid(linspace(-(imsize(1)-1)/2,(imsize(1)-1)/2,imsize(1)),...
    linspace(-(imsize(2)-1)/2,(imsize(2)-1)/2,imsize(2)));
R2 = X.^2+Y.^2;
gaussfn = exp(-R2/smoothingKernelWidth/2);
gaussfn = gaussfn/sum(gaussfn(:));

offxsample = zeros(1,frameend-framestart+1);
offysample = zeros(1,frameend-framestart+1);

% Calculate offset for the channel(s) of interest
for j=framestart:frameend
    waitbar((j-framestart)/nframes, w, ['Determining offsets: Frame ' num2str(j-framestart) ...
           ' of ' num2str(nframes)])
    fcurr = double(bfGetPlane(reader,reader.getIndex(0,channel(1),j-1)+1)); 
    if numel(channel) > 1
        for k=1:numel(channel)
            fcurr = fcurr + double(bfGetPlane(reader,reader.getIndex(0,channel(k),j-1)+1)); 
        end
    end
    fcurr = imnorm(fcurr);

%     clf
%     imagesc(fcurr)
%     drawnow

    currcorr = imcorr2(fcurr,refim,gaussfn);
    offxsample(j-framestart+1) = currcorr(1);
    offysample(j-framestart+1) = currcorr(2); 
    
%     hold on;
%     plot(currcorr(1),currcorr(2))
    
end
close(w);


offx = round(movmean(offxsample,maw));
offy = round(movmean(offysample,maw));

figure; plot(offxsample-offx(1),'bo');
hold on; 
plot(offysample-offy(1),'ro');

offx = offx - offx(1);
offy = offy - offy(1);

plot(offx,'b-','LineWidth',2);
plot(offy,'r-','LineWidth',2);

end

function normed = imnorm(im)
    normed = (im-min(im(:)))/(max(im(:))-min(im(:)));
end

function bestcorr = imcorr2(im1, im2, filteringKernel)
% usage: bestcorr = imcorr2(im1, im2, filteringKernel)
% 
% bestcorr - translation in x and y that gives the best pixel-wise
% correlation between im1 and im2
% 
% im1 and im2 - input images, which must be the same size
% 
% filteringKernel - optional kernel (e.g., Gaussian function) to use for
% filtering the two images; must be an array the same size as im1 and im2

sim1 = size(im1);
sim2 = size(im2);
assert(~sum(sim1~=sim2),...
 ['Error in function bestcorr: The size of the two images must be equal.',...
 '\nimage 1 size = [%d,%d]\nimage 2 size = [%d,%d]'],sim1(1),sim1(2),sim2(1),sim2(2));

bestcorr = [0 0];

corr = [];

im1 = double(im1);
im1 = im1-mean(im1(:));
im2 = double(im2);
im2 = im2-mean(im2(:));

if exist('filteringKernel','var')
    %filteringKernel = fftshift(filteringKernel);
    corr = ifft2(conj(fft2(im1)).*fft2(im2));
    corr = real(ifft2(conj(fft2(filteringKernel)).*fft2(corr)));
else
    corr = real(fftshift(ifft2(conj(fft2(im1)).*fft2(im2))));
end

imagesc(corr);

[maxes, besti] = max(corr);
[~, bestj] = max(maxes);

bestcorr = [besti(bestj), bestj];

end