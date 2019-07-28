function out = FRAPanalyze(fname,bleachFrame,filterWidth,channelNum,outfname)
% FRAPanalyze
% 
% Automatically detects the region in a FRAP movie that was photobleached,
% and sums the total intensity inside and outside of this region for each
% frame in the movie.
% 
% usage:
% out = FRAPanalyze(fname,bleachFrame,filterWidth)
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

% TODO
% Add drift correction

inparams.fname = fname;
inparams.bleachFrame = bleachFrame;
inparams.filterWidth = filterWidth;
inparams.channelNum = channelNum;
inparams.outfname = outfname;

out.inparams = inparams;


r = bfGetReader(fname);

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

for j=1:nframes
    waitbar(j/nframes,w,['Analyzing frame ' num2str(j) ' of ', num2str(nframes)]);
	for k=1:nchan
        im = double(bfGetPlane(r,r.getIndex(0,k-1,j-1)+1));
        fin(j,k) = sum(sum(im.*mask));
        fout(j,k) = sum(sum(im.*(invmask)));
    end
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

