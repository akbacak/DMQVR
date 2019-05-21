


mov=VideoReader('myDataset/Videofolder/1.mp4');
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  %imshow(imread(fname)); axis image;
  imshow(videoFrame);

end
