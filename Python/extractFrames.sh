
# Removes voice 
ffmpeg -i Python/q1.mp4  -c copy -an Python/video_1.mp4
ffmpeg -i Python/q2.mp4  -c copy -an Python/video_2.mp4

# Extract frames
ffmpeg -i Python/video_1.mp4 -vf fps=1 Python/Frames/Frames_1/video_1_%02d.jpg -hide_banner
ffmpeg -i Python/video_2.mp4 -vf fps=1 Python/Frames/Frames_2/video_2_%02d.jpg -hide_banner


# Resize for VGG16
for f in "Python/Frames/Frames_*/*.jpg"
do
     mogrify $f -resize 224x224! $f
done

