
# Removes voice 
ffmpeg -i Python/q1.mp4  -c copy -an Python/Outputs/video_1.mp4
ffmpeg -i Python/q2.mp4  -c copy -an Python/Outputs/video_2.mp4
ffmpeg -i Python/q3.mp4  -c copy -an Python/Outputs/video_3.mp4

# Extract frames
ffmpeg -i Python/Outputs/video_1.mp4 -vf fps=1 Python/Frames/video_1/video_1_%02d.jpg -hide_banner
ffmpeg -i Python/Outputs/video_2.mp4 -vf fps=1 Python/Frames/video_2/video_2_%02d.jpg -hide_banner
ffmpeg -i Python/Outputs/video_3.mp4 -vf fps=1 Python/Frames/video_3/video_3_%02d.jpg -hide_banner

# Resize for VGG16
for f in "Python/Frames/video_*/*.jpg"
do
     mogrify $f -resize 224x224! $f
done

