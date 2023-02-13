#https://stackoverflow.com/questions/25569180/ffmpeg-convert-without-loss-quality
ffmpeg -framerate 10 -i frame%04dfig0.png -qscale 0 ssha_Model_B.mp4
echo "Done."
