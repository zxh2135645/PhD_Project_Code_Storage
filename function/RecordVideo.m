function RecordVideo(movieVector, fname, FrameRate)

 myWriter = VideoWriter(fname);
 myWriter.FrameRate = FrameRate;
 
 % Open the videowriter object, write the movie and close the file
 open(myWriter);
 writeVideo(myWriter, movieVector);
 close(myWriter);
 
end