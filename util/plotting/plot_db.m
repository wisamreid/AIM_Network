function plot_db(input,thresh)%plot for debugging
temp = db(input);
temp(temp<-thresh) = -thresh;
imagesc(temp);
end