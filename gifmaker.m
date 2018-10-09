function []=gifmaker(H,i,filename)
IMAGE=frame2im(H);
[imind,cm]=rgb2ind(IMAGE,256);
if i==1;
imwrite(imind,cm,filename,'gif','Loopcount',inf);
else
imwrite(imind,cm,filename,'gif','WriteMode','append');   
end
end