# Interactive-Pool-Table
A new and unique way to help in learning the game of pool.
# Purpose
To make a fun and better learning environment for pool or snooker.
# How it works
wbsfunc.m is the main file here.
The function firstly takes input an overhead image of a pool table. It is then resized and converted to double data type for faster and better calculations. After this, a mask is applied on the image to remove the background color, i.e. the color of the table . This is done to make it easier to separate out the balls and cue stick in the image. The PictureMask function applied here is an auto-generated function in MATLAB using colorThresholder, an inbuilt toolbox in MATLAB which has been done using some sample images of the table and can be generated for any kind of table similarly. After this, MATLAB function imfindcircles is applied on the masked out image to find the position of the balls on the table. Once these balls are identified, a new image of same size is created with only the circles(or balls) inserted in the image using the function insertShape already available in MATLAB. This is shown below.
![alt text](https://github.com/himanshudua/Interactive-Pool-Table/blob/master/Pool%20Photos/1.jpg)
As this is done, the median of the color in the selected circular regions is calculated in the original image and then is applied to the whole circle. This is done so that we can easily find the white ball on the table by applying threshold values and selecting only the circle(s) satisfying the threshold conditions. In case more than one balls are selected, the one closer to the cue stick is chosen as the white ball( as it is the logical choice in case the player is about to take a shot). Now, the image where we have found out the balls is subtracted from the masked out image in order to separate out the remaining region which mainly contains the cue stick. From this image, the median of all the white pixel positions is calculated and the distance of the centers of all the balls selected as white balls from this point and the ball with minimum distance is selected as the white ball. This is as shown below.
![alt text](https://github.com/himanshudua/Interactive-Pool-Table/blob/master/Pool%20Photos/2.jpg)
As this happens, the point closest to the center of the selected white ball is selected as the tip of the cue stick( as we need the code to help play shot from cue stick and it needs to work correctly in that case). Then, a straight line is drawn from the selected point to the center of white ball, which is further extended and is also shown to reflect from the edges of tables two times, in order to show the predicted path of the white ball. This path is shown in the original image and is given as the output as shown below.
![alt text](https://github.com/himanshudua/Interactive-Pool-Table/blob/master/Pool%20Photos/3.jpg)
# Acknowledgements
[Soumay Gupta](https://github.com/iamguptaji?tab=repositories) for Guidance.
Google, Mathworks for solving issues.
