function [image, coinValue, debugImage] = CoinReader(fileName, format)
%COINREADER:[image, coinValue, debugImage] = CoinReader(fileName, format)  
%returns the total value of coins in the picture. Requirements:
%picture has a solid black background and no glare. Lighting must be
%decent.
%-----------------------------------------------
%The following algorithms are used to calculate the total coin value:
%
%Algorithm 1: Average Background Color:
%Sum the pixels in an edge boarder of a specified boarder width in two passes. 
%Pass one covers the vertical edges.
%Pass two covers the middle segments of the top and bottom edges, making
%sure not to count any square twice.
%
%Algorithm 2: Coin Finder:
%Iterate over the entire image, minus the edges taken out by boarders.
%Compare each pixel to the average background color (backgroundColor). If
%the difference is greater than the allowable threshold
%(pixelCompareThreshold), then step into the triangle
%algorithm. Step out of algorithm when all coins have been accounted for.
%
%Algorithm 3: Triangle Comparison:
%From the initial given position, travel down by TriangleDownStep. From
%that position, iteratively "draw" an isosceles triangle of height
%triangleSize. Take the average of the color values of all of these pixels,
%then compare it to the backgroundColor. If the difference is greater than 
%triangleAverageVsBoarderAverageThreshold, then step into the diameter
%algorithm.
%
%Algorithm 4: Diameter:
%From the initial given position, travel downwards iteratively while
%incrementing diameterCount. At each position, compare the current pixel to
%backgroundColor. If the difference is greater than pixelCompareThreshold,
%then append diameterCount to the end of the array "diameters." If the
%diameter count goes above 100 (e.g. something is wrong), then retry the
%algoritm from the pixel one down and one to the right. If that doesn't work,
%try one down and to the left of the original point. If that still doesn't 
%work, abandon ship on this triangle. Once (or if) the diameter
%has successfully been found, step into Coin Cleaner.
%
%Algorithm 5: Penny Finder
%Starting from the end of the Diameter algorithm (e.g. the last recorded
%picture is at the bottom center of a coin), travel diagonally up and to
%the right by a small, configurable value. Then, "draw" a box down and to the
%left such that the box is centered around the original pixel. Take the
%average red, green, and blue components of the pixels in this box.
%Considering these color components as a vector, subtract the blue value.
%Cut the matrix down to the first two elements. Repeat, this time
%subtracting the new green value. Finally, cut the resulting vector down to
%one element. This final value has been experimentally been determined to
%be 20-30 for gray coins and 32+ for pennnies.
%        
%
%Algorithm 6: Coin Cleaner
%From the initial given position, travel left to approximately the top left
%corner of the coin's bounding box. Then, replace all the
%pixels from the top left corner to the the approximate bottom right corner 
%of the coin's bounding box with backgroundColor.
%
%Algorithm 7: Diameter Compare:
%Take the coins larger than a given percentage of the largest coin's
%diameter and consider them as multiple coins of the same value. Remove
%these coins from the array and repeat. The algorithm ends when there are
%no more coins to consider.
%
%Algorithm 8: Final Sum
%From the array containing the amounts of each type of coin, map each value 
%A corrosponding coin value. Finally, take the sum of these values.


%Read in the image
img = imread(fileName, format);

%For now, we only need the grayscale. We'll need RGB version later.
imgColor=img;
img=rgb2gray(img);

%Storage Variables
tempCol=0;
tempRow=0;
tempSum=0;
singleTriangleAverageTotal=uint64(0);
diameterCount=0;
diameters=[];
pennyMatrix=[];
pennyMatricies=[];
abandonShip=0;
largestCoinDiameter=0;
numPennies=0;
coinValue=[];
coinDenominations=[1.00,0.25,0.05,0.10,0.01];

%Config Variables
pixelCompareThreshold=10;
triangleAverageVsBoarderAverageThreshold=30;                               %Too high messes with diameters!
boarderWidth=5;                                                            %Width of boarder used to calculate background color.
eliminationBoxConstant=1.07;                                               %Used to determine how large of a box to clear after finding a coin.
pennyAccuracy=10;
pennyCompareThreshold=32;
histTol=.07;


%Calculated Variables
[totalRows,totalCols]=size(img);
triangleSize=uint64(3+min(totalRows,totalCols)/100);                       %Defined as a fraction of the smaller image dimension
triangleBoarder=triangleSize;                                              %Extra boarder needed so as not to draw triangles outside of the image.
triangleDownStep=round(triangleSize/3);                                    %When triangles are placed, they are actually placed a little below the given position
triangleBottomBoarder = triangleDownStep+1;                                %Extra boarder for the bottom of the board needed because of the downwards travel of Triangle Comparison
triangleAveragePixels=triangleSize.^2;                                     %Sum of first n integers
backgroundColor=uint64(0);                                                 %Used as storage for the sum of boarder color values, then acted on to find the average background color
maxDiameter=max(totalCols,totalRows)./3;

%Debug
debug=0;
majorDebug=0;
graphicDebug=0;
triangleDebug=0;
oneRun=0;
hitIt=0;

%Set up debug image
debugImage=img;


%Algorithm 1: Vectorized background average color calculator. Can be rewritten as two, 
%double-nested, "for" loops.
backgroundColor=(uint64(sum(sum(img(:, horzcat(1:boarderWidth,totalCols-boarderWidth+1:totalCols))))) + uint64(sum(sum(img(horzcat(1:boarderWidth,totalRows-boarderWidth+1:totalRows),boarderWidth+1:totalCols-boarderWidth)))))./(2.*totalCols.*boarderWidth + 2.*totalRows.*boarderWidth - 4.*boarderWidth.^2);

if(debug)
    display(horzcat('backgroundColor: ', num2str(backgroundColor)));
end

%Algorithm 2: Non-vectorized Coin Finder. 
for(row=(triangleBoarder+boarderWidth):totalRows-(triangleBoarder+boarderWidth+triangleBottomBoarder+3)) %Extra 3 is to account for diameter drawing
    for(col=(triangleBoarder+boarderWidth):totalCols-(triangleBoarder+boarderWidth))
           
        if(abs(int8(img(row,col))- int8(backgroundColor))> pixelCompareThreshold)
           if(graphicDebug)
               debugImage(row,col) = 150;
           end
           if(majorDebug)
               display(horzcat('Row: ',num2str(row),', Col: ',num2str(col), ', Color: ', num2str(img(row,col)),', Diff: ',num2str(abs(int8(img(row,col))- int8(backgroundColor)))));
           end
            
           singleTriangleAverageTotal=0;
           
%Algorithm 3: Non-vectorized Triangle Calculator
           tempCol=col;
           tempRow=row+triangleDownStep;
           for(rr=1:triangleSize)
               for(cc=1:(2*rr-1))
                   tempSum=tempSum+img(tempRow+rr,tempCol-rr+cc);
                   
                   singleTriangleAverageTotal = singleTriangleAverageTotal + uint64(img(tempRow+rr,tempCol-rr+cc));
                   if(graphicDebug)
                   debugImage(tempRow+rr,tempCol-rr+cc)=200;
                   end
                   
               end
           end

           singleTriangleAverageTotal = singleTriangleAverageTotal./triangleAveragePixels; 
     
%Algorithm 4: Non-vectorized Diameter Calculator
           if(abs(int8(singleTriangleAverageTotal)-int8(backgroundColor)) > triangleAverageVsBoarderAverageThreshold )
               if(majorDebug)
                   display(horzcat('Row: ',num2str(row),', Col: ',num2str(col), ', Color Average: ', num2str(singleTriangleAverageTotal),', Diff: ',num2str(abs(int8(singleTriangleAverageTotal)-int8(backgroundColor)))));
               end
               tempCol=col;
               tempRow=row;
               if(graphicDebug)
                   debugImage(tempRow-triangleDownStep-1,tempCol)=255;
               end
               while((~(tempRow>totalRows))&&(abs(int8(img(tempRow,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+1,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+2,tempCol))-int8(backgroundColor)) >pixelCompareThreshold))
                    if(graphicDebug)
                        debugImage(tempRow,tempCol)=255;
                    end
                    tempRow=tempRow+1;
                    diameterCount=diameterCount+1;
               %This algorithm should end here. However, in case of pixel glitches, the following code attempts to perform the above algorithm on a different starting pixel    
                    if(diameterCount==maxDiameter) %Problem
                        diameterCount=0;
                        tempCol=col+1;
                        tempRow=row+1;
                        while(abs(int8(img(tempRow,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+1,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+2,tempCol))-int8(backgroundColor)) >pixelCompareThreshold)
                            tempRow=tempRow+1;
                            diameterCount=diameterCount+1;
                            if(diameterCount==maxDiameter) %Bigger Problem
                                diameterCount=0;
                                tempCol=col-1;
                                tempRow=row+1;
                                while(abs(int8(img(tempRow,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+1,tempCol))-int8(backgroundColor)) >pixelCompareThreshold || abs(int8(img(tempRow+2,tempCol))-int8(backgroundColor)) >pixelCompareThreshold)
                                    tempRow=tempRow+1;
                                    diameterCount=diameterCount+1;
                                    if(diameterCount==maxDiameter) %Uh oh
                                        display('Abandon Ship!');
                                        diameterCount=0;
                                        abandonShip=1;
                                        break;
                                    end
                                end
                            end
                            if(abandonShip)
                                break;
                            end
                        end     
                    end
                    if(abandonShip)
                        abandonShip=0;
                        break;
                    end
               end                
               if(diameterCount>10)
                   
%Algorithm 5: Penny Finder (Part A)
                    tempRow=tempRow-(diameterCount./2);
                    tempCol=tempCol-(round(pennyAccuracy./2));
                    pennyMatrix=imgColor(tempRow:tempRow+pennyAccuracy,tempCol:tempCol+pennyAccuracy,1:3);
                    pennyMatrix=[mean(mean(pennyMatrix(:,:,1))), mean(mean(pennyMatrix(:,:,2))),mean(mean(pennyMatrix(:,:,3)))];
                    pennyMatricies=vertcat(pennyMatricies,pennyMatrix);
                   
                   diameters = horzcat(diameters, diameterCount);
                   if(debug)
                       display(horzcat('Diameters: ', num2str(diameters)));
                   end
                   
%Algorithm 6: Coin Cleaner
                   tempRow=row-(diameterCount*(eliminationBoxConstant-1));
                   tempCol=col-(diameterCount.*(eliminationBoxConstant-1) + diameterCount./2);
                   rowEnd=round(tempRow + diameterCount.*(2.*eliminationBoxConstant - 1));
                   colEnd=col+(diameterCount.*(2.*(eliminationBoxConstant-1)) + diameterCount./2);
                   
                   if(tempRow<(boarderWidth+1))
                       tempRow=boarderWidth+1;
                   end
                   if(tempRow>(totalRows-triangleBottomBoarder-triangleBoarder-boarderWidth))
                       tempRow=totalRows-triangleBottomBoarder-triangleBoarder-boarderWidth;
                   end
                   
                   if(tempCol<(boarderWidth+triangleBoarder))
                       tempCol=boarderWidth+triangleBoarder;
                   end
                   if(tempCol>(totalCols-triangleBoarder-boarderWidth))
                       tempCol=totalCols-triangleBoarder-boarderWidth;
                   end            
                   
                   if(rowEnd<(boarderWidth+1))
                       rowEnd=boarderWidth+1;
                   end
                   if(rowEnd>(totalRows-triangleBottomBoarder-triangleBoarder-boarderWidth))
                       rowEnd=totalRows-triangleBottomBoarder-triangleBoarder-boarderWidth;
                   end
                   
                   if(colEnd<(boarderWidth+triangleBoarder))
                       colEnd=boarderWidth+triangleBoarder;
                   end
                   if(colEnd>(totalCols-triangleBoarder-boarderWidth))
                       colEnd=totalCols-triangleBoarder-boarderWidth;
                   end
                   
                   if(graphicDebug)
                        debugImage(tempRow,tempCol)=150;
                        debugImage(rowEnd,colEnd)=150;
                   end
                   

                   img(tempRow:rowEnd,tempCol:colEnd)=backgroundColor;
                   %debugImage(tempRow:rowEnd,tempCol:colEnd)=200;


                   if(debug)
                       display('Coin removed');
                   end                                   
                   diameterCount=0;         
               end
               diameterCount=0;
           end
            %End of diameters
            
            if(oneRun)
               hitIt=1;
               break;
            end
           
        end
    end
    
    if(oneRun && hitIt)
        break;
    end
    
end

diameterSize=size(diameters);
values=zeros(1,diameterSize(2));
largestCoinDiameter=max(diameters);
pennyMatriciesSize=size(pennyMatricies);

%Algorithm 5: (Part B)
for(r=1:pennyMatriciesSize(1))
pennyMatricies(r,:)=pennyMatricies(r,:)-min(pennyMatricies(r,:));
end
pennyMatricies=pennyMatricies(:,1:2);
for(r=1:pennyMatriciesSize(1))
pennyMatricies(r,:)=pennyMatricies(r,:)-min(pennyMatricies(r,:));
end
pennyMatricies=pennyMatricies(:,1);
pennyMatricies=pennyMatricies>pennyCompareThreshold;
pennyMatricies=rot90(pennyMatricies);
numPennies=sum(pennyMatricies);
pennyMatricies=pennyMatricies<1;
diameters=diameters.*pennyMatricies;
diameters=diameters(diameters~=0);

%Algorithm 7: Diameter Compare
index=1;
while(length(diameters)>0)
    coinValue(index)=sum(diameters>(1-histTol).*max(diameters));
    diameters(diameters>(1-histTol).*max(diameters))=[];
    index=index+1;    
end

%Algorithm 8: Final Sum
coinValue=horzcat(coinValue, numPennies);
coinValue=coinValue.*(coinDenominations(6-length(coinValue):5));
             
if(graphicDebug)
    imshow(debugImage);
end

coinValue=sum(coinValue);
image=img;


end