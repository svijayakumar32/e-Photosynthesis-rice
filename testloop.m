list = ([1,4,0.5;4,6,2;5,4,1;23,20,25;89,87,90;3,2,1;5,7,2;31,33,25])';
for j = 1:3
if list(1)<1 
    disp('list item 1 is below the limit.')
           continue
elseif list(2)<4
    disp('list item 2 is below the limit.')
           continue
elseif list(3)<5
    disp('list item 3 is below the limit.')
           continue
elseif list(4)<23
    disp('list item 4 is below the limit.')
           continue 
elseif list(5)<89
    disp('list item 5 is below the limit.')
           continue
elseif list(6)<3
    disp('list item 6 is below the limit.')
           continue
elseif list(7)<5
    disp('list item 7 is below the limit.')
           continue
elseif list(8)<31
    disp('list item 8 is below the limit.')
           continue
else 
end   
end
