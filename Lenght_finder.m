function L = Lenght_finder(acc,Zout,inlet_window,latOut,Zoutlets,Zin)

%First horizontal distance: row and colum of outlet
[rOutlet, cOutlet] =  ind2sub(size(acc),Zout);

%Row and colum of inlets
[rInlet, cInlet] =  ind2sub(size(acc),inlet_window);

a = ((rOutlet - rInlet) * 450) * cos(degtorad(latOut)); % Longitudial length, corrected for latitude
b = (cOutlet - cInlet) * 450; % Latitudial length
hL = sqrt(a^2+b^2); %Horizontal distance

%Second correcting for vertical distance
L = sqrt(hL^2 + (Zoutlets - Zin)^2); %Total pipelenght

end