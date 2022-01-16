%Second order Lane-Emden equation into system of first-order ODEs

function [dy] = lane_emden(t,y,n)
    dy=zeros(2,1);
    dy(1)=y(2);
    dy(2)=-(y(2)*2/(t))-y(1)^n;
end