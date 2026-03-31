function [ interHD ] = interpolateHDData( HD, TPOS, dt, t1, t2 )

scalingVector = 0:dt:((t2-1) + (1-dt)); 

interHD = mod(interp1(TPOS, unwrap(HD), scalingVector,'spline'), 2*pi);

i1 = (t1)/dt + 1;
i2 = t2/dt;

interHD = interHD(i1:i2);

end


