function result = Function2(x1,x2)
    result = -cosd(x1)*cosd(x2)  * exp(-((x1-pi).^2) - ((x2-pi).^2));
end
