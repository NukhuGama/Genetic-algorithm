function fitnessValue = Fitness2(x1, x2)
    h = Function2(x1,x2);
    
    fitnessValue = 2.^-h;
      
end