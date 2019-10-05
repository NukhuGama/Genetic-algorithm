function fitnessValue = Fitness1(x1,x2)
    h = Function1(x1,x2);
    fitnessValue = 2.^-h;
      
end