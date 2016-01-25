function para = trans(para,trspec);
% this procedure transforms variables from max to model

  npara = size(para,1);

  for i = 1:npara

    a = trspec(i,2);
    b = trspec(i,3);
    c = trspec(i,4);
    
    if trspec(i,1) == 1;
      
      para(i) = (a+b)/2 + 0.5*(b-a)*c*para(i)/sqrt(1+c^2*para(i)^2);
      
      %      para(i) = a + (b-a)*cdfn(c*para(i));
      
    elseif trspec(i,1) ==2;
      para(i) = a + exp(c*(para(i)-b));
    end

  end