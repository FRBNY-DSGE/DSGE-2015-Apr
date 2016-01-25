function para = invtrans(para,trspec);
% this procedure transforms variables from model to max
% Note that the max parameter do not contain para 8 and 9

  npara = size(para,1);

  for i = 1:npara

    a = trspec(i,2);
    b = trspec(i,3);
    c = trspec(i,4);
    
    if trspec(i,1) == 1;
      cx = 2*(para(i)-(a+b)/2)/(b-a);
      para(i) = (1/c)*cx/sqrt(1-cx^2);
    elseif trspec(i,1) == 2;
      para(i) = b + (1/c)*log(para(i)-a);
    end

  end
