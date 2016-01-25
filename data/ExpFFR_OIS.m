
%% Loads in expected FFR derived from OIS quotes


function [ExpFFR,peachdata_FFR] = ExpFFR_OIS(nant,antlags,psize,zerobound,peachflag)

% 2008-Q4 expectations (from Jan-2009 BCFF survey, conducted mid-/end- of Dec-2008)
ExpFFR(1,:) = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.1 2.2 2.4]; 

% % 2009-Q1 expectations (from Apr-2009 BCFF survey, conducted mid-/end- of Mar-2009)
ExpFFR(2,:) = [0.2 0.3 0.5 0.6 0.8 1.1 1.3 1.5 1.7 2.0 2.2 2.4 2.5];

% 2009-Q2 expectations (from Jul-2009 BCFF survey, conducted mid-/end- of Jun-2009)
ExpFFR(3,:) = [0.2 0.4 0.7 1.1 1.4 1.8 2.2 2.6 2.9 3.2 3.5 3.7 4.0];

% 2009-Q3 expectations (from Oct-2009 BCFF survey, conducted mid-/end- of Sep-2009)
ExpFFR(4,:) = [0.2 0.4 0.7 1.0 1.4 1.8 2.1 2.4 2.7 2.9 3.1 3.3 3.6];

% 2009-Q4 expectations (from Jan-2010 BCFF survey, conducted mid-/end- of Dec-2009)
ExpFFR(5,:) = [0.2 0.4 0.8 1.2 1.6 2.0 2.4 2.7 3.0 3.3 3.5 3.7 4.0];

% 2010-Q1 expectations (from Apr-2010 BCFF survey, conducted mid-/end- of Mar-2010)
ExpFFR(6,:) = [0.2 0.4 0.6 0.9 1.2 1.6 1.9 2.3 2.6 2.9 3.1 3.4 3.6];

% 2010-Q2 expectations (from Jul-2010 BCFF survey, conducted mid-/end- of Jun-2010)
ExpFFR(7,:) = [0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.7 1.9 2.2 2.2];

% 2010-Q3 expectations (from Oct-2010 BCFF survey, conducted mid-/end- of Sep-2010)
ExpFFR(8,:) = [0.2 0.2 0.2 0.3 0.4 0.5 0.6 0.8 0.9 1.1 1.3 1.5 1.6];

% 2010-Q4 expectations (from Jan-2011 BCFF survey, conducted mid-/end- of Dec-2010)
ExpFFR(9,:) = [0.2 0.2 0.3 0.4 0.6 0.8 1.0 1.3 1.6 1.8 2.1 2.4 2.6];

% 2011-Q1 expectations (from Apr-2011 BCFF survey, conducted mid-/end- of Mar-2011)
ExpFFR(10,:) = [0.1 0.2 0.4 0.7 0.9 1.2 1.5 1.8 2.1 2.4 2.7 2.9 3.2];

% 2011-Q2 expectations (from Jul-2011 BCFF survey, conducted mid-/end- of Jun-2011)
ExpFFR(11,:) = [0.1 0.2 0.3 0.4 0.5 0.7 0.9 1.2 1.4 1.7 2.0 2.2 2.4];

% 2011-Q3 expectations (from Oct-2011 BCFF survey, conducted mid-/end- of Sep-2011)
ExpFFR(12,:) = [0.1 0.1 0.1 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.8 0.9 1.0];

% 2011-Q4 expectations (from 12/30/2011 OIS data)
ExpFFR(13,:) = [0.1 0.1 0.1 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.7 0.9 1.0];
 
% 2012-Q1 expectations (from 3/30/2012 OIS data)
ExpFFR(14,:) = [0.1 0.1 0.2 0.2 0.3 0.3 0.4 0.5 0.6 0.7 0.9 1.0 1.2];
%  
% % 2012-Q2 expectations (from 6/30/2012 OIS data)
ExpFFR(15,:) = [0.2 0.2 0.2 0.2 0.2 0.3 0.3 0.4 0.4 0.5 0.5 0.6 0.8];
 
% % 2012-Q3 expectations (from 9/30/2012 OIS data)
ExpFFR(16,:) = [0.1 0.1 0.1 0.1 0.1 0.1 0.2 0.2 0.3 0.4 0.5 0.6 0.7];
%
% % 2012-Q4 expectations (from 12/31/2012 OIS data)
ExpFFR(17,:) = [0.1 0.1 0.1 0.1 0.1 0.2 0.2 0.3 0.3 0.4 0.5 0.6 0.8];

% 2013-Q1 expectations (from 3/29/2013 OIS data)
ExpFFR(18,:) = [0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.3 0.3 0.4 0.5 0.6 0.7];

% 2013-Q2 expectations (from 6/29/2013 OIS data)
ExpFFR(19,:) = [0.1 0.2 0.2 0.3 0.4 0.5 0.7 0.8 1.0 1.2 1.3 1.5 1.7];

% 2013-Q3 expectations (from 10/3/2013 OIS data)
ExpFFR(20,:) = [0.1 0.1 0.1 0.2 0.2 0.3 0.4 0.6 0.7 0.9 1.1 1.3 1.5]; 

% 2013-Q4 expectations (from 12/31/2013 OIS data)
ExpFFR(21,:) = [0.1 0.1 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.7 1.9];

% 2014-Q1 expectations (from 3/31/2014 OIS data)
ExpFFR(22,:) = [0.1 0.1 0.1 0.2 0.4 0.5 0.7 0.9 1.1 1.4 1.7 1.9 2.2];

% 2014-Q2 expectations (from 6/30/2014 OIS data)
ExpFFR(23,:) = [0.1 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.7 1.9 2.2];

% 2014-Q3 expectations (from 9/30/2014 OIS data)
ExpFFR(24,:) = [0.1 0.1	0.2 0.4	0.6 0.9	1.1 1.4	1.6 1.9	2.1 2.3	2.5];

% 2014-Q4 expectations (from 12/31/2014 OIS Data)
ExpFFR(25,:) = [0.1 0.2 0.4 0.5 0.7 0.9 1.2 1.4 1.6 1.8 1.9 2.1 2.3];
               
if nant <= size(ExpFFR,2)
    ExpFFR = ExpFFR(:,1:nant);
else
    ExpFFR = [ExpFFR,NaN(antlags+1,nant-size(ExpFFR,2))];
end

if peachflag
       xpeachdata_FFR(1,:) = [0.16 0.28 0.44 0.61 0.79 0.98 1.17 1.34 1.51 1.67 NaN];
       
% This should start with the OIS value for the first quarter AFTER the
% one you are forecasting. So if you have data until Q2 and you are
% forecating Q3, the first element of the array should be for Q4. The
% conditional peachdata file will supply the Q3 value.

    if nant <= size(xpeachdata_FFR,2)
        for t = 1:psize
            peachdata_FFR(t,:) = [xpeachdata_FFR(t,1:nant-t),NaN(1,t)];
        end
    else
        peachdata_FFR = [xpeachdata_FFR,NaN(psize,nant-size(xpeachdata_FFR,2))];
    end
    
    
else
    peachdata_FFR = [];
end
