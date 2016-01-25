function [ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = getmeasur(mspec,TTT,RRR,valid,params,nvar,nlags,npara,coint,cointadd, nant)

if exist('nant','var')
    eval(strcat('[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT,RRR,valid,params,nvar,nlags,',num2str(mspec),',npara,coint,cointadd,nant);'));
else
    eval(strcat('[ZZ,DD,DDcointadd,QQ,EE,MM,retcode] = measur',num2str(mspec),'(TTT,RRR,valid,params,nvar,nlags,',num2str(mspec),',npara,coint,cointadd);'));
end
