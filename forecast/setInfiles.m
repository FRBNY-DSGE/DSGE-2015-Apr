function [infile,numb,num] = setInfiles(indataType,...
                                        spath,nburn,nstate,nshocks,npara,nblocks,nsim)
for iDT = 1:length(indataType)
    switch indataType{iDT}
        case 'params'
            infile.params = [spath,'/params'];
            numb.params = nburn*npara*4;
            num.params = nblocks*nsim*npara*4;
        case 'TTT'
            infile.TTT = [spath,'/TTT'];
            numb.TTT = nburn*nstate^2*4;
            num.TTT = nblocks*nsim*nstate^2*4;
        case 'RRR'
            infile.RRR = [spath,'/RRR'];
            numb.RRR = nburn*nstate*nshocks*4;
            num.RRR = nblocks*nsim*nstate*nshocks*4;
        case 'zend'
            infile.zend = [spath,'/zend'];
            numb.zend = nburn*nstate*4;
            num.zend = nblocks*nsim*nstate*4;
    end
end
