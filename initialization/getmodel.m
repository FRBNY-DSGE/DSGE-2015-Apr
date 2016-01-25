function [TTT,RRR,zend] = getmodel(TTTsim,RRRsim,zendsim,j2,nstate,nshocks)

TTT = reshape(TTTsim(j2,:),nstate,nstate);
RRR = reshape(RRRsim(j2,:),nstate,nshocks);
zend = reshape(zendsim(j2,:),nstate,1);