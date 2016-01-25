%                           OVERVIEW
% mom.m: Computes moments, tabulates parameter moments, and plots parameter
%        draws from the prior and posterior distribution. 
%           
%                       IMPORTANT VARIABLES
% PLOTDRAWS: Specify whether or not you want to plot the parameter draws
%           (runs momPlot.m)
% percent: In the output laTex table, we report <percent> bands for
%          parameter draws from the posterior 
%       
%                            INPUTS
% infile1: This refers to a an output file (mhparam*) from gibb.m, 
%          which holds draws from the posterior.  
% infile4: This refers to a an output file (post*) from gibb.m, which
%          holds the value of the posterior, for each posterior draw. 
% theta: (ndraws x npara) matrix holding the posterior draws
%        (stored in mhpara*, output from gibb.m)
% post: (ndraws x 1) matrix holding the posterior value 
%       (stored in post*, output from gibb.m)
% priotheta: (ndraws x npara) matrix holding the prior draws 
%            (stored in mhNIpriopara*, output from prio.m) 
%            (optional input)
%
%                           OUTPUTS
% From momTable.m: 1: (mhpara*Mean): holds the mean parameter across
%                     draws from the posterior.
%                  2: (*Mom_MainParams*): laTex table that lists the
%                     moments for important parameters. 
%                  3: (*Mom_PeriphParams*): laTex table that lists the
%                     moments for less important parameters. 
%                  4: (*PrioPostMean*) that lists the prior and
%                     posterior means 

close all;
keepVars; 
initializePrograms;

%% Settings
Verbose = 0;
percent = .90; % Bands

%% Input files

% Posterior draws (output from gibb.m)
infile1 = [spath, 'params'];
num1 = nblocks*nsim*npara*4;  % number of bites per file
numb1 = nburn*npara*4;        % number of bites to dicard
disp(['Infile: ' infile1]);

% Posterior values file (output from gibb.m)
infile4 = [spath,'/post'];
num4 = nblocks*nsim*1*4;            
numb4 = nburn*1*4;                  

%% Load in theta (postrior draws), priotheta (prior draws), and post (posterior value)
fid1 = fopen(infile1,'r');		
status = fseek(fid1,numb1,'bof');

fid4 = fopen(infile4,'r');		
status = fseek(fid4,numb4,'bof');

theta = [];
priotheta = [];
post = [];

% Read blocks of size nsim
while ( ftell(fid1) < num1 )
    
    thetadd = fread(fid1,[npara,nsim],'single')';
    theta = [theta;thetadd];
    clear thetadd;
    
    postadd = fread(fid4,[nsim,1],'single');
    post = [post;postadd];
    clear postadd;
    
end;

fclose(fid1);
fclose(fid4);


ndraws = size(theta,1);

%% Produce table (Tex) of moments 
momTable;

infile1 = [spath,'params'];
num1 = nblocks*nsim*npara*4;    % number of bites per file
numb1 = nburn*npara*4;          % number of bites to dicard
fid1 = fopen(infile1,'r');
status = fseek(fid1,numb1,'bof');

theta = [];

while ( ftell(fid1) < num1 ) % Read blocks of size nsim
    theta_add = fread(fid1,[npara,nsim],'single')';
    theta = [theta; theta_add];
end

cov_theta = cov(theta);

outfile100 = [spath,'/cov'];
fid100 = fopen(outfile100,'w');

fwrite(fid100,cov_theta(:),'single');

fclose('all');
    
