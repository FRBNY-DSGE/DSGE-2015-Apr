%                       OVERVIEW
% momTable.m: tabulates parameter moments in 3 laTex tables: 
%             (1) For our MAIN parameters, a list of the prior means, prior
%                 standard deviations, posterior means and 90% bands for posterior
%                 draws.
%             (2) For LESS IMPORTANT parameters, a list of the prior means, prior
%                 standard deviations, posterior means and 90% bands for posterior
%                 draws.
%             (3) A list of prior means and posterior means
%
%             momTable.m is called on in mom.m
% 
%                       INPUTS
% theta: (ndraws x npara) matrix holding the posterior draws
%        (stored in mhpara*, output from gibb.m)
% post: (ndraws x 1) matrix holding the posterior value 
%       (stored in post*, output from gibb.m)
% priotheta: (ndraws x npara) matrix holding the prior draws 
%            (stored in mhNIpriopara*, output from prio.m)
%
%                       OUTPUTS
% outfile0: This refers to a file (mhpara*Mean) that holds the mean
%           parameter across draws from the posterior.
% outfile01: This refers to a latex table (*Mom_MainParams*) that lists the
%            moments for important parameters. 
% outfile02: This refers to a latex table (*Mom_PeriphParams*) that lists the
%            moments for less important parameters. 
% outfile03: This refers to a latex table (*PrioPostMean*) that lists
%            the prior and posterior means 

%% Compute mean, st dev, and 90% bands across draws from the posterior

% pmean and pstdd (prior mean and st dev are computed in
% initializePrograms.m

thetahat = mean(theta,1)';
thetasig = cov(theta,1);
thetabands = hpdint(theta,percent,1)';

%% Save posterior mean

outfile0 = [spath,'/params_mean'];
fid0 = fopen(outfile0,'w');
fwrite(fid0,thetahat','single');
fclose(fid0);

%% Open outfiles for writing

% Moments for main parameters
outfile01 = [fpath,'Mom_MainParams.tex'];
fid01 = fopen(outfile01,'w');

% Moments for peripheral parameters
outfile02 = [fpath,'Mom_PeriphParams.tex'];
fid02 = fopen(outfile02,'w');

% Parameter mean across prior draws and posterior draws
outfile03 = [fpath,'PrioPostMean.tex'];
fid03 = fopen(outfile03,'w');

%% Create variables used for writing to the output table

colnames = {'Parameter ',...
    'Prior Mean ',...
    'Prior Stdd ',...
    'Post Mean ',...
    [num2str(100*percent),'\% {\tiny Lower Band}'],...
    [num2str(100*percent),'\% {\tiny Upper Band}'] };

outmat = [pmean,pstdd,thetahat,thetabands]; % prior mean and standard devaition, posterior mean and bands
outmat2 = [pmean thetahat]; % prior mean and posterior mean

%% Write to Table 1: Prior mean, st dev and posterior mean, bands for IMPORTANT parameters

fprintf(fid01,'\n \\documentclass[12pt]{article}');
fprintf(fid01,'\n \\usepackage[dvips]{color}');
fprintf(fid01,'\n \\begin{document}');
fprintf(fid01,'\n \\pagestyle{empty}');
fprintf(fid01,'\n \\begin{table}[h] \\centering');
fprintf(fid01,'\n \\caption{Parameter Estimates}');
fprintf(fid01,'\n \\vspace*{.5cm}');
fprintf(fid01,'\n {\\small \n');
fprintf(fid01,'\n \\begin{tabular}{lllllll}\\hline \n');

fprintf(fid01,' %4.99s  & ',colnames{:});

% Keep track of the indices for these important parameters within para_names
important_para = [];
for i = 1:length(para_names)
    if isempty(regexp(para_names{i},'rho_','once')) && ...
            isempty(regexp(para_names{i},'zeta_','once')) && ...
            isempty(regexp(para_names{i},'psi_','once')) && ...
            isempty(regexp(para_names{i},'nu_l','once')) && ...
            isempty(regexp(para_names{i},'pi^*','once')) && ...
            (mspec ~= 16 || isempty(regexp(para_names{i},'u^*','once')))
        continue;
    end
    if strcmp(para_names{i},'\rho_{\chi}') || (isequal(subspec,7) && strcmp(para_names{i},'\rho_{b}'))
        continue;
    end
    fprintf(fid01,'\\\\ \n $%4.99s$ &  ',para_names{i});
    fprintf(fid01,' %8.3f & ',outmat(i,:));
    important_para = [important_para, i];
end

fprintf(fid01,'\\\\  \\hline');
fprintf(fid01,'\n \\end{tabular}}');
fprintf(fid01,'\n \\end{table} \n\n');
fprintf(fid01,'\n \\end{document}');
fclose(fid01);

%% Write to Table 2: Prior mean, st dev and posterior mean, bands for OTHER parameters

fprintf(fid02,'\n \\documentclass[12pt]{article}');
fprintf(fid02,'\n \\usepackage[dvips]{color}');
fprintf(fid02,'\n \\begin{document}');
fprintf(fid02,'\n \\pagestyle{empty}');
fprintf(fid02,'\n \\begin{table}[h] \\centering');
fprintf(fid02,'\n \\caption{Parameter Estimates}');
fprintf(fid02,'\n \\vspace*{.2cm}');
fprintf(fid02,'\n {\\small \n');
fprintf(fid02,'\n \\begin{tabular}{lllllll}\\hline \n');

fprintf(fid02,' %4.99s  & ',colnames{:});

% Counter for parameters to track length of table, number of tables in
% excess of default '1'
other_para = 1;
table_count = 0;
for i = 1:length(para_names)
    if ismember(i,important_para)
        continue;
    end
    
    if mod(other_para,25) == 0 && ~(i == length(para_names))
        fprintf(fid02,'\\\\  \\hline');
        fprintf(fid02,'\n \\end{tabular}}');
        fprintf(fid02,'\n \\end{table} \n\n');
        fprintf(fid02,'\n \\end{document}');
        
        table_count = table_count + 1;
        fclose(fid02);
        
        filename = sprintf('Mom_PeriphParams_%i.tex',table_count);
        outfile02 = [fpath,filename];
        fid02 = fopen(outfile02,'w');

        fprintf(fid02,'\n \\documentclass[12pt]{article}');
        fprintf(fid02,'\n \\usepackage[dvips]{color}');
        fprintf(fid02,'\n \\begin{document}');
        fprintf(fid02,'\n \\pagestyle{empty}');

        fprintf(fid02,'\n \\begin{table}[h] \\centering');
        fprintf(fid02,'\n \\caption{Parameter Estimates}');
        fprintf(fid02,'\n \\vspace*{.2cm}');
        fprintf(fid02,'\n {\\small \n');
        fprintf(fid02,'\n \\begin{tabular}{lllllll}\\hline \n');

        fprintf(fid02,' %4.99s  & ',colnames{:});
    end
    
    fprintf(fid02,'\\\\ \n $%4.99s$ &  ',para_names{i});
    fprintf(fid02,' %8.3f & ',outmat(i,:));
    other_para = other_para + 1;
end
fprintf(fid02,'\\\\  \\hline');
fprintf(fid02,'\n \\end{tabular}}');
fprintf(fid02,'\n \\end{table} \n\n');
fprintf(fid02,'\n \\end{document}'); 
fclose(fid02);
clear parasimm postsim;

%% Write to Table 3: Prior mean and posterior mean for all parameters

fprintf(fid03,'\n \\documentclass[12pt]{article}');
fprintf(fid03,'\n \\usepackage[dvips]{color}');
fprintf(fid03,'\n \\begin{document}');
fprintf(fid03,'\n \\pagestyle{empty}');
fprintf(fid03,'\n \\begin{table}[h] \\centering');
fprintf(fid03,'\n \\caption{Parameter Estimates: Prior and Posterior Mean}');
fprintf(fid03,'\n \\vspace*{.5cm}');
fprintf(fid03,'\n \\begin{tabular}{ccc');
fprintf(fid03,'}\\hline \n');
fprintf(fid03,' Parameter & Prior ');
fprintf(fid03,'\\\\ \\hline \n');

for i = 1:length(para_names)
    fprintf(fid03,' \n $%4.99s$   ',para_names{i});
    fprintf(fid03,' & %8.3f  ',outmat2(i,:));
    fprintf(fid03,'\\\\ ');
end
fprintf(fid03,'\\\\  \\hline');
fprintf(fid03,'\n \\end{tabular}');
fprintf(fid03,'\n \\end{table} \n\n');
fprintf(fid03,'\n \\end{document}'); 
fclose(fid03);

close all;

fprintf('Tables are in %s \n',fpath);
