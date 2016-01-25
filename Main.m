%% Main.m
%% This script produces all output associated with the FRBNY DSGE model. 
% 
% Copyright Federal Reserve Bank of New York.  You may reproduce, use, modify,
% make derivative works of, and distribute and this code in whole or in part
% so long as you keep this notice in the documentation associated with any
% distributed works.   Neither the name of the Federal Reserve Bank of New
% York (FRBNY) nor the names of any of the authors may be used to endorse or
% promote works derived from this code without prior written permission.
% Portions of the code attributed to third parties are subject to applicable
% third party licenses and rights.  By your use of this code you accept this
% license and any applicable third party license.  
% 
% THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR
% CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT
% LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT,
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT
% THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID.  FRBNY IS NOT, UNDER
% ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR
% IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT
% LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR
% EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY,
% TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH
% DAMAGES OR LOSS IS FORESEEABLE.

%% Initialization
clear
close all
spec_990  % sets important variables and flags
set_paths % adds necessary paths

%% Estimation
% In this stage we draw from the distribution for the parameters. The modal
% parameters as well as the draws of parameters, are outputted in the /save
% folder.
gibb


%% Forecasting
% Here we produce forecasts for our observable variables, one associated
% with each draw of the parameters, and saves them in the /save folder. 
forecast_parallel_est_ant

%% Plotting
forplot          % produces series to be plotted
plotPresentation % produces plots of series outputted from forplot, which 
                 % are saved in the /graphs folder
