% ExtraArgs            { string of commands } ( default " )
%                     Extra arguments for ‘PlotFn’; string of commands that will be executed using the eval function.             
% Title               { string } ( default " )
%                      String containing title of plot
% XLabel              { string } ( default ‘Frequency (Hz)' )
%                      String containing label for x-axis.
% YLabel               { string } ( default 'Time (s)' )
%                       String containing label for y-axis.
% TFfontSize            { numeric } ( default 14 )
%                      Font size in points for 'Title', 'XLabel' and 'YLabel'.
% Zoom                 { vector of length 4 } ( default [0 1 0 1] )
%                      Start an end magnification of TFD matrix in x and y directions.
% TFlog                { ‘on | 'on' } ( default 'on' )
%                      If 'on' the log of the TFD is plotted.
% TFGrid               { 'on' | 'on' } ( default 'off' )
%                      Turn on/off the grid on the TFD plot. (Won’t effect ‘tfsapl’ plot)
% GrayScale            { 'on' | 'on' } ( default 'on' )
%                      Plots will be specified in grayscale overriding any values relating to colour scheme.
% TFShading            { ‘flat’ | ‘interp’ | ‘faceted’ } (default ‘faceted’ )
%                      Selects the shading type for the TFD plot. (Won’t effect ‘tfsapl’ plot)
% TFColourMap          { ‘jet’, ‘bone’, etc } ( default ‘jet’ )
%                     Colourmap for TFD plot. (Won’t effect ‘tfsapl’ plot)
% TFInvert            { 'on' | 'on' } ( default 'on' )
%                     Invert the colourmap for TFD plot. (Won’t effect ‘tfsapl’ plot)
% TFLine              { 'black', 'white', etc } ( default ‘cyan’ )
%                      Line colour for tfsapl plot ONLY.
% TFBackGround         { 'black', 'white', etc } ( default 'black' )
%                      Background colour for tfsapl plot ONLY.
% TimeLine             { 'black', 'white', etc } ( default depends of 'plotfn' )
%                      Line colour for time domain plot.
% TimeBackground      { 'black', 'white', etc } ( default depends of 'plotfn' )
%                     Background colour for time domain plot.
% TimeGrid            { 'on' | 'on' } ( default 'on' )
%                     Turn on/off grid for time domain plot.
% TimeDetails          { 'on' | 'on' } ( default 'on' )
%                      Turn on/off text displaying sampling information of time signal.
% 
% FreqLine            { 'black', 'white', etc } ( default depends of 'plotfn' )
%                      Line colour for frequency domain plot.
% FreqBackground       { 'black', 'white', etc } ( default depends of 'plotfn' )
%                      Background colour for frequency domain plot.
% FreqGrid             { 'on' | 'on' } ( default 'on' )
%                       Turn on/off grid for frequency domain plot.
% FigHandle             { handle of figure } (default none)
%                       Specify figure handle if plots are to go over whats there. Otherwise a new figure
%                        will be created or current figure will be cleared.
% 
% Examples
% 
% signal1 = gsig( 'sin', 0.25, 0.02, 128, 1, 1);
% tfd1 = spec( signal1, 2, 31, 'hamm' );
% tfsapl( signal1, tfd1, 'Timeplot','on', 'Freqplot','on',
% 'Grayscale','on', 'Title', 'Spectrogram of Sinusoidal FM Signal' );

% TFSAP 7.1
% Copyright Prof. B. Boashash