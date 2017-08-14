classdef recon < handle
%RECON - Create a reconstruction object
% The recon class is used to train and test linear reconstruction
% filters for an isetbio object.
% 
% rc = recon(isetbioObject, 'PARAM1', val1, 'PARAM2', val2,...) creates the
% recon object. Optional parameter name/value pairs are listed below.
% 
% The recon object computes the optimal linear decoder for an isetbio
% object based on its response and the image stimulus. The recon object was
% primarily tested on RGC responses, but may also be used for cone
% absorptions, cone current and bipolar responses.
%
%       pRecon.innerRetina = innerRetina;
%       pRecon.mosaicFile = 'mosaicAll_1.mat';
%       pRecon.filterFile = 'filter_mosaicAll_1.mat';
%       reconHealthy = recon(pRecon);
% 
% Input: the response of an isetbio object (absorptions, current, or
% spikes) and the image input.
% 
% Output: the decoding filters, which can be applied to a set of responses
% to generate a stimulus reconstruciton.
% 
%     cellLocation;                    % location of bipolar RF center
% 
%  ISETBIO wiki: <a href="matlab:
%  web('https://github.com/isetbio/isetbio/wiki/reconstruction','-browser')">reconstruction</a>.
%   
% 5/2016 JRG (c) isetbio team

%% Define object
% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties.
properties (SetAccess = protected, GetAccess = public)
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

properties (Access = public)    
    %MOSAICFILE - name of file used to generate inner retina mosaic
    mosaicFile;
    
    %BUILDFILE - name of file with original stim and resp values
    buildFile;
    
    %STIMFILE - the processed stimulus file for training
    stimFile
    
    %RESPFile - the processed response (spikes) file for training
    respFile
    
    %FILTERFILE - name of filter trained on inner retina mosaic
    filterFile;
    
    %WINDOWSIZE - size of the temporal decoding window
    windowSize;
    
    %SHIFTTIME - the length of the temporal shift for the decoding window
    shiftTime;
    
    %PERCENTSV - the percentage of singular values kept for the decoding
    %   filters
    percentSV;

end

% Public methods
methods
    
    % Constructor
    function obj = recon(varargin)     
        % Initialize the recon class
        %   reconHealthy = recon('mosaicHealhty','naturalSceneBuild');
        
        p = inputParser;
        addParameter(p,  'mosaicFile','',@ischar);% ['mosaic_' num2str(round(cputime))], @ischar);
        addParameter(p, 'buildFile',  '',@ischar);
        addParameter(p, 'stimFile',   '',@ischar);
        addParameter(p, 'respFile',   '',@ischar);
        addParameter(p, 'filterFile', '',@ischar);
        addParameter(p, 'windowSize', 1, @isnumeric);
        addParameter(p, 'shiftTime',  0, @isnumeric);
        addParameter(p, 'percentSV',  100, @isnumeric);        
        p.KeepUnmatched = true;
        p.parse(varargin{:});  
        
        obj.mosaicFile = p.Results.mosaicFile;
        obj.buildFile = p.Results.buildFile;
        obj.stimFile = p.Results.stimFile;
        obj.respFile = p.Results.respFile;
        obj.filterFile = p.Results.filterFile;
        
        obj.windowSize = p.Results.windowSize;
        obj.shiftTime = p.Results.shiftTime;
        obj.percentSV = p.Results.percentSV;
       
    end
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = build(obj, varargin);    
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = buildLandolt(obj, varargin);    
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = buildPrima(obj, varargin);
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = buildPrimaLandolt(obj, varargin);
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = buildHallway(obj, varargin);
    
    % Declare the method for building the training data set    
    [mosaicFile, saveFile] = buildPrimaHallway(obj, varargin);
    
    % Declare the method for learning the filters from the training set    
    [filterFile] = train(obj, varargin);
    
    % Declare the method for testing the accuracy of the recon filters
    [testAccuracy] = test(obj, varargin); 
    
    obj = testCV(obj,stim,spikeResp,varargin);
        
    obj = testImagenet(obj,varargin);
    
    % Declare the method for plotting properties from the recon object
    h = plot(obj, varargin);
    
    % Declare the movie method
    obj = movie(obj);
    
    % Declare the publish method
    obj = publish(obj);
    
%     function window(obj)
%         obj.figureHandle = reconWindow(obj);
%         % Tip: Retrieve guidata using
%         %    gui = guidata(obj.figureHandle);
%         %
%     end
    
end


properties (Constant)
end

% Methods that must only be implemented (Abstract in parent class).

 % Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
end

end