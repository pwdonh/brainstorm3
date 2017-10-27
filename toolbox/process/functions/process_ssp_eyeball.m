function varargout = process_ssp_eyeball( varargin )
% PROCESS_SSP2_EOG: Reject cardiac artifact for a group of recordings file (calculates SSP from FilesA and applies them to FilesB)
%
% USAGE:  OutputFiles = process_ssp2_eog('Run', sProcess, sInputs)

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, Elizabeth Bock, 2011-2016

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SSP: Eyeball';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Artifacts';
    sProcess.Index       = 2000;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsSsp?highlight=%28Process2%29#Troubleshooting';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
    % Time point
    sProcess.options.coord1.Comment = 'SCS coordinates 1st eye [x,y,z]:  ';
    sProcess.options.coord1.Type    = 'value';
    sProcess.options.coord1.Value   = {[76 -30 0], 'mm', 1};
    sProcess.options.coord2.Comment = 'SCS coordinates 2nd eye [x,y,z]:  ';
    sProcess.options.coord2.Type    = 'value';
    sProcess.options.coord2.Value   = {[76 30 0], 'mm', 1};    
    % === Number of scan and project iterations
    sProcess.options.ncomp.Comment = 'Number of components to select';
    sProcess.options.ncomp.Type    = 'value';
    sProcess.options.ncomp.Value   = {2, '', 0};
    % Split up
    sProcess.options.fixeyes.Comment = 'Fix the two eyes..';
    sProcess.options.fixeyes.Type    = 'checkbox';
    sProcess.options.fixeyes.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    if isfield(sProcess.options, 'eventname') && ~isempty(sProcess.options.eventname.Value)
        Comment = ['SSP EOG: ' sProcess.options.eventname.Value];
    else
        Comment = sProcess.Comment;
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    % Input arguments
    coords(:,1) = sProcess.options.coord1.Value{1}./1000;
    coords(:,2) = sProcess.options.coord2.Value{1}./1000;
    nSelect = sProcess.options.ncomp.Value{1};
    isFixEyes = sProcess.options.fixeyes.Value;
    
    % Load some files
    SubjectMat = bst_get('Subject', sInputs.SubjectFile);
    CortexFile = SubjectMat.Surface(SubjectMat.iCortex).FileName;
    ChannelMat = in_bst_channel(sInputs.ChannelFile);
    iMeg = channel_find(ChannelMat.Channel,'MEG');
    Channels = ChannelMat.Channel(iMeg);
    
    % Compute gain matrix
    sSurfInner = tess_envelope(file_fullpath(CortexFile), 'convhull', 1082, .003);
    sph = bst_os(Channels, double(sSurfInner.Vertices), double(sSurfInner.Faces));
    G = bst_meg_sph(coords, Channels,  sph);
    
    % SVD of gain matrix
    if isFixEyes
        GG(:,1) = G(:,1) + G(:,4);
        GG(:,2) = G(:,2) + G(:,5);
        GG(:,3) = G(:,3) + G(:,6);
        nComp = 3;
    else
        GG = G;
        nComp = 6;
    end
    [U, S] = svd(GG, 'econ');
    s = diag(S)';
    
    % Add to channelfile
    nProj = numel(ChannelMat.Projector);
    ChannelMat.Projector(nProj+1).Comment = 'SSP: Eyeball';
    ChannelMat.Projector(nProj+1).Components = zeros(numel(ChannelMat.Channel),nComp);
    ChannelMat.Projector(nProj+1).Components(iMeg,:) = U;
    ChannelMat.Projector(nProj+1).CompMask = zeros(1,nComp);
    ChannelMat.Projector(nProj+1).CompMask(1:nSelect) = 1;
    ChannelMat.Projector(nProj+1).Status = 1;
    ChannelMat.Projector(nProj+1).SingVal = s;
    
    % Save channelfile
    bst_save(file_fullpath(sInputs.ChannelFile), ChannelMat);
    
    OutputFiles = {sInputs.FileName};
    
end



