function c0 = cfg_mlbatch_root
% function c = cfg_mlbatch_root
% The root node of a matlabbatch configuration. This file is called by
% cfg_util when initialising its internal data structure.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_mlbatch_root.m 2308 2008-10-06 16:05:26Z volkmar $

rev = '$Rev: 2308 $'; %#ok
c0      = cfg_repeat;
c0.name = 'MATLAB Batch System';
c0.tag  = 'matlabbatch';
c0.num  = [1 Inf];
c0.forcestruct = true;
c0.help = {
                '%Matlabbatch User Interface'
                ''
                '%* Menu and Toolbar'
                '/*\subsection*{Menu and Toolbar}*/'
                'The "File" and "Edit" menu offer options to load, save and run a job and to modify the configuration of the batch system. For each application which is known to the batch system, a separate pulldown menu lists the available modules. Depending on the application, these modules may be grouped into submenus. Application specific defaults can be edited by choosing "Edit Defaults" from the application menu. The toolbar offers some shortcuts to frequently used operations (e.g. load, save, run).'
                'Jobs are saved as MATLAB .m files. These files contain a MATLAB script, which can be executed in MATLAB to recreate the job variable. Multiple jobs can be loaded at once. This allows to concatenate parts of a job.'
                'Most actions are accessible via toolbar or buttons in the main figure window. Two important exceptions are:'
                '/*\begin{description}\item[Run Batch (serial)]*/'
                '%* Run Batch (serial)'
                'If a batch fails because modules depend on each other, but default execution order is wrong, then the batch can be run module by module instead.'
                '/*\item[Clear Value]*/'
                '%* Clear Value'
                'In the value editing dialog boxes, it is not possible to clear a value without entering a replacement. To unset an input, "Clear Value" from the Edit Menu has to be used.'
                '/*\end{description}*/'
                ''
                '%* Top Left Panel'
                '/*\subsection*{Top Left Panel}*/'
                'The current job, which is represented as a list of executable modules. Modules marked with DEP depend on the successful execution of other modules in the job. Modules marked with X still require some values to be set before the job can be run, although an incompletely specified job can still be saved and loaded.'
                ''
                '%* Top Right Panel'
                '/*\subsection*{Top Right Panel}*/'
                'These are the configuration details for the currently selected module. Items marked with DEP depend on the successful execution of other modules in the job. Items marked with X still require some values to be set before the job can be run. Depending on the kind of detail currently selected, a choice of buttons appears below the Centre Right Panel to manipulate the current value.'
                ''
                '%* Centre Right Panel'
                '/*\subsection*{Centre Right Panel}*/'
                'This panel shows the current value of the highlighted item (where relevant).'
                ''
                '%* Edit Buttons'
                '/*\subsection*{Edit Buttons}*/'
                'Depending on the type of configuration item, different edit buttons appear.'
                '/*\begin{description}*/'
                '/*\item[Files]*/'
                '%* Files'
                '"Select Files" opens a file selection dialog box to select multiple files. "Edit Value" opens a generic value edit dialog to edit the list of files. "Dependencies" offers a list of outputs from other modules that can be used as an input to this item.'
                '/*\item[Generic Value]*/'
                '%* Generic Value'
                '"Edit Value" opens a generic value edit dialog to edit the list of files. "Dependencies" offers a list of outputs from other modules that can be used as an input to this item.'
                '/*\item[Menu]*/'
                '%* Menu'
                '"Edit Value" opens a selection dialog showing allowed menu options.'
                '/*\item[Choice]*/'
                '%* Choice'
                '"Edit Value" opens a selection dialog showing allowed menu options. Depending on the choosen option the module configuration may change.'
                '/*\item[Repeat]*/'
                '%* Repeat'
                '"Add Item", "Replicate Item", "Delete Item" allow to add new repeated items, to replicate or to delete items from the list. If more than one item or item type exists, a dialog popup will appear listing the available options. Multiple selections are allowed.'
                '/*\end{description}*/'
                ''
                '%* Bottom Panel'
                '/*\subsection*{Bottom Panel}*/'
                'This panel provides information about the meaning of the current item.'
                '/*\begin{figure} \begin{center} \includegraphics[width=70mm]{images/batch_ui1} \includegraphics[width=70mm]{images/batch_ui2} \includegraphics[width=70mm]{images/ui3} \includegraphics[width=70mm]{images/ui4}\end{center} \caption{The SPM5 user interface. \emph{Top left:} The usual user-interface.  \emph{Top right:} The Defaults user-interface. \emph{Bottom left:} The file selector (click the (?) button for more information about filtering filenames, or selecting individual volumes within a 4D file). \emph{Bottom right:} more online help can be obtained via the main help button.} \end{figure} */'
    };