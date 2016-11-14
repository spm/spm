function sendmail = spm_cfg_sendmail
% SPM Configuration file for sendmail
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_sendmail.m 6929 2016-11-14 13:07:31Z guillaume $

% ---------------------------------------------------------------------
% Recipient
% ---------------------------------------------------------------------
recipient         = cfg_entry;
recipient.tag     = 'recipient';
recipient.name    = 'Recipient';
recipient.help    = {'User to receive mail.'};
recipient.strtype = 's';
recipient.num     = [1 Inf];
% ---------------------------------------------------------------------
% Subject
% ---------------------------------------------------------------------
subject         = cfg_entry;
subject.tag     = 'subject';
subject.name    = 'Subject';
subject.val     = {['[SPM] [%DATE%] On behalf of ' spm('Ver')]};
subject.help    = {'The subject line of the message. %DATE% will be replaced by a string containing the time and date when the email is sent.'};
subject.strtype = 's';
subject.num     = [1 Inf];
% ---------------------------------------------------------------------
% Message
% ---------------------------------------------------------------------
message         = cfg_entry;
message.tag     = 'message';
message.name    = 'Message';
message.val     = {'Hello from SPM!'};
message.help    = {'A string containing the message to send.'};
message.strtype = 's';
message.num     = [1 Inf];
% ---------------------------------------------------------------------
% Attachments
% ---------------------------------------------------------------------
attachments         = cfg_files;
attachments.tag     = 'attachments';
attachments.name    = 'Attachments';
attachments.val     = {{}};
attachments.help    = {'List of files to attach and send along with the message.'};
attachments.filter  = '.*';
attachments.ufilter = '.*';
attachments.num     = [0 Inf];
% ---------------------------------------------------------------------
% SMTP Server
% ---------------------------------------------------------------------
smtp         = cfg_entry;
smtp.tag     = 'smtp';
smtp.name    = 'SMTP Server';
smtp.help    = {'Your SMTP server. If not specified, look for sendmail help.'};
smtp.strtype = 's';
try
    smtp.val = {getpref('Internet','SMTP_Server')};
end
smtp.num     = [1 Inf];
% ---------------------------------------------------------------------
% E-mail
% ---------------------------------------------------------------------
email         = cfg_entry;
email.tag     = 'email';
email.name    = 'E-mail';
email.help    = {'Your e-mail address. Look in sendmail help how to store it.'};
email.strtype = 's';
try
    email.val = {getpref('Internet','E_mail')};
end
email.num     = [1 Inf];
% ---------------------------------------------------------------------
% Zip attachments
% ---------------------------------------------------------------------
zip         = cfg_menu;
zip.tag     = 'zip';
zip.name    = 'Zip attachments';
zip.val     = {'No'};
zip.help    = {'Zip attachments before being sent along with the message.'};
zip.labels  = {'Yes' 'No'}';
zip.values  = {'Yes' 'No'}';
% ---------------------------------------------------------------------
% Parameters
% ---------------------------------------------------------------------
params         = cfg_branch;
params.tag     = 'params';
params.name    = 'Parameters';
params.val     = { smtp email zip};
params.help    = {'Preferences for your e-mail server (Internet SMTP server) and your e-mail address. MATLAB tries to read the SMTP mail server from your system registry. This should work flawlessly. If you encounter any error, identify the outgoing mail server for your electronic mail application, which is usually listed in the application''s preferences, or, consult your e-mail system administrator, and update the parameters. Note that this function does not support e-mail servers that require authentication.'};
% ---------------------------------------------------------------------
% Sendmail
% ---------------------------------------------------------------------
sendmail       = cfg_exbranch;
sendmail.tag   = 'sendmail';
sendmail.name  = 'Sendmail';
sendmail.val   = { recipient subject message attachments params};
sendmail.help  = {'Send a mail message (attachments optionals) to an address.'};
sendmail.prog  = @spm_sendmail;
%_______________________________________________________________________

%_______________________________________________________________________
function spm_sendmail(job)

try
    setpref('Internet','SMTP_Server',job.params.smtp);
    setpref('Internet','E_mail',job.params.email);

    subj = strrep(job.subject,'%DATE%',datestr(now));
    mesg = strrep(job.message,'%DATE%',datestr(now));
    mesg = [mesg 10 10 '-- ' 10 10 'Statistical Parametric Mapping']; 

    if ~isempty(job.attachments)
        if strcmpi(job.params.zip,'Yes')
            zipfile = fullfile(tempdir,'spm_sendmail.zip');
            zip(zipfile,job.attachments);
            job.attachments = {zipfile};
        end
        sendmail(job.recipient,subj,mesg,job.attachments);
    else
        sendmail(job.recipient,subj,mesg);
    end
catch
    %- not an error to prevent an analysis to crash because of just that...
    fprintf('Sendmail failed...\n');
end
