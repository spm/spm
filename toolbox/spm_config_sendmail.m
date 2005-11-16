function c = spm_config_sendmail(varargin)
% Configuration file for sending emails
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Guillaume Flandin
% $Id: spm_config_sendmail.m 302 2005-11-16 21:42:45Z guillaume $

to.type = 'entry';
to.name = 'Recipient';
to.tag  = 'recipient';
to.strtype = 's';
to.num  = [1 1];
to.help = {'User to receive mail.'};

subject.type = 'entry';
subject.name = 'Subject';
subject.tag  = 'subject';
subject.strtype = 's';
subject.num  = [1 1];
subject.val  = {['[SPM] [%DATE%] On behalf of ' spm('Ver')]};
subject.help = {['The subject line of the message. %DATE% will be ',...
'replaced by a string containing the time and date when the email is sent.']};

message.type = 'entry';
message.name = 'Message';
message.tag  = 'message';
message.strtype = 's';
message.num  = [1 1];
message.val  = {'Hello from SPM!'};
message.help = {'A string containing the message to send.'};

attachments.type = 'files';
attachments.name = 'Attachments';
attachments.tag  = 'attachments';
attachments.filter = '.*';
attachments.num  = [0 Inf];
attachments.val  = {{}};
attachments.help = {'List of files to attach and send along with the message.'};

smtpserver.type = 'entry';
smtpserver.name = 'SMTP Server';
smtpserver.tag  = 'smtp';
smtpserver.strtype = 's';
smtpserver.num  = [1 1];
try
	smtpserver.val = {getpref('Internet','SMTP_Server')};
end
smtpserver.help = {'Your SMTP server. If not specified, look for sendmail help.'};

email.type = 'entry';
email.name = 'E-mail';
email.tag  = 'email';
email.strtype = 's';
email.num  = [1 1];
try
	email.val = {getpref('Internet','E_mail')};
end
email.help = {'Your e-mail address. Look in sendmail help how to store it.'};

zipattach.type = 'menu';
zipattach.name = 'Zip attachments';
zipattach.tag  = 'zip';
zipattach.labels = {'Yes' 'No'};
zipattach.values = {'Yes' 'No'};
zipattach.val  = {'No'};
zipattach.help = {'Zip attachments before being sent along with the message.'};

parameters.type = 'branch';
parameters.name = 'Parameters';
parameters.tag  = 'params';
parameters.val  = {smtpserver,email,zipattach};
parameters.help = {['Preferences for your e-mail server (Internet SMTP ',...
'server) and your e-mail address. MATLAB tries to read the SMTP mail ',...
'server from your system registry. This should work flawlessly. If you '...
'encounter any error, identify the outgoing mail server for your ',...
'electronic mail application, which is usually listed in the ',...
'application''s preferences, or, consult your e-mail system administrator, ',...
'and update the parameters. Note that this function does not support ',... 
'e-mail servers that require authentication.']};

c.type = 'branch';
c.name = 'Sendmail';
c.tag  = 'sendmail';
c.val  = {to,subject,message,attachments,parameters};
c.prog = @spm_sendmail;
c.help = {['Send a mail message (attachments optionals) to an ',...
'address.']};

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
