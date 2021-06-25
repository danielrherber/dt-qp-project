%--------------------------------------------------------------------------
% DTQP_commandWindowTasks.m
% Common command window tasks
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
switch flag
%--------------------------------------------------------------------------
case 'line'
% disp('_______________________________________________________________')
% disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% disp('***************************************************************')
% disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
% disp('---------------------------------------------------------------')
disp(char(repmat(9644,1,63)))
%--------------------------------------------------------------------------
case 'banner'
disp('                         __ ___ __  __                         ')
disp('                        |  \ | /  \|__)                        ')
disp('                        |__/ | \_\/|                           ')
disp('                                                               ')
% disp("                        █▀▄ ▀█▀ █▀█ █▀█                         ")
% disp("                        █▄▀  █  ▀▀█ █▀▀                         ")
%--------------------------------------------------------------------------
case 'link'
disp('Link: <a href = "https://github.com/danielrherber/dt-qp-project">https://github.com/danielrherber/dt-qp-project</a>')
%--------------------------------------------------------------------------
case 'info'
str = strings(0);
str(1) = strcat("Machine: ",string(java.net.InetAddress.getLocalHost.getHostName));
str(2) = strcat("Date: ",datestr(now,'yyyy/mm/dd HH:MM'));
str = strjoin(str,' | ');
disp(str)
%--------------------------------------------------------------------------
case 'contributor'
disp('<strong>Primary contributor</strong>: Daniel Herber (danielrherber on GitHub)')
%--------------------------------------------------------------------------
case 'creation-time'
disp(strcat(string(char(9658))," Creation time: ",string(ctime)," s"))
%--------------------------------------------------------------------------
case 'solver-time'
disp(strcat(string(char(9658))," Solver time: ",string(opts.timer.qpsolver)," s"))
%--------------------------------------------------------------------------
case 'total-time'
disp(strcat(string(char(9658))," Total time: ",string(opts.timer.total)," s"))
%--------------------------------------------------------------------------
end