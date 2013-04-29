To have new code actually run on the server, you have to copy it to the proper directory and then re-launch Rserve. To do this from the datelife directory with R, it's

cp *.R /var/FastRWeb/web.R/

then kill Rserve or restart the server

then 

/bin/sh /var/FastRWeb/code/start > /Users/bomeara/Dropbox/recentFastRWebStart.Rout

The command above is automatically run by a cron job two minutes after rebooting (delay so that the server has time to get on the network)
