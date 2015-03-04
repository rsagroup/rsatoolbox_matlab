%FJ 02/2014

function setupEmail(mailto)

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
msg =['Your experiment was completed at ', datestr(clock)];
% Send the email. Note that the first input is the address you are sending the email to
sendmail(mailto,'Your results are ready now!', msg);

end
