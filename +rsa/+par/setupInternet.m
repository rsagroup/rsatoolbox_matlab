% This function initialises internet settings 
%FJ Feb-2014
function seupInternet()
name='RSA_Toolbox';
mail = 'rsatoolbox@gmail.com'; %Your GMail email address
password = 'Neurolex'; %Your GMail password


setpref('Internet','E_mail',name);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
setpref('Internet', 'From', name);
end