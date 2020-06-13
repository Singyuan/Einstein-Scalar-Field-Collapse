function logtxt()
t = datetime('now');
nowday = [datestr(t) '\r\n\r\n'];
load('para')
fid=fopen('mylog.txt','a');
fprintf(fid, nowday);
fprintf(fid, 'Initial Condition Prameter\r\n');
fprintf(fid, 'Amp:\t%6.4f\r\n', phi0);
fprintf(fid, 'Power:\t%d\r\n\r\n', pow);
fprintf(fid, 'Spacetime Parameter\r\n');
fprintf(fid, 'Sapce:\t%4.2f\r\n', r(end));
fprintf(fid, 'Total Time:\t%4.2f\r\n', T);
fprintf(fid, 'fg:\t%4.2f\r\n', fg);
fclose(fid);