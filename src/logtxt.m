function logtxt()
t = datetime('now');
nowday = [datestr(t) '\r\n\r\n'];
load('para')
fid=fopen('../result/data/mylog.txt','a');
fprintf(fid, nowday);
fprintf(fid, 'Initial Condition Prameter\r\n');
fprintf(fid, 'Amp:\t%6.4f\r\n', phi0);
fprintf(fid, 'Spacetime Parameter\r\n');
fprintf(fid, 'Sapce boundary:\t%4.2f\r\n', rend);
fprintf(fid, 'Total Time:\t%4.2f\r\n', T);
fprintf(fid, 'fg:\t%4.2f\r\n', fg);
fprintf(fid, '\r\n');
fclose(fid);