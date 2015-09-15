function [] = cdinfosend2file(filename,cdinfo)

coordarray = cell2mat(cdinfo.coord);
allinfo = [cdinfo.dim ,[coordarray.l], [coordarray.r],[coordarray.n],[coordarray.p],cell2mat(cdinfo.gridspam),cell2mat(cdinfo.nofgrid)];
send2file(filename,allinfo);
