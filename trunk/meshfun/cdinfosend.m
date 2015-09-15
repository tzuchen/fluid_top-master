function [] = cdinfosend(socketp,cdinfo)
% This function sends the structure cdinfo to Petsc program.
% In Petsc side, it has to be received by function CdinfoReceive
%

coordarray = cell2mat(cdinfo.coord);


send(socketp,cdinfo.dim);
send(socketp,[coordarray.l]);
send(socketp,[coordarray.r]);
send(socketp,[coordarray.n]);
send(socketp,[coordarray.gsize]);
send(socketp,[coordarray.p]);
send(socketp,cell2mat(cdinfo.gridspam));
send(socketp,cell2mat(cdinfo.nofgrid));







