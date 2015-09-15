function meshstruc = BCsetting(meshstruc,info)

meshstruc = feval(info.BC,meshstruc,info);