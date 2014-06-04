function [mat,arch,load] = MOC_read(inputfile)

    fid = fopen(inputfile,'r');

    s = fgetl(fid);
    while ~feof(fid)
        if ~isempty(strfind(s,'%')) % comment line
            s = fgetl(fid);
            continue
       
        elseif ~isempty(strfind(s,'*CONSTITUENTS'))
            mat = constituents(fid);
               
        elseif ~isempty(strfind(s,'*CELL'))
            arch = architecture(fid);
                
        elseif ~isempty(strfind(s,'*LOADING'))
            load = loading(fid);
                
        elseif ~isempty(strfind(s,'*END'))
            return
        
        else
            fprintf(2,'LINE IGNORED!\n\n\t%s',s);
        end
        
        s = fgetl(fid);
    end
end

function mat = constituents(fid)

    nmats = fscanf(fid,'NMATS=%i\n',1);
    mat = cell(nmats,1); % size the material database
    
    matnum = fscanf(fid,'MAT=%i,CMOD=%i\n',2);
    cmod = matnum(2); matnum(2) = [];
    
    mat{matnum}.cmod = cmod;
    if cmod == 1
        mat{matnum}.E11 = fscanf(fid,'EA=%e\n',1);
        mat{matnum}.E22 = fscanf(fid,'ET=%e\n',1);
        mat{matnum}.V11 = fscanf(fid,'NUA=%f\n',1);
        mat{matnum}.V22 = fscanf(fid,'NUT=%f\n',1);
        mat{matnum}.G23 = fscanf(fid,'GA=%e\n',1);
        mat{matnum}.G12 = fscanf(fid,'GT=%e\n',1);
    elseif cmod == 2
        mat{matnum}.E = fscanf(fid,'E=%e\n',1);
        mat{matnum}.V = fscanf(fid,'NU=%f\n',1);
    end
    s = fgetl(fid);
    
    clear matnum
    matnum = fscanf(fid,'MAT=%i,CMOD=%i\n',2);
    cmod = matnum(2); matnum(2) = [];

    mat{matnum}.cmod = cmod;
    if cmod == 1
        mat{matnum}.E11 = fscanf(fid,'EA=%e\n',1);
        mat{matnum}.E22 = fscanf(fid,'ET=%e\n',1);
        mat{matnum}.V11 = fscanf(fid,'NUA=%f\n',1);
        mat{matnum}.V22 = fscanf(fid,'NUT=%f\n',1);
        mat{matnum}.G23 = fscanf(fid,'GA=%e\n',1);
        mat{matnum}.G12 = fscanf(fid,'GT=%e\n',1);
    elseif cmod == 2
        mat{matnum}.E = fscanf(fid,'E=%e\n',1);
        mat{matnum}.V = fscanf(fid,'NU=%f\n',1);
    end
   
end

function arch = architecture(fid)

    arch.amod = fscanf(fid,'AMOD=%i\n',1);
    
    if ismember(arch.amod,[1 2 3])
        arch.vf = fscanf(fid,'VF=%f\n',1);
        arch.df = fscanf(fid,'DF=%e\n',1);
        
    elseif arch.amod == 4
        dim = fscanf(fid,'DIM=%i,%i,\n',2);
        arch.h(1) = fscanf(fid,'H=%e,',1);
        for i = 2:dim(1)
            arch.h(i) = fscanf(fid,'%e,',1);
        end
        fgets(fid);
        arch.l(1) = fscanf(fid,'L=%e,',1);
        for i = 2:dim(2)
            arch.l(i) = fscanf(fid,'%e,',1);
        end
        fgets(fid);
        arch.sm(1) = fscanf(fid,'SM=%e,',1);
        for i = 2:dim(1)*dim(2)
            arch.sm(i) = fscanf(fid,'%i,',1);
        end
        fgets(fid);
        
    end
end

function load = loading(fid)

    load.lmod = fscanf(fid,'LMOD=%i\n',1);
    
    if ismember(load.lmod,[1 2 3])
        load.nl = fscanf(fid,'NL=%i\n',1);
        
        load.l(1) = fscanf(fid,'L=%f,',1);
        for i = 2:load.nl
            load.l(i) = fscanf(fid,'%f,',1);
        end
        fgets(fid);
    else
        % add general loading here
    end
end

