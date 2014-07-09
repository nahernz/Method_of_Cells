function [mat,arch,load,out] = MOC_read(inputfile)

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
                
        elseif ~isempty(strfind(s,'*MATLAB'))
            out = output(fid);
            
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
        mat{matnum}.V23 = fscanf(fid,'NUA=%f\n',1);
        mat{matnum}.V12 = fscanf(fid,'NUT=%f\n',1);
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
    
    if arch.amod == 1
        % 4 cell square
        vf = fscanf(fid,'VF=%f\n',1);
        df = fscanf(fid,'DF=%e\n',1);
        
        arch.h(1) = sqrt(pi/4*df^2);
        arch.h(2) = arch.h(1)/sqrt(vf)-arch.h(1);
        arch.l = arch.h;
        
        arch.sm = [1 2; 2 2];
        
    elseif arch.amod == 2
        % square symmetric element

         
        vf = fscanf(fid,'VF=%f\n',1);
        df = fscanf(fid,'DF=%e\n',1);
        
        hn = [0, 0.1229, 0.1229, 0.4916, 0.1229, 0.1229, 0]; 
        % fiber with unit area, scale by fiber area 
        H  = (pi/4*df^2)*hn;
        L = H;
        
        SM = [2,2,2,2,2,2,2;
              2,2,2,1,2,2,2;
              2,2,1,1,1,2,2;
              2,1,1,1,1,1,2;
              2,2,1,1,1,2,2;
              2,2,2,1,2,2,2;
              2,2,2,2,2,2,2];  

        fiber  = 1;
        matrix = 2;

        % initialize the volumes
        Vfiber  = 0;
        Vmatrix = 0;

        % calculate the inner square length
        dx = sum(L);
        dy = sum(H);

        % loop over elements
        for x = 1:size(L,2);
            for y = 1:size(H,2);
                if     SM(y,x) == fiber % calculate the volume of the fiber
                    Vfiber = Vfiber + L(x)*H(y);
                elseif SM(y,x) == matrix % calculate the volume of the surrounding matrix
                    Vmatrix = Vmatrix + L(x)*H(y);
                end
            end
        end
        Vfiber/(Vfiber+Vmatrix);
        



            % Calculate the needed extra matrix volume            
            Vm_Needed = Vfiber*(1-vf)/vf-Vmatrix;

            % calculate the needed height and lengths of the surrounding matrix
            
            T = sqrt(dx^2+Vm_Needed)/2 - dx/2; % the thickness of the surrounding membrane
            % Complete the H and L vectors
            H(1) = T;
            H(end) = T;
            L(1) = T;
            L(end) = T;
            
        arch.h = H;
        arch.l = L;
        arch.sm = SM;
        
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
        

        fscanf(fid,'SM=',1);
        for i = 1:dim(1)
            for j = 1:dim(2)
                arch.sm(i,j) = fscanf(fid,'%i,',1);
            end            
        end
        arch.sm
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

function matlab = output(fid)

    matlab.ns = fscanf(fid,'NSE=%i\n',1);  
    matlab.s(1) = fscanf(fid,'SE=%i,',1);
    for i = 2:matlab.ns
        matlab.s(i) = fscanf(fid,'%i,',1);
    end
    fgets(fid);
    
    matlab.ne = fscanf(fid,'NSA=%i\n',1);   
    matlab.e(1) = fscanf(fid,'SA=%i,',1);
    for i = 2:matlab.ne
        matlab.e(i) = fscanf(fid,'%i,',1);
    end
    fgets(fid);
    
    matlab.nE = fscanf(fid,'NSF=%i\n',1);   
    matlab.E(1) = fscanf(fid,'SF=%i,',1);
    for i = 2:matlab.nE
        matlab.E(i) = fscanf(fid,'%i,',1);
    end
    fgets(fid);
    
end

