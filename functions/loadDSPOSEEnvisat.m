function rso = loadDSPOSEEnvisat()
    [fid, errormsg] = fopen('inputfiles/dspose_envisat.txt', 'r');
    if fid<0
        disp(errormsg);
        return;
    end
    for i = 1:5; fgets(fid); end % ignore first 5 lines
    
    A = fscanf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[18 inf]);

    facets_body(1:12) = clTriangle();
    for i = 1:12
        facets_body(i) = clTriangle();
        facets_body(i).normal = -A(1:3,i);
        facets_body(i).base = A(4:6,i);
        facets_body(i).edge1 = A(7:9,i) - A(4:6,i);
        facets_body(i).edge2 = A(10:12,i) - A(4:6,i);
        
        facets_body(i).init();
    end
    facets_body = simplifyFacets(facets_body, 0.001, 0.01);
    
    facets_panel(1:4) = clTriangle();
    for i = 1:4
        j=i+12;
        facets_panel(i) = clTriangle();
        facets_panel(i).normal = -A(1:3,j);
        facets_panel(i).base = A(4:6,j);
        facets_panel(i).edge1 = A(7:9,j) - A(4:6,j);
        facets_panel(i).edge2 = A(10:12,j) - A(4:6,j);
        
        facets_panel(i).init();
    end
    facets_panel = simplifyFacets(facets_panel, 0.001, 0.01);
    
    rso = clCompoundTarget();
    rso.name = 'EnviSat';
    
    tar_input = clCompoundTargetInput();
    
    inputstr(1:2) = clCompoundSegmentsInput();
    %% BODY
    seg_body = clCompoundSegment();
    seg_body.facets = facets_body;
    seg_body.segmentReferenceCube = getBoundingBoxfromFacets(facets_body);
    
    inputstr(1).bSTL = false;
    inputstr(1).stlFile = "";
    inputstr(1).surfAtr = clSA_wilken2015();
    a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
    t = [15.04 1.19 2.41 23.38];
    inputstr(1).surfAtr.init(a, t);
    inputstr(1).surfAtr.name = 'wilken2015_aluminium';
    inputstr(1).offset = [0; 0; 0];
    inputstr(1).name = 'DSPOSE EnviSat Body';
    inputstr(1).density = 3.980832860420303e+02;
    inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
    inputstr(1).resolution = 5;
    inputstr(1).solid = true;
    inputstr(1).twosided = false;
    inputstr(1).thickness = 0.003;
    inputstr(1).segobj = seg_body;
    
    %% PANEL
    seg_panel = clCompoundSegment();
    seg_panel.facets = facets_panel;
    seg_panel.segmentReferenceCube = getBoundingBoxfromFacets(facets_panel);
    
    inputstr(2).bSTL = false;
    inputstr(2).stlFile = "";
    inputstr(2).surfAtr = clSA_wilken2015();
    a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
    t = [15.04 1.19 2.41 23.38];
    inputstr(2).surfAtr.init(a, t);
    inputstr(2).surfAtr.name = 'wilken2015_aluminium';
    inputstr(2).offset = [0; 0; 0];
    inputstr(2).name = 'DSPOSE EnviSat Panel';
    inputstr(2).density = 1.933992372478913e+03;
    inputstr(2).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
    inputstr(2).resolution = 5;
    inputstr(2).solid = false;
    inputstr(2).twosided = false;
    inputstr(2).thickness = 0.003;
    inputstr(2).segobj = seg_panel;
    
    %% INIT
    tar_input.bUseProvidedParam = true;
    tar_input.inertia = [17023.3 397.1 -2171.4; 397.1 124825.7 344.2; -2171.4 344.2 129112.2];
    tar_input.mass = 7827.867;
    tar_input.barycenter = [0; 0; 0];
    tar_input.magnetic_tensor = diag([931500,1059000 1059000]);    
    rso.init(inputstr, tar_input);
    
end
    
    