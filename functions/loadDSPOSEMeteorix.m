function rso = loadDSPOSEMeteorix()
    [fid, errormsg] = fopen('inputfiles/dspose_meteorix.txt', 'r');
    if fid<0
        disp(errormsg);
        return;
    end
    for i = 1:5; fgets(fid); end % ignore first 5 lines
    
    A = fscanf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[15 inf]);

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
    
    facets_panel(1:16) = clTriangle();
    for i = 1:16
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
    rso.name = 'D-SPOSE Meteorix';
    
    tar_input = clCompoundTargetInput();
    
    inputstr(1:5) = clCompoundSegmentsInput();
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
    inputstr(1).name = 'DSPOSE Meteorix Body';
    inputstr(1).density = 3.980832860420303e+02;
    inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
    inputstr(1).resolution = 5;
    inputstr(1).solid = true;
    inputstr(1).twosided = false;
    inputstr(1).thickness = 0.003;
    inputstr(1).segobj = seg_body;
    
    %% PANEL
    for i = 1:4
        seg_panel = clCompoundSegment();
        i_s = i*2-1;
        i_e = i*2;
        seg_panel.facets = facets_panel(i_s:i_e);
        seg_panel.segmentReferenceCube = getBoundingBoxfromFacets(seg_panel.facets);
      
        inputstr(i+1).bSTL = false;
        inputstr(i+1).stlFile = "";
        inputstr(i+1).surfAtr = clSA_wilken2015();
        a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
        t = [15.04 1.19 2.41 23.38];
        inputstr(i+1).surfAtr.init(a, t);
        inputstr(i+1).surfAtr.name = 'wilken2015_aluminium';
        inputstr(i+1).offset = [0; 0; 0];
        inputstr(i+1).name = 'DSPOSE Meteorix Panel';
        inputstr(i+1).density = 1.933992372478913e+03;
        inputstr(i+1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
        inputstr(i+1).resolution = 5;
        inputstr(i+1).solid = false;
        inputstr(i+1).twosided = false;
        inputstr(i+1).thickness = 0.003;
        inputstr(i+1).segobj = seg_panel;
    end
    
    %% INIT
    tar_input.bUseProvidedParam = true;
    tar_input.inertia = zeros(3);
    tar_input.inertia(1,1) = 0.0586;
    tar_input.inertia(2,2) = 0.0589;
    tar_input.inertia(3,3) = 0.0432;
    tar_input.mass = 3.78375;
    tar_input.barycenter = [0; 0; 0];
    rso.init(inputstr, tar_input);
    
end
    
    