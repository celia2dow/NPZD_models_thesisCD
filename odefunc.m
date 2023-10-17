function dxdt = odefunc(t,x,params)
% Edited ODE system for NPZD model
    % We require that solutions are ordered: long1 lat1, long1 lat2, ..., 
    % long1 latEnd, long2 lat1, long2 lat2, ..., long2 latEnd, ...
    lattice_size = [params.num_lats params.num_longs];

    % Convert state variables from vector-form to matrix-form
    % If the site is a landmass, maintain 0 concentration and 0 gradient
    % N
    Nvec = x(1 : params.num_sites);
    N = reshape(Nvec, lattice_size);
    N(logical(params.landmass)) = 0;
    % P
    Pvec = x(params.num_sites+1 : 2*params.num_sites);
    P = reshape(Pvec, lattice_size);
    P(logical(params.landmass)) = 0;
    % Z
    Zvec = x(2*params.num_sites+1 : 3*params.num_sites);
    Z = reshape(Zvec, lattice_size);
    Z(logical(params.landmass)) = 0;
    if params.detritus_layer
        % D
        Dvec = x(3*params.num_sites+1 : end);
        D = reshape(Dvec, lattice_size);
        D(logical(params.landmass)) = 0;
    end
    
    % Rename params
    a=params.a; b=params.b; c=params.c; d=params.d; e=params.e; k=params.k; 
    q=params.q; s=params.s; N_0=params.N_0; 
    alpha=params.alpha; beta=params.beta; gamma=params.gamma; 
    lambda=params.lambda; mu=params.mu; r=params.r; omega = params.omega; 
    phi = params.phi; psi=params.psi; delta_x = params.delta_x; 
    delta_y = params.delta_y; dif_N = params.dif_N; dif_P = params.dif_P; 
    dif_Z = params.dif_Z; dif_D = params.dif_D;
    
    dNdt = zeros(params.num_lats, params.num_longs);
    dPdt = zeros(params.num_lats, params.num_longs);
    dZdt = zeros(params.num_lats, params.num_longs);
    dDdt = zeros(params.num_lats, params.num_longs);
    dxdt = zeros(length(x),1);

    % Day within year, needed for repeated data
    day = mod(t,params.num_days);
    if day == 0
        day = params.num_days;
    end

    % Create spatial cubic splines for easterly and northerly currents
    if params.current
        vN_grid = zeros(params.num_lats, params.num_longs);
        vE_grid = zeros(params.num_lats, params.num_longs);
        vE_CS = cell(params.num_lats);          % a cubic spline with respect to j*delta_x
        vN_CS = cell(params.num_longs);         % a cubic spline with respect to i*delta_y
        lats = (1:params.num_lats).*delta_y;
        longs = (1:params.num_longs).*delta_x;
        for i = 1:params.num_lats
            for j = 1:params.num_longs
                if ~params.landmass(i,j)        % If there is no landmass at (i,j)
                    vN_grid(i,j) = CSfunction(params.v_north_CS{i,j},day);
                    vE_grid(i,j) = CSfunction(params.v_east_CS{i,j},day);
                else                            % If there is a landmass at (i,j)
                    vN_grid(i,j) = 0;           % Set the current values to 0
                    vE_grid(i,j) = 0;
                end
            end
            vE_CS{i}= spline(longs,vE_grid(i,:));   % For easterly currents 
        end
        for j = 1:params.num_longs
            vN_CS{j} = spline(lats,vN_grid(:,j)');  % For northerly currents
        end
    end

    for i = 1:params.num_lats
        for j = 1:params.num_longs
            % Mixed-layer depth (LEV), temperature stratification index 
            % (TSI), gradient of depth (h), and gradient of currents 
            % (dvNdt and dvEdt)

            %delta_x = params.measures(1,i,j);
            %delta_y = params.measures(2,i,j);
            if params.temp_vert_mix
                if ~params.landmass(i,j)    % if site (i,j) is not a landmass
                    T = CSfunction(params.temp_CS{i,j},day);
                    TSI = CSfunction(params.TSI_CS{i,j},day);
                    strat = params.strat_func(TSI); % stratification factor 
                end
            end

            if params.change_mix_depth
                if ~params.landmass(i,j)    % if site (i,j) is not a landmass
                    LEV = CSfunction(params.MLD_CS{i,j},day);
                    h = CSgradient(params.MLD_CS{i,j},day);
                else                        % if site (i,j) is a landmass
                    LEV = 1;
                    h = 0;
                end
            end

            if params.current
                dvNdy = CSgradient(vN_CS{j},i.*delta_y);
                dvEdx = CSgradient(vE_CS{i},j.*delta_x);
            end
            
            if params.temp_vert_mix && params.change_mix_depth
                % a, k, s and psi need to be updated with the new depth
                % Switch for temperature-dependent vertical mixing/sinking
                a = params.a*params.assumed_MLD/LEV;
                k = params.k*strat*params.assumed_MLD/LEV;
                s = params.s*strat*params.assumed_MLD/LEV;
                psi = params.psi*strat*params.assumed_MLD/LEV;
            elseif params.change_mix_depth
                % Change with mixed layer depth
                % a, k, s and psi need to be updated with the new depth
                a = params.a*params.assumed_MLD/LEV;
                k = params.k*params.assumed_MLD/LEV;
                s = params.s*params.assumed_MLD/LEV;
                psi= params.psi*params.assumed_MLD/LEV;
            elseif params.temp_vert_mix
                % Switch for temperature-dependent vertical mixing/sinking
                k = params.k*strat;
                s = params.s*strat;
                psi = params.psi*strat;
            end

            % ODE SYSTEM equations 1-3
            if params.paper == 1999
                h_Z = q*Z(i,j);
            elseif params.paper == 1996
                h_Z = d*Z(i,j)^2;
            end

            % Switch for detritus layer inclusion
            if params.detritus_layer % NPZD model
                dNdt(i,j) = -N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) + ...
                    beta*lambda*(P(i,j)^2+omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) + gamma*h_Z + ...
                    phi*D(i,j) + k*(N_0-N(i,j));
                dPdt(i,j) = N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) - ...
                    r*P(i,j) - lambda*P(i,j)^2*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - (s+k)*P(i,j);
                dZdt(i,j) = alpha*lambda*(P(i,j)^2+omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - h_Z;
                dDdt(i,j) = r*P(i,j) + ...
                    lambda*((1-alpha-beta)*P(i,j)^2-(alpha+beta)*omega*D(i,j)^2)*Z(i,j)/...
                    (mu^2+P(i,j)^2+omega*D(i,j)^2) - (phi+psi+k)*D(i,j);
            else % NPZ model
                dNdt(i,j) = -N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) + ...
                    r*P(i,j) + beta*lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) ...
                    + gamma*h_Z + k*(N_0-N(i,j));
                dPdt(i,j) = N(i,j)*P(i,j)*a/((e+N(i,j))*(b+c*P(i,j))) - ...
                    r*P(i,j) - lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) - ...
                    (s+k)*P(i,j);
                dZdt(i,j) = alpha*lambda*P(i,j)^2*Z(i,j)/(mu^2+P(i,j)^2) - ...
                    h_Z;
            end
            
            % Change with mixed layer depth
            if params.change_mix_depth == 1
                dNdt(i,j) = dNdt(i,j) + max(h,0) * (N_0-N(i,j))/LEV;
                dPdt(i,j) = dPdt(i,j) - max(h,0) * P(i,j)/LEV;
                dZdt(i,j) = dZdt(i,j) - h*Z(i,j)/LEV;
                if params.detritus_layer
                    dDdt(i,j) = dDdt(i,j) - max(h,0) * D(i,j)/LEV;
                end
            end

            % Change with current
            if params.current == 1
                % Determine indices of nearest neighbours with periodic
                % boundary conditions
                sd = 1; su = 1; sl = 1; sr = 1;
                if i > 1
                    down = i-1;
                else
                    down = params.num_lats;
                    if params.BC == 1
                        sd = 0;
                    end
                end

                if i < params.num_lats
                    up = i+1;
                else
                    up = 1;
                    if params.BC == 1
                        su= 0;
                    end
                end

                if j > 1 
                    left = j-1;
                else
                    left = params.num_longs;
                    if params.BC == 1
                        sl = 0;
                    end
                end

                if j < params.num_longs
                    right = j+1;
                else
                    right = 1;
                    if params.BC == 1
                        sr = 0;
                    end
                end

                % Use land as a boundary too
                if params.landmass(i,left)
                    left = params.num_longs;
                    if params.BC == 1
                        sl = 0;
                    end
                end
                if params.landmass(i,right)
                    right = 1;
                    if params.BC == 1
                        sr = 0;
                    end
                end
                if params.landmass(down,j)
                    down = params.num_lats;
                    if params.BC == 1
                        sd = 0;
                    end
                end
                if params.landmass(up,j)
                    up = 1;
                    if params.BC == 1
                        su= 0;
                    end
                end
                
                % For species u, the time derivative is:
                %   du/dt = f(u) - u(dvE/dx + dvN/dy) - vE*du/dx - vN*du/dy
                % where f(u) is the previous expression for du/dt
                dNdt(i,j) = dNdt(i,j) - N(i,j)*(dvEdx+dvNdy) - ...
                    (max(vE_grid(i,j),0)*N(i,j) - sl*max(vE_grid(i,left),0)*N(i,left))/delta_x - ...
                    (sr*min(vE_grid(i,right),0)*N(i,right) - min(vE_grid(i,j),0)*N(i,j))/delta_x - ...
                    (max(vN_grid(i,j),0)*N(i,j) - sd*max(vN_grid(down,j),0)*N(down,j))/delta_y - ...
                    (su*min(vN_grid(up,j),0)*N(up,j) - min(vN_grid(i,j),0)*N(i,j))/delta_y;
                dPdt(i,j) = dPdt(i,j) - P(i,j)*(dvEdx+dvNdy) - ...
                    (max(vE_grid(i,j),0)*P(i,j) - sl*max(vE_grid(i,left),0)*P(i,left))/delta_x - ...
                    (sr*min(vE_grid(i,right),0)*P(i,right) - min(vE_grid(i,j),0)*P(i,j))/delta_x - ...
                    (max(vN_grid(i,j),0)*P(i,j) - sd*max(vN_grid(down,j),0)*P(down,j))/delta_y - ...
                    (su*min(vN_grid(up,j),0)*P(up,j) - min(vN_grid(i,j),0)*P(i,j))/delta_y;
                dZdt(i,j) = dZdt(i,j) - Z(i,j)*(dvEdx+dvNdy) - ...
                    (max(vE_grid(i,j),0)*Z(i,j) - sl*max(vE_grid(i,left),0)*Z(i,left))/delta_x - ...
                    (sr*min(vE_grid(i,right),0)*Z(i,right) - min(vE_grid(i,j),0)*Z(i,j))/delta_x - ...
                    (max(vN_grid(i,j),0)*Z(i,j) - sd*max(vN_grid(down,j),0)*Z(down,j))/delta_y - ...
                    (su*min(vN_grid(up,j),0)*Z(up,j) - min(vN_grid(i,j),0)*Z(i,j))/delta_y;
                if params.detritus_layer
                    dDdt(i,j) = dDdt(i,j) - D(i,j)*(dvEdx+dvNdy) - ...
                        (max(vE_grid(i,j),0)*D(i,j) - sl*max(vE_grid(i,left),0)*D(i,left))/delta_x - ...
                        (sr*min(vE_grid(i,right),0)*D(i,right) - min(vE_grid(i,j),0)*D(i,j))/delta_x - ...
                        (max(vN_grid(i,j),0)*D(i,j) - sd*max(vN_grid(down,j),0)*D(down,j))/delta_y - ...
                        (su*min(vN_grid(up,j),0)*D(up,j) - min(vN_grid(i,j),0)*D(i,j))/delta_y;
                end
            end
            
            % Change with diffusion
            if params.diffus
                dNdt(i,j) = dNdt(i,j) +  dif_N * (...
                    (N(down,j)+N(up,j)-2*N(i,j))/(delta_y^2) + ...
                    (N(i,left)+N(i,right)-2*N(i,j))/(delta_x^2));
                dPdt(i,j) = dPdt(i,j) +  dif_P * (...
                    (P(down,j)+P(up,j)-2*P(i,j))/(delta_y^2) + ...
                    (P(i,left)+P(i,right)-2*P(i,j))/(delta_x^2));
                dZdt(i,j) = dZdt(i,j) +  dif_Z * (...
                    (Z(down,j)+Z(up,j)-2*Z(i,j))/(delta_y^2) + ...
                    (Z(i,left)+Z(i,right)-2*Z(i,j))/(delta_x^2));
                if params.detritus_layer
                    dDdt(i,j) = dDdt(i,j) +  dif_D * (...
                        (D(down,j)+D(up,j)-2*D(i,j))/(delta_y^2) + ...
                        (D(i,left)+D(i,right)-2*D(i,j))/(delta_x^2));
                end
            end

            % If the site is a landmass, maintain 0 concentration and 0
            % gradient
            if params.landmass(i,j)
                dNdt(i,j) = 0;
                dPdt(i,j) = 0;
                dZdt(i,j) = 0;
                if params.detritus_layer
                    dDdt(i,j) = 0;
                end
            end
        end
    end
    
    % Convert state variables from matrix-form to vector-form
    % N
    dNdt_vect = reshape(dNdt,[1,params.num_sites]);
    % P
    dPdt_vect = reshape(dPdt,[1,params.num_sites]);
    % Z
    dZdt_vect = reshape(dZdt,[1,params.num_sites]);
    if params.detritus_layer
        % D
        dDdt_vect = reshape(dDdt,[1,params.num_sites]);
        dxdt = [dNdt_vect, dPdt_vect, dZdt_vect, dDdt_vect]';
    else
        dxdt = [dNdt_vect, dPdt_vect, dZdt_vect]';
    end
end