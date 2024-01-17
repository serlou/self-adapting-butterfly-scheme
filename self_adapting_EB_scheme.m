function new_P = self_adapting_EB_scheme(P,max_lev,trunc,w)

tol = 10^-6;
if nargin<3
trunc = [0.5,1.1];
end

% apply the self-adapting EB scheme to data P on a regular grid
% max_lev times

if nargin < 4
    w = 0;
end

if nargin < 2
    max_lev = 3;
end

% Extended Butterfly rule
a = 1/2-w;
b = 1/8+2*w;
c = -1/16-w;
d = w;
EB_rule = [ a, a, b, c, d, c, b, c, d, c ];

% stencils template: beware that we have the Xs as columns and
% Ys as rows
h_stencil = fliplr( [
    -1 0;       %p0
    1 0;        %p1
    -1 -2;      %p01
    -3 -2;      %p02
    -3 0;       %p03
    -1 2;       %p04
    1 2;        %p10
    3 2;        %p11
    3 0;        %p12
    1 -2        %p13
    ] );

v_stencil = fliplr( [
    0 -1;       %p0
    0 1;        %p1
    2 1;        %p01
    2 -1;       %p02
    0 -3;       %p03
    -2 -3;      %p04
    -2 -1;      %p10
    -2 1;       %p11
    0 3;        %p12
    2 3         %p13
    ] );

d_stencil = fliplr( [
    -1 -1;      %p0
    1 1;        %p1
    1 -1;       %p01
    -1 -3;      %p02
    -3 -3;      %p03
    -3 -1;      %p04
    -1 1;       %p10
    1 3;        %p11
    3 3;        %p12
    3 1         %p13
    ] );

% subdivision level 0
new_P = P;
[new_row,new_col] = size(new_P);
new_flags = ones(new_row,new_col);

count_hv = 0;
count_d = 0;
count_EB = 0;

% further subdivision step(s)
for j = 1 : max_lev
    
    old_P = new_P;
    old_row = new_row;
    old_col = new_col;
    old_flags = new_flags;
    
    new_row = 2*old_row - 1;
    new_col = 2*old_col - 1;
    new_P = zeros( new_row, new_col );
    new_flags = new_P;
    
    % interpolation
    new_P(1:2:end,1:2:end) = old_P;
    new_flags(1:2:end,1:2:end) = old_flags;
    
    for row_idx = 1 : new_row
        
        for col_idx = 1 : new_col
            
            temp = mod([row_idx,col_idx],2);
            check_orientation = temp * [2,1]';
            if check_orientation ~= 3
                
                switch check_orientation
                    
                    case 2 % horizontal
                        stencil = h_stencil;
                        
                    case 1 % vertical
                        stencil = v_stencil;
                        
                    case 0 % diagonal
                        stencil = d_stencil;
                        
                end
                
                
                idxs = [ row_idx, col_idx ] + stencil;
                if min(min(idxs)) < 1 || max(idxs(:,1)) > new_row || max(idxs(:,2)) > new_col
                    continue
                end
                
                p_loc = zeros(10,1);
                flag_loc = 1;
                for k = 1 : 10
                    
                    p_loc(k) = new_P(idxs(k,1),idxs(k,2));
                    flag_loc = flag_loc * new_flags(idxs(k,1),idxs(k,2));
                    
                end
                
                % if all the points involved are eligible we check for gamma
                % otherwise we use the extended butterfly
                
                if flag_loc
                    if abs( p_loc(7) - p_loc(3) ) > tol
                        phi_1st = ( p_loc(6) - p_loc(4) + p_loc(8) - p_loc(10) ) / ( 2*p_loc(7) - 2*p_loc(3) );
                    else
                        phi_1st = 1;
                    end
                    
                    if abs( p_loc(1) - p_loc(2) ) > tol
                        
                        phi_2nd = ( p_loc(7) - p_loc(6) - p_loc(3) + p_loc(10) ) / ( 2*p_loc(2) - 2*p_loc(1) );
                        
                        phi_1 = ( p_loc(10) - p_loc(3) + p_loc(7) - p_loc(6) ) / ( 2 * p_loc(2) - 2 * p_loc(1) );
                        phi_2 = ( p_loc(8) - p_loc(7) + p_loc(3) - p_loc(4) ) / ( 2 * p_loc(2) - 2 * p_loc(1) );
                        
                    else
                        
                        phi_2nd = 1;
                        
                        phi_1 = 1;
                        phi_2 = 1;
                        
                    end
                    
                    %                                Horizontal
                    %                         	    p(6)    p(7)    p(8)
                    %                         p(5)	p(1)    p(2)    p(9)
                    %                         p(4)	p(3)    p(10)
                    
                    check_hv = abs([
                        p_loc(1) - p_loc(2) - p_loc(5) + p_loc(9) + ( 2*p_loc(1) - 2*p_loc(2) ) * phi_1st %(univariate)
                        p_loc(2) - p_loc(4) + p_loc(5) - p_loc(10) + ( -2*p_loc(1) + 2*p_loc(3) )* phi_1st
                        p_loc(1) - p_loc(4) + p_loc(9) - p_loc(10) + ( -2*p_loc(2) + 2*p_loc(3) ) * phi_1st
                        p_loc(6) - p_loc(5) - p_loc(2) + p_loc(8) + ( 2*p_loc(1) - 2*p_loc(7) ) * phi_1st
                        p_loc(7) - p_loc(6) - p_loc(3) + p_loc(10) + ( 2*p_loc(1) - 2*p_loc(2) ) * phi_2nd
                        p_loc(6) - p_loc(1) + p_loc(8) - p_loc(9) + ( 2*p_loc(2) - 2*p_loc(7) ) * phi_1st
                        p_loc(6) - p_loc(4) + p_loc(8) - p_loc(10) + ( 2*p_loc(3) - 2*p_loc(7) ) * phi_1st %(standard Butterfly)
                        ]);
                    check_hv = sort(check_hv);
                    
                    %                                Diagonal
                    %                                       p(8)    p(9)
                    %                               p(7)    p(2)    p(10)
                    %                         p(6)  p(1)    p(3)
                    %                         p(5)  p(4)
                    
                    check_d = abs([
                        p_loc(3) - p_loc(4) - p_loc(7) + p_loc(8) + ( 2*p_loc(1) - 2*p_loc(2) ) * phi_2
                        p_loc(7) - p_loc(6) - p_loc(3) + p_loc(10) + ( 2*p_loc(1) - 2*p_loc(2) ) * phi_1
                        p_loc(5) - p_loc(2) - p_loc(1) + p_loc(9) + ( p_loc(3) - p_loc(4) + p_loc(7) - p_loc(8) ) * phi_1 + ( p_loc(3) - p_loc(6) + p_loc(7) - p_loc(10) ) * phi_2 %(extended Butterfly)
                        p_loc(6) - p_loc(1) - p_loc(8) + p_loc(9) + ( 2*p_loc(3) - 2*p_loc(1) ) * phi_1 + ( p_loc(2) - p_loc(4) + p_loc(5) - p_loc(10) ) * phi_2 + ( 2*p_loc(1) - 2*p_loc(6) ) * phi_2^2 + ( 2*p_loc(7) - 2*p_loc(2) ) * phi_1 * phi_2
                        p_loc(4) - p_loc(1) + p_loc(9) - p_loc(10) + ( p_loc(2) + p_loc(5) - p_loc(6) - p_loc(8) ) * phi_1 + ( 2*p_loc(7) - 2*p_loc(1) ) * phi_2 + ( 2*p_loc(1) - 2*p_loc(4) ) * phi_1^2 + ( 2*p_loc(3) - 2*p_loc(2) ) * phi_1 * phi_2
                        ]);
                    check_d = sort(check_d);
                    
                    scaling = max( abs( [
                        p_loc(1) - p_loc(2);
                        p_loc(1) - p_loc(3);
                        p_loc(1) - p_loc(4);
                        p_loc(1) - p_loc(5);
                        p_loc(1) - p_loc(6);
                        p_loc(1) - p_loc(7);
                        p_loc(2) - p_loc(3);
                        p_loc(2) - p_loc(7);
                        p_loc(2) - p_loc(8);
                        p_loc(2) - p_loc(9);
                        p_loc(2) - p_loc(10);
                        ] ) );
                    
                    
                    [ check_min, check_rule ] = min( [check_hv(end), check_d(end)] );
                    if check_min >= tol * scaling
                        [ check_min, check_rule ] = min( [check_hv(2), check_d(2)] );
                    end
                    if check_min < tol * scaling
                        
                        switch check_rule
                            
                            case 1 %horizontal/vertical
                                phi_1st = max(min(phi_1st,trunc(2)),trunc(1));
                                phi_2nd = max(min(phi_2nd,trunc(2)),trunc(1));
                                if min( [phi_1st,phi_2nd] + 1 ) > tol
                                    
                                    new_P(row_idx,col_idx) = hv_rule(w,phi_1st) * p_loc;
                                    new_flags(row_idx,col_idx) = 1;
                                    count_hv = count_hv + 1;
                                    
                                end
                                
                            case 2 %diagonal
                                phi_1 = max(min(phi_1,trunc(2)),trunc(1));
                                phi_2 = max(min(phi_2,trunc(2)),trunc(1));
                                if min( [phi_1, phi_2] + 1 ) > tol && abs( phi_1 + phi_2 - 2 ) > tol
                                    
                                    new_P(row_idx,col_idx) = d_rule(w,phi_1,phi_2) * p_loc;
                                    new_flags(row_idx,col_idx) = 1;
                                    count_d = count_d + 1;
                                    
                                end
                                
                        end
                        
                        if not( new_flags(row_idx,col_idx) ) || isnan( new_P(row_idx,col_idx) )
                            new_flags(row_idx,col_idx) = 0;
                            flag_loc = 0;
                        end
                        
                    else
                        flag_loc = 0;
                    end
                end
                
                if not( flag_loc )
                    
                    new_P(row_idx,col_idx) = EB_rule * p_loc;
                    count_EB = count_EB + 1;
                    
                    if check_min < 1e-12
                        new_flags(row_idx,col_idx) = 1;
                    end
                end
            end
            
        end
        
    end
    
    % cutting the outer unnecessary parts of new_P
    new_P = new_P(3:end-2,3:end-2);
    new_flags = new_flags(3:end-2,3:end-2);
    new_row = new_row - 4;
    new_col = new_col - 4;
    
    trunc = sqrt((1+trunc)/2);
end
end


%% New Rules

function HV = hv_rule(w,phi)

hva = w * ( 1 - 2*phi ) + 1 / ( 2 * sqrt((1+phi)/2) );
hvb = phi * ( 16*w + ...
    1 / ( sqrt((1+phi)/2) * (1+sqrt((1+phi)/2))/2 ) ) / 8;
hvc = - w - 1 / ( 16 * sqrt((1+phi)/2) * (1+sqrt((1+phi)/2))/2 );
hvd = w;
HV = [ hva, hva, hvb, hvc, hvd, hvc, hvb, hvc, hvd, hvc ];

end

function D = d_rule(w,phi_1,phi_2)

da = ( - 2 * w * ( -4*phi_1*phi_2 + phi_1 + phi_2 + 2*phi_1^2-1 + 2*phi_2^2-1 ) ...
    - ( phi_1 + phi_2 + 2 ) / ( 2 * sqrt((1+phi_1)/2) * sqrt((1+phi_2)/2) )...
    + phi_1 + phi_2 ) / ( 2 * ( phi_1 + phi_2 - 2 ) );
db = ( 4 * w * ( 2 - 2*phi_1 - 2*phi_2 + 2*phi_1^2-1 + 2*phi_2^2-1 )...
    + ( phi_1 + phi_2 ) / ( sqrt((1+phi_1)/2) * sqrt((1+phi_2)/2) ) ...
    - 2 )/ ( 4 * ( phi_1 + phi_2 - 2 ) );
dc = ( 4 * w * ( -2 * phi_1 * phi_2 + phi_1 + phi_2 ) ...
    + 1 / ( sqrt((1+phi_1)/2) * sqrt((1+phi_2)/2) ) - 1 ) ...
    / ( 4 * ( phi_1 + phi_2 - 2 ) );
dd = w;
D = [ da, da, db, dc, dd, dc, db, dc, dd, dc ];

end