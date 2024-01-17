clear all
close all
clc

max_lev = 5;
scheme = @self_adapting_EB_scheme;

nu = 13;
nv = 31;

a = 2;
b = 2.2;
c = 4;

k = 1; %choose between 1 and 12               
trunc = [];
switch k
   
    case 1 % Ellipsoid
        fx = @(theta,phi) a * cos( phi ) .* cos( theta );
        fy = @(theta,phi) b * cos( phi ) .* sin( theta );
        fz = @(theta,phi) c * sin( phi ) .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi/2,pi/2];
        
    case 2 % (N) Elliptic Paraboloid 
        fx = @(theta,phi) a * phi .* cos( theta );
        fy = @(theta,phi) b * phi .* sin( theta );
        fz = @(theta,phi) ( phi .^ 2 ) .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [0,2*pi];
    
    case 3 % Hyperbolic Paraboloid
        fx = @(theta,phi) cos( 0.*phi ) .* theta;
        fy = @(theta,phi) phi .* cos( 0.*theta );
        fz = @(theta,phi) phi .* theta;
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
        
    case 4 % Elliptic Hyperboloid (one sheet)
        fx = @(theta,phi) a * cosh( phi ) .* cos( theta );
        fy = @(theta,phi) b * cosh( phi ) .* sin( theta );
        fz = @(theta,phi) c * sinh( phi ) .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
        
    case 5 % Elliptic Hyperboloid (two sheet)
        fx = @(theta,phi) a * sinh( phi ) .* cos( theta );
        fy = @(theta,phi) b * sinh( phi ) .* sin( theta );
        fz = @(theta,phi) c + cosh( phi ) .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [0,pi];
        
    case 6 % (N) Elliptic Cone 
        fx = @(theta,phi) a * phi .* cos( theta );
        fy = @(theta,phi) b * phi .* sin( theta );
        fz = @(theta,phi) c * phi .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
        
    case 7 % Elliptic Cylinder
        fx = @(theta,phi) a * cos( 0.*phi ) .* cos( theta );
        fy = @(theta,phi) b * cos( 0.*phi ) .* sin( theta );
        fz = @(theta,phi) phi .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
        
    case 8 % Hyperbolic Cylinder
        fx = @(theta,phi) a * cos( 0.*phi ) .* cosh( theta );
        fy = @(theta,phi) b * cos( 0.*phi ) .* sinh( theta );
        fz = @(theta,phi) phi .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
         
    case 9 % Parabolic Cylinder
        fx = @(theta,phi) cos( 0.*phi ) .* theta;
        fy = @(theta,phi) cos( 0.*phi ) .* ( a * theta.^2 + b * theta + c );
        fz = @(theta,phi) phi .* cos( 0.*theta );
        range_1 = [-pi,pi];
        range_2 = [-pi,pi];
         
    case 10 % Mix1 = Elliptic Hyperboloid + Elliptic Cylinder + Sphere
        g1 = 2*pi;
        g2 = pi;
        
        fx=@(X,Y) sin(g2*(Y-1/2)).*cos(g1*X-pi/4).*(Y<=1+1e-4) + cos(g1*X-pi/4).*(Y>1+1e-4).*(Y<2-1e-4) + cosh(2*(-Y+2)).*cos(g1*X-pi/4).*(Y>=2-1e-4);
        fy=@(X,Y) sin(g2*(Y-1/2)).*sin(g1*X-pi/4).*(Y<=1+1e-4) + sin(g1*X-pi/4).*(Y>1+1e-4).*(Y<2-1e-4) + cosh(2*(-Y+2)).*sin(g1*X-pi/4).*(Y>=2-1e-4);
        fz=@(X,Y) - (0*X+(0.2+cos(g2*(Y-1/2))+1).*(Y<=1+1e-4) + (-2*(Y-1)+1).*(Y>1+1e-4).*(Y<2-1e-4) + (-0.2+sinh(2*(-Y+2))-1).*(Y>=2-1e-4));

        range_1 = [0,2];
        range_2 = [0.5,3];
         
    case 11 % Elliptic Paraboloid + Elliptic Cone + pseudo-Ellipsoid
        phi0 = 0.8;
        phi1 = 1.5;
        delay = phi1-atan(-1/phi1);
        delay2 = pi/8;
        fx = @(theta,phi) a * phi .* cos( theta+delay2 ).*(phi<phi0+1e-4)    + a * phi .* cos( theta+delay2 ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + -a/sin(phi1-delay)*cos( phi-delay ).* cos( theta+delay2 ).*(phi>phi1+1e-4);
        fy = @(theta,phi) b * phi .* sin( theta+delay2 ).*(phi<phi0+1e-4)     +b * phi .* sin( theta+delay2 ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + -b/sin(phi1-delay)*cos( phi-delay ).* sin( theta+delay2 ).*(phi>phi1+1e-4);
        fz = @(theta,phi) ((2*phi0*c - c*phi0^2) + c*( phi .^ 2 ) .* cos( 0.*theta )).*(phi<phi0+1e-4)  +2*c * phi .* cos( 0.*theta ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + (2*c*phi1-2*c*sin(phi1-delay)/cos(phi1-delay)+2*c/cos(phi1-delay) * sin( phi-delay )) .* cos( 0.*theta ).*(phi>phi1+1e-4);
        
        range_1 = [0,4*pi];
        range_2 = [0,phi1+pi/4];
        
        trunc(1) = cos(1*(range_1(2)-range_1(1))/(nu-1));
        trunc(2) = 1;
         
    case 12 % Elliptic Hyperboloid (two sheet) + Elliptic Cone + pseudo-Ellipsoid
        phi0 = 0.8;
        phi1 = 1.5;
        fac = 3;
        delay = phi1-atan(-1/phi1);
        delay2 = pi/8;
        fx = @(theta,phi) a*(phi0-tanh(fac*phi0)/fac+ 1/cosh(fac*phi0)/fac* sinh( fac*phi )) .* cos( theta+delay2 ).*(phi<phi0+1e-4)    + a * phi .* cos( theta+delay2 ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + -a/sin(phi1-delay)*cos( phi-delay ).* cos( theta+delay2 ).*(phi>phi1+1e-4);
        fy = @(theta,phi) b*(phi0-tanh(fac*phi0)/fac+ 1/cosh(fac*phi0)/fac* sinh( fac*phi )) .* sin( theta+delay2 ).*(phi<phi0+1e-4)     +b * phi .* sin( theta+delay2 ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + -b/sin(phi1-delay)*cos( phi-delay ).* sin( theta+delay2 ).*(phi>phi1+1e-4);
        fz = @(theta,phi) 2*c*(phi0-1/tanh(fac*phi0)/fac + 1/sinh(fac*phi0)/fac*cosh( fac*phi )) .* cos( 0.*theta ).*(phi<phi0+1e-4)  +2*c * phi .* cos( 0.*theta ).*(phi>=phi0+1e-4).*(phi<=phi1+1e-4) + (2*c*phi1-2*c*sin(phi1-delay)/cos(phi1-delay)+2*c/cos(phi1-delay) * sin( phi-delay )) .* cos( 0.*theta ).*(phi>phi1+1e-4);
        
        range_1 = [0,4*pi];
        range_2 = [0,phi1+pi/4];
        
        trunc(1) = cos(1*(range_1(2)-range_1(1))/(nu-1));
        trunc(2) = cosh(fac*(range_2(2)-range_2(1))/(nv-1));
end

u = linspace(range_1(1),range_1(2),nu);
v = linspace(range_2(1),range_2(2),nv)';

X = fx(u,v);
Y = fy(u,v);
Z = fz(u,v);

new_X = scheme(X,max_lev);
new_Y = scheme(Y,max_lev);
new_Z = scheme(Z,max_lev);

plotMesh(new_X,new_Y,new_Z)