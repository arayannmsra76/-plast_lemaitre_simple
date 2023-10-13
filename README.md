# -plast_lemaitre_simple
function [D] = plast_lemaitre_simple(i, ig)
%% Parameters 
global mat;
global eldat;
tol = 1e-7;
dim = eldat.compute_dim;
if dim == 3
% array for hydraulic component retrieval
Ih = [1;1;1;0;0;0];
Ih2 = [0; 0;0;1;1;1];
% 2nd order identity tensor in voigt notation
I = [1;1;1;0;0;0];
% conversion arrays
enToPhys = [1;1;1;.5;.5;.5];
physToEn = [1;1;1;2;2;2];
% number of unique stresses
ncomp = 6;
%linear elastic material tangent
Ce = lin_material (i,0);

% 2nd order identity tensor in Voigt notation
I = [1; 1; 1; 0; 0; 0];
elseif dim == 2
% array for hydraulic component retrieval
Ih = [1;1;0;1];
Ih2 = [0; 0;1;0];
%conversion arrays
enToPhys = [1;1;.5;1];
physToEn = [1;1;2;1];
% number of unique stresses
ncomp = 4;
% linear elastic material tangent
Ce = lin_material (i,1);
% 2nd order identity tensor in Voigt notation
I = [1; 1; 0; 1];
end
% tensor product of two 2nd order identity tensors in Voigt notation
IxI = I*I';
% 4th order symmetric projection tensor in Voigt notation
Is = 0.5.*(diag(I) + eye( ncomp));
% deviatoric projection tensor in Voigt notation
Id = Is - (1./3) .* IxI;

% material parameters (see input_lemaitre_lmat.m)
matrnr = eldat. mat (i);
sig_y0 = mat. sig( matnr ); % initial yield stress
G = mat. G( matnr);         % shear modulus
K = mat.K( matrnr );          % bulk modulus
r = mat.r ( matrnr );         % isotropic hardening
s = mat.s ( matnr);         % isotropic hardening
Rinf = mat. Rinf ( matnr);  % isotropic hardening
gamma = mat. gamma (matnr); % isotropic hardening

% some factors
K2 = 2.*K;
G2 = 2.*G;
G3 = 3.*G;
G6 = 6.*G;
%% Get Strain, Damage, Hardening
eps = eldat. epsilon(i, :, ig)';        % total strain
eps_pl = eldat. eps_pl(i, :, ig)';      % plastic strain
Dam0 = eldat. damage (i, ig);            % initial damage
Int0 = 1 - Dam0;                        % initial integrity
R0 = eldat. R(i, ig);                   % initial hardening internal variable
%% Trial State
eps_e_tr = eps - eps_pl;                % elastic trial strain 
eps_e_hyd_tr = sum(eps_e_tr.*Ih);       % hydrostatic strain

sig_hyd_tr = K.*eps_e_hyd_tr;                        % hydrostatic effective stress
eps_e_dev_tr = eps_e_tr - (1./3).*eps_e_hyd_tr.* Ih;  % deviatoric strain
% convert engineering shear strain to physical shear strain
eps_e_dev_tr = eps_e_dev_tr .* enToPhys;
% compute effective trial von Mises stress
temp = eps_e_dev_tr.^2;
J_2 = G2.^2.* (.5*sum(temp.*Ih) + sum (temp.*Ih2) );
q_tr = sqrt (3*J_2);
% compute yield stress
fsig0 = sig_y0 + Rinf.*(1-exp (-gamma.*R0));
Phi = q_tr - fsig0;
%% Check if Yield Criterion is met 
if Phi >= 0
    % return mapping
    % initial guess for the plastic multiplier
    plasticMult = Int0.*Phi./ (3.*G);

    % inital guess for the hardening variable
    R = R0 + plasticMult;

    % initial values
    norm_F = 1;
    sig_hyd_tr2 = sig_hyd_tr.^2;

    % Newton Raphson iteration for finding the true plastic multiplier
    iter = 0;
    maxiter = 100;
while norm_F >= tol && iter <= maxiter
if iter == maxiter
    disp('Fatal: return mapping reached maximum iterations!');
    break;
end

%current yield stress
fsig = sig_y0 + 3300.*(1-exp (-.4.*R));

% integrity & strain energy release rate function
f1 = (G3 ./ (q_tr - fsig));
Int = f1 .* plasticMult;
Y = -(fsig.^2)./G6 - (sig_hyd_tr2)./K2;
f2 = -Y./r;

% compute residual
F = Int - Into + (-Y./r).^s ./f1;
norm_F = abs (F);
% derivatives

dfsig = 1320 .* exp(-.4.*(R));
dY = -(fsig.*dfsig) ./ G3;

% residual derivative
f= f1 + f1.*plasticMult.*dfsig./(q_tr - fsig) - ...
(dfsig./G3).*f2.^s - (s.*dY./(f1.*r)).*f2.^(s-1) ;

% get next plastic multiplier
plasticMult = plasticMult - F./f;

% update hardening variable
R = R0 + plasticMult;
iter = iter + 1;

end

% having now received the true plastic multiplier, update 
% effective yield stress
fsig = sig_y0 + Rinf.*(1-exp(-gamma.*R));
% hardening slope
dfsig = gamma.* Rinf .* exp(-gamma.*R) ;
% integrity
f1 = (G3 ./ (q_tr - fsig));
Int = f1 .* plasticMult;
dInt = (G3+Int.*dfsig)./(q_tr-fsig);
% strain energy release rate
Y =-(fsig.^2)./G6 - (sig_hyd_tr2)./K2;
dY = -(fsig.*dfsig) ./ G3;
f2 = -Y./r;
% residual derivative
f= f1 + f1.*plasticMult.*dfsig./(q_tr - fsig) - ...
    (dfsig./G3).*f2.^s - (s.*dY./(f1.*1)).*f2.^(s-1);
% check if NR yielded an acceptable damage variable 
if (Int < 1e-20)
    disp('GP integrity too small!');
end
% update damage
Dam = 1-Int;

% update true von Mises stress
q = Int .* fsig;
% update true stresses, strains
sig_hyd = Int .* sig_hya_tr;
sig_dev = G2.*(q./q_tr).*eps_e_dev_tr;
sig = sig_dev + sig_hyd .* Ih;

% plastic corrector
plCor = G3.*plasticMult./(Int.* q_tr);
% restore engineering strain and update total elastic strain
eps_e = (1-plCor) .* eps_e_dev_tr .* physToEn + (1./3).*eps_e_hyd_tr .* Ih;

% increase in plastic strain
deps_pl = eps_e_tr - eps_e;

% update model
eldat. sigma (i, :,ig) = sig;
eldat.eps_pl_n(i,:, ig) = eps_pl + deps_pl;
eldat. damage_n (i, ig) = Dam;
eldat. R_n(i, ig) = R;
eldat. q_n(i, ig) = q;
eldat. triax_n(i, ig) = sig_hyd./q;
eldat. sigy_n(i, ig) = fsig;

%norm of deviatoric stress
snorm = sqrt (sum (physToEn.*(sig-dev.^2)));

% compute coefficients for elastoplastic tangent
f3 = q_tr - fsig;
a1 = (1./f).*(Int./13 - (1./G3).*f2.^s );
a2 = -s.*sig_hyd_tr.*f3./(G3.*r.*K.*f).*f2.^(s-1);
a3 = a2.*dInt;
a4 = a1.*dInt - Int./f3;
a = G2.*Int*fsig./q_tr;
b= G2.*(a1.*dfsig.*Int + a4.*fsig - Int.*fsig./q_tr );
b = b./(snorm.^2);
c = K.*sqrt (2./3).*(a2.*dfsig.*Int + a3.*fsig );
c = c./ snorm;
d = G2.*sqrt(3./2).*sig_hyd_tr.*a4;
d = d./ snorm;
e = K.*(Int + a3.*sig_hyd_tr);
% Consistent Elastoplastic Tangent Operator
% (minor symmetric, not major symmetric)
D= a.*Id + b.*sig_dev*sig_dev' + c.sig_dev*I' + d.*I*sig_dev' + e.*IxI;

% If yield criterion is not met, material is behaving elastically
else
 D = Int0.*Ce;
% true stress components
sig_hyd = Int0 .* sig_hyd_tr;
sig_dev = G2.* Int0.* eps_e_dev_tr;

% true stress
sig = sig_dev + sig_hyd .* Ih;
% update model
    eldat.sigma (i, :,ig) = sig;
%plastic strain
    eldat. eps_pl_n(i,:, ig) = eps_pl;
% initial damage
    eldat. damage_n (i, ig) = Dam0;
% isotropic hardening internal variable
    eldat. R_n(i, ig) = R0;
% true vM stress
    eldat. q_n(i, ig) = Int0.*q_tr;
% stress triaxiality

    eldat. triax_n(i, ig) = sig_hyd./(Int0.*q_tr);
% effective yield stress
    eldat. sigy_n(i, ig) = fsig0;
end




