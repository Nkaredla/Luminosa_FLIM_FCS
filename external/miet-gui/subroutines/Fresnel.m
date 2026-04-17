function [rp,rs,tp,ts] = Fresnel(w1,n1,n2)

% [rp,rs,tp,ts] = Fresnel(w1,n1,n2) calculates Fresnel reflection (r) and
% transmission (t) coefficients for p-wave and s-wave, either for crossing
% one interface or the effective coefficients when crossing a stack of 
% media. Since these coefficients depend on the angle of incidence of the
% wave onto the interface, the component w1 of the wave vector perpendicular
% to the interface in the first medium (e.g. the medium where the emitter 
% is situated) is also needed. The component q parallel to the interface is
% the same in all media. We assume that the length scale is chosen in such 
% a way that the wave vector in vacuum,k_0=omega/c, has the length 1. Then,
% k_1=n_1, k_2=n_2, n_1^1=q^2+w_1^2 and so on.
% 
% CAREFUL: TWO COMPLETELY DIFFERENT USES OF THE INPUT VARIABLES!
%
% CASE 1 -- n1 and n2 are scalars: n1 = n_1 and n2 = n_2
%
% CASE 2 -- n1 and n2 are vectors with length(n1)==length(n2)+2:
% n1 = vector of indices of refraction for all the layers
% n2 = thicknesses of all layers except the top and bottom one (these two
%      are assumed to be infinitely thick)

w1 = w1(:).';                   % make sure w1 is a row vector [x x x ... x]
n1 = n1(:).'; n2 = n2(:).';     % make sure n1 and n2 are row vectors

% CASE 1: one interface between media 1 and 2
if length(n1)==1 && length(n2)==1   
    w2 = sqrt(n2^2-n1^2+w1.^2);     % w_2=sqrt(n_2^2-q^2)=sqrt(n_2^2-(n_1^1-w_1^2))
    w2(imag(w2)<0) = conj(w2(imag(w2)<0));  % should not be necessary: imag(w2)>=0 anyway
    rp = (w1*n2^2-w2*n1^2)./(w1*n2^2+w2*n1^2);
    rs = (w1-w2)./(w1+w2);
    tp = 2*n1*n2*w1./(w1*n2^2+w2*n1^2);
    ts = 2*w1./(w1+w2);
    ind= (w1==0 & w2==0);   % only happens if n1==n2, would give lots of NANs
    rp(ind)=0; rs(ind)=0; tp(ind)=1; ts(ind)=1; % set to physically meaningful values
% CASE 2: stack of different materials
elseif length(n1)==length(n2)+2     
    if length(n1)==2    % same situation as in case 1: calculate same way
        [rp,rs,tp,ts] = Fresnel(w1,n1(1),n1(2));
    else    % more than one interface
        n = n1;         % better names for variables: refr. indices in n,
        d = [0 n2 0];   % layer thicknesses in d
        % matrix w: one row for each layer, one column for each propagation direction
        w(1,:) = w1;   
        for j = 2:length(n)    
            w(j,:) = sqrt(n(j)^2-n(1)^2+w1.^2); % q^2=n(1)^2-w1^2=n(j)^2-wj^2
            w(j,imag(w(j,:))<0) = conj(w(j,imag(w(j,:))<0)); % should not be necessary: imag(w2)>=0 anyway       
        end  
        
    % p-wave: calculate transition matrix M for traversing the stack
        j = length(n);  % start at layer furthest away from the dipole 
        % interface between layer j and layer j-1
        M11 = (w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
        M21 = (-w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;
        inftyalert=(ones(1,length(w1))<0);    % initialise vector as "false"
        for j=length(n)-1:-1:2
            % phase shift aquired when traversing layer j ~exp(+/-i*w*d);
            % if w is imaginary, phase shift dampens or amplifies matrix
            % elements, that can lead to numerical problems
            % => realmin=2.2251e-308 and realmax=1.7977e+308 or
            %    log(realmin)=-708.3964 and log(realmax)=709.7827
            % => choose 708 as critical exponent (small safety margin)
            % => if exp(-iw*d)>largest number MATLAB can cope with, results
            % will be useless: set inftyalert=true so tp will be zero 
            % & divide both matrix elements by the exponent because they
            % will cancel in rp anyway
            % => if exp(iw*d)<smallest number MATLAB can cope with,
            % results will be useless: set M21 to zero
            % => within numerical error, the same is true if M11*exp(-iw*d)>realmax
            inftyIndex=((imag(w(j,:))*d(j)>708)|(imag(w(j,:))*d(j)+log(M11)>708));
            M11(~inftyIndex) = exp(-1i*w(j,~inftyIndex)*d(j)).*M11(~inftyIndex);
            M21(~inftyIndex) = exp(1i*w(j,~inftyIndex)*d(j)).*M21(~inftyIndex);
            M21(inftyIndex)=0;
            inftyalert(inftyIndex)=true;
            % transfer matrix for interface between layers j and j-1
            N11 = (w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2;  % N22=N11
            N21 = (-w(j,:)./w(j-1,:)*n(j-1)/n(j)+n(j)/n(j-1))/2; % N12=N21
            % effective transfer matrix for getting from layers length(n) to j-1
            tmp11 = N11.*M11+N21.*M21;
            tmp21 = N21.*M11+N11.*M21;
            M11 = tmp11;
            M21 = tmp21;
            
        end
     % p-wave: get reflection and transmission coefficients from transfer matrix
        rp = M21./M11;
        tp = 1./M11;
        tp(inftyalert)=0;
        
    % s-wave: calculate transition matrix M for traversing the stack    
        j = length(n);  % start at layer furthest away from the dipole 
        % interface between layer j and layer j-1
        M11 = (w(j,:)./w(j-1,:)+1)/2;
        M21 = (-w(j,:)./w(j-1,:)+1)/2;
        inftyalert=(ones(1,length(w1))<0);
        for j=length(n)-1:-1:2
            % phase shift aquired when traversing layer j
            inftyIndex=((imag(w(j,:))*d(j)>708)|(imag(w(j,:))*d(j)+log(M11)>708));
            M11(~inftyIndex) = exp(-1i*w(j,~inftyIndex)*d(j)).*M11(~inftyIndex);
            M21(~inftyIndex) = exp(1i*w(j,~inftyIndex)*d(j)).*M21(~inftyIndex);
            M21(inftyIndex)=0;
            inftyalert(inftyIndex)=true;
            % transfer matrix for interface between layers j and j-1
            N11 = (w(j,:)./w(j-1,:)+1)/2;
            N21 = (-w(j,:)./w(j-1,:)+1)/2;
            % effective transfer matrix for getting from layers length(n) to j-1
            tmp11 = N11.*M11+N21.*M21;
            tmp21 = N21.*M11+N11.*M21;
            M11 = tmp11;
            M21 = tmp21;
        end
    % s-wave: get reflection and transmission coefficients from transfer matrix
        rs = M21./M11;
        ts = 1./M11;
        ts(inftyalert)=0;
    end
else
    error('Wrong input');
end

