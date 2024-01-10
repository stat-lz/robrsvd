function [u, s, v, diagout]=robfsvdls(data, paramstruct);
%robfsvdls robust functional singular value decomposition 
%
% Usage:
%   [u, s, v, diagout]=robfsvdls(data, paramstruct);
%
% Inputs:
%
%   data        an input data matrix
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    irobust           1 use huber function to robustify the
%                         singular value decomposition
%
%                      0 (default) use power algorithm to calculate
%                         singular value decomposition
%
%    huberk            when irobust=1, huberk provides the parameter for
%                      the huberfunction. Default value is 1.345, which is
%                      very robust in most case
%
%                      Note that when irobust=0, this option won't work
%
%    iinituv           0 if there is no initial estimate for (s, u, v), we
%                        will use the standard SVD to estimate sigma
%                        (default)
%
%                      1  you need to specify s, u, v by your self, make
%                         the following inits, initu and initv are
%                         specified.
%
%    inits             the initial estimate of s
%   
%    initu             the initial estimate of u
%
%    initv             the initial estimate of v
%
%    niter             the maximum iteration for the iterative reweight
%                      least square (default value is 1000)
%
%    tol               the tolerance value for the iteration, default value
%                      is 1e-5;
%
%    istabilize        1 (default), the program will first scale the input
%                         matrix, so that the largest absolute value of the
%                         entries is 1. This may reduce the numerical
%                         instability in many cases
%
%                      0 directly work on the input data. Please make sure
%                        that your input data is well scaled, so that the
%                        program will yield Inf or NaN
%
%    uspar             the smoothing parameter to smooth the left singular
%                      vector (or singular column). The default value is 0,
%                      i.e., do not smooth u.
%
%    vspar             the smoothing parameter to smooth the right singular
%                      vector (or singular row. The default value is 0,
%                      i.e., do not smooth v. 
%
%                      these two smoothing parameters should be nonnegative
%
%    iugcv             to select the smoothing parameter for u by gcv
%                      (default value is 0, i.e. not select the smoothing
%                      parameter by gcv. Please specify 1 if gcv is needed)
%
%    ivgcv             to select the smoothing parameter for v by gcv
%                      (default value is 0, i.e. not select the smoothing
%                      parameter by gcv. Please specify 1 if gcv is needed)
%
%
%    usparmax          when iugcv is 1, a maximum u value may be needed.
%                      the default value is 1000
%
%    usparmin          when iugcv is 1, a minimum u value may be needed.
%                      the default value is 1e-10
%
%    nuspar            number of smoothing parameters of u to search by gcv,
%                      when iugcv is 1, this argument may be needed. The
%                      default value is 14.
%
%    iloguspar         whether the number of grids in gcv is on log scale
%                      or linear scale. default value is 1. The default
%                      uspar vector varies from 10.^[-10:4]'; If specify by
%                      0, the value will be linspace(1e-10, 10000, 14)';
%
%
%    vsparmax          when ivgcv is 1, a maximum v value may be needed.
%                      the default value is 1000
%
%    vsparmin          when ivgcv is 1, a minimum v value may be needed.
%                      the default value is 1e-10
%
%    nvspar            number of smoothing parameters of v to search by gcv,
%                      when ivgcv is 1, this argument may be needed. The
%                      default value is 14.
%
%    ilogvspar         whether the number of grids in gcv is on log scale
%                      or linear scale. default value is 1. The default
%                      vspar vector varies from 10.^[-10:4]'; If specify by
%                      0, the value will be linspace(1e-10, 10000, 14)';
%
%   
% Outputs:
%
%    u                 the singular column or the left singular vector
%
%    s                 the singular value
%
%    v                 the singular row or the right singular vector
%
%    diagout	       a data structure which contains several diagnosis metrics
%
%                      ugcvscore contains the last step gcv scores for u
%
%                      vgcvscore contains the last step gcv scores for v
%
%(c) Copyright Lingsong Zhang (lingsong@purdue.edu)
%
%  Please make sure the following functions are available in the same
%  directory of robfsvdls.m
%
%    ssmatls.m
%    huberweightls.m
%
%


irobust=0;
huberk=1.345;
niter=1000;
tol=1e-5;
istabilize=1;
uspar=0;
vspar=0;
iinituv=0;
iugcv=0;
ivgcv=0;
usparmin=1e-10;
usparmax=1e4;
nuspar=14;
iloguspar=1;
vsparmin=1e-10;
vsparmax=1e4;
nvspar=14;
ilogvspar=1;

ugcvmat=[];
vgcvmat=[];

if nargin>1;
    if isfield(paramstruct, 'irobust');
        irobust=getfield(paramstruct, 'irobust');
    end;
    
    if isfield(paramstruct, 'huberk');
        huberk=getfield(paramstruct, 'huberk');
    end;
    
    if isfield(paramstruct, 'niter');
        niter=getfield(paramstruct, 'niter');
    end;
    
    if isfield(paramstruct, 'tol');
        tol=getfield(paramstruct, 'tol');
    end;
    
    if isfield(paramstruct, 'uspar');
        uspar=getfield(paramstruct, 'uspar');
    end;
    
    if isfield(paramstruct, 'vspar');
        vspar=getfield(paramstruct, 'vspar');
    end;
    
    if isfield(paramstruct, 'iinituv');
        iintuv=getfield(paramstruct, 'iinituv');
        if isfield(paramstruct, 'inits');
            inits=getfield(paramstruct, 'inits');
        else;
            error('Please specify an initial value of s');
        end;
        
        if isfield(paramstruct, 'initu');
            initu=getfield(paramstruct, 'initu');
        else;
            error('Please specify an initial value of u');
        end;
        
        if isfield(paramstruct, 'initv');
            initv=getfield(paramstruct, 'initv');
        else
            error('Please specify an initial value of v');
        end;      
    end;
    
    if isfield(paramstruct, 'iugcv');
        iugcv=getfield(paramstruct, 'iugcv');
    end;
    
    if isfield(paramstruct, 'ivgcv');
        ivgcv=getfield(paramstruct, 'ivgcv');
    end;
    
    if isfield(paramstruct, 'usparmin');
        usparmin=getfield(paramstruct, 'usparmin');
    end;
    
    if isfield(paramstruct, 'usparmax');
        usparmax=getfield(paramstruct, 'usparmax');
    end;
    
    if isfield(paramstruct, 'nuspar');
        nuspar=getfield(paramstruct, 'nuspar');
    end;
    
    if isfield(paramstruct, 'iloguspar');
        iloguspar=getfield(paramstruct, 'iloguspar');
    end;
    
    if isfield(paramstruct, 'vsparmin');
        vsparmin=getfield(paramstruct, 'vsparmin');
    end;
    
    if isfield(paramstruct, 'vsparmax');
        vsparmax=getfield(paramstruct, 'vsparmax');
    end;
    
    if isfield(paramstruct, 'nvspar');
        nvspar=getfield(paramstruct, 'nvspar');
    end;
    
    if isfield(paramstruct, 'ilogvspar');
        ilogvspar=getfield(paramstruct, 'ilogvspar');
    end;

end;


[m, n]=size(data);
if istabilize==1;
    myscale=1.4785*mad(mat2vec(data), 1);
    localdata=data./myscale;
else;
    myscale=1;
    localdata=data;
end;


if iinituv==0;
    [uold, sold, vold]=svds(localdata, 1);
else;
    uold=initu;
    sold=inits;
    vold=initv;
    if istabilize==1;
        sold=sold/myscale;
    end;
end;
%uold;
%vold;
% for uniqueness, we always keep v as norm 1
uold=sold.*uold;



Appold=uold*vold';
Rmat=localdata-Appold;
Rvec=zeros(m*n, 1);
Rvec(:)=Rmat(:, :);
mysigma=median(abs(Rvec))/0.675;
clear Rvec;
iter=1;
localdiff=9999;
diffvec=[];
uspar_current=uspar;
vspar_current=vspar;

ugcvscore=[];
vgcvscore=[];

while (localdiff>tol &iter<niter); %starting the iterative reweighting
    if (irobust==1);
        % we will use the Huber's function
        Wmat=huberweightls(Rmat./mysigma, huberk);
    else
        Wmat=ones(m, n); %then this algorithm becomes the power algorithm
    end;


    
    %to update u first
    
    if iugcv==0;
        uterm1=diag(sum(diag(vold.^2)*Wmat'))+(2*mysigma^2)*(vold'*(speye(n)+vspar.*ssmatls(n))*vold*(speye(m)+uspar.*ssmatls(m))-norm(vold)^2*speye(m));
        uterm2=(Wmat.*localdata)*vold;
        unew=inv(uterm1)*uterm2;
    elseif iugcv==1;
        if nuspar<0;
            error('number of smoothing parameter can not be negative');
        else 
            if iloguspar==1;
%                usparvec=10.^linspace(log10(usparmin), log10(usparmax), nuspar)';
                usparvec=logspace(log10(usparmin), log10(usparmax), nuspar)';
				
            else
                usparvec=linspace(usparmin, usparmax, nuspar)';
            end;
            %end of specify uspar vector
            
            ugcvvec=[];
			ugcvmat=[];
            for (iter_uspar=1:nuspar);
                u_nsrobust=inv(diag(sum(diag(vold.^2)*Wmat')))*(Wmat.*localdata)*vold;
                usterm1=diag(sum(diag(vold.^2)*Wmat'))+(2*mysigma^2)*(vold'*(speye(n)+vspar_current.*ssmatls(n))*vold*(speye(m)+usparvec(iter_uspar).*ssmatls(m))-norm(vold)^2*speye(m));
                usterm2=(Wmat.*localdata)*vold;
                u_srobust=inv(usterm1)*usterm2;
                smooth_u=inv(usterm1)*diag(sum(diag(vold.^2)*Wmat'));
                gcv_ut=m*norm(u_nsrobust-u_srobust)^2/(m-trace(smooth_u))^2;
                ugcvvec=[ugcvvec; gcv_ut];
				ugcvmat=[ugcvmat, u_srobust./norm(u_srobust)];
            end;
            uspar_current=max(usparvec(ugcvvec==min(ugcvvec)));        
            
	    ugcvscore=[usparvec, ugcvvec];

            uterm1=diag(sum(diag(vold.^2)*Wmat'))+(2*mysigma^2)*(vold'*(speye(n)+vspar_current.*ssmatls(n))*vold*(speye(m)+uspar_current.*ssmatls(m))-norm(vold)^2*speye(m));
            uterm2=(Wmat.*localdata)*vold;
            unew=inv(uterm1)*uterm2;            
        end;
    else
        error('iugcv should be either 0 or 1!');
    end;
    
    %to update v secondly
    
    if ivgcv==0;
        vterm1=diag(sum(diag(unew.^2)*Wmat))+(2*mysigma^2)*(unew'*(speye(m)+uspar.*ssmatls(m))*unew*(speye(n)+vspar.*ssmatls(n))-norm(unew)^2*speye(n));
        vterm2=(Wmat.*localdata)'*unew;
        vnew=inv(vterm1)*vterm2;
    elseif ivgcv==1;
        if nvspar<0;
            error('number of smoothing parameter can not be negative');
        else 
            if ilogvspar==1;
%                vsparvec=10.^linspace(log10(vsparmin), log10(vsparmax), nvspar)';
                vsparvec=logspace(log10(vsparmin), log10(vsparmax), nvspar)';
            else
                vsparvec=linspace(vsparmin, vsparmax, nvspar)';
            end;
            %end of specify uspar vector
            
            vgcvvec=[];
			vgcvmat=[];
            for (iter_vspar=1:nvspar);
                v_nsrobust=inv(diag(sum(diag(unew.^2)*Wmat)))*(Wmat.*localdata)'*unew;
                vsterm1=diag(sum(diag(unew.^2)*Wmat))+(2*mysigma^2)*(unew'*(speye(m)+uspar_current.*ssmatls(m))*unew*(speye(n)+vsparvec(iter_vspar).*ssmatls(n))-norm(unew)^2*speye(n));
                vsterm2=(Wmat.*localdata)'*unew;
                v_srobust=inv(vsterm1)*vsterm2;
                smooth_v=inv(vsterm1)*diag(sum(diag(unew.^2)*Wmat));
                gcv_vt=n*norm(v_nsrobust-v_srobust)^2/(n-trace(smooth_v))^2;
                vgcvvec=[vgcvvec; gcv_vt];
				vgcvmat=[vgcvmat, v_srobust./norm(v_srobust)];
            end;
            vspar_current=max(vsparvec(vgcvvec==min(vgcvvec)));        

	    vgcvscore=[vsparvec, vgcvvec];
            
            vterm1=diag(sum(diag(unew.^2)*Wmat))+(2*mysigma^2)*(unew'*(speye(m)+uspar_current.*ssmatls(m))*unew*(speye(n)+vspar_current.*ssmatls(n))-norm(unew)^2*speye(n));
            vterm2=(Wmat.*localdata)'*unew;
            vnew=inv(vterm1)*vterm2;            
        end;       
    else
        error('ivgcv should be either 0 or 1!');
    end;
        
    
    Appnew=unew*vnew';
    Rmat=localdata-Appnew;
    
    localdiff=max(max(abs(Appnew-Appold)));
    Appold=Appnew;
    uold=norm(vnew).*unew;
    vold=vnew./norm(vnew);
    %this one makes v is always with norm 1
    
    iter=iter+1;
    diffvec=[diffvec; localdiff];
end;

v=vold;
s=myscale*norm(uold);
u=uold./norm(uold);

if iugcv==1;
    uspar=uspar_current;
end;
if ivgcv==1;
    vspar=vspar_current;
end;

diagout=struct('ugcvscore', ugcvscore, 'vgcvscore', vgcvscore, 'ugcvmat', ugcvmat, 'vgcvmat', vgcvmat);

