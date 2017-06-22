function [omega_set, lim, timesupport] = BSS_STRE(tfd, np, W, thr1, thr2)
omega_set = [];
lim = [];
timesupport = [];

%% Inputs Checkup
if(nargin == 0), fprintf(2,'ERROR (No Inputs)\n'); return;
elseif(nargin == 1)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    np = 3; W = 25; thr1 = 0.1; thr2 = 0.03;
elseif(nargin == 2)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(np)), np = 3;
    elseif(length(np)>1), fprintf(2,'ERROR (Number of Components must be a 1x1 number)\n'); return;
    elseif(~isnumeric(np)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(mod(np,1)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(np<=0), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    end
    W = 25; thr1 = 0.1; thr2 = 0.03;
elseif(nargin == 3)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(np)), np = 3;
    elseif(length(np)>1), fprintf(2,'ERROR (Number of Components must be a 1x1 number)\n'); return;
    elseif(~isnumeric(np)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(mod(np,1)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(np<=0), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    end
    if(isempty(W)), W = 25;
    elseif(length(W)>1), fprintf(2,'ERROR (Analysis Width must be a 1x1 integer)\n'); return;
    elseif(~isnumeric(W)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(mod(W,1)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(W<=0), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    end
    thr1 = 0.1; thr2 = 0.03;
elseif(nargin == 4)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(np)), np = 3;
    elseif(length(np)>1), fprintf(2,'ERROR (Number of Components must be a 1x1 number)\n'); return;
    elseif(~isnumeric(np)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(mod(np,1)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(np<=0), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    end
    if(isempty(W)), W = 25;
    elseif(length(W)>1), fprintf(2,'ERROR (Analysis Width must be a 1x1 integer)\n'); return;
    elseif(~isnumeric(W)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(mod(W,1)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(W<=0), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    end
    if(isempty(thr1)), thr1 = 0.1;
    elseif(length(thr1)>1), fprintf(2,'ERROR (Threshold Ceil must be a 1x1 number)\n'); return;
    elseif(~isnumeric(thr1)), fprintf(2,'ERROR (Threshold Ceil must must be a 1x1 number)\n'); return;
    elseif(thr1<=0), fprintf(2,'ERROR (Threshold Ceil must must be a 1x1 number)\n'); return;
    end
    thr2 = 0.03;
elseif(nargin == 5)
    if(isempty(tfd)), fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    if(isempty(np)), np = 3;
    elseif(length(np)>1), fprintf(2,'ERROR (Number of Components must be a 1x1 number)\n'); return;
    elseif(~isnumeric(np)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(mod(np,1)), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    elseif(np<=0), fprintf(2,'ERROR (Number of Components must must be a 1x1 number)\n'); return;
    end
    if(isempty(W)), W = 25;
    elseif(length(W)>1), fprintf(2,'ERROR (Analysis Width must be a 1x1 integer)\n'); return;
    elseif(~isnumeric(W)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(mod(W,1)), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    elseif(W<=0), fprintf(2,'ERROR (Analysis Width must must be a 1x1 integer)\n'); return;
    end
    if(isempty(thr1)), thr1 = 0.1;
    elseif(length(thr1)>1), fprintf(2,'ERROR (Threshold Ceil must be a 1x1 number)\n'); return;
    elseif(~isnumeric(thr1)), fprintf(2,'ERROR (Threshold Ceil must must be a 1x1 number)\n'); return;
    elseif(thr1<=0), fprintf(2,'ERROR (Threshold Ceil must must be a 1x1 number)\n'); return;
    end
    if(isempty(thr2)), thr2 = 0.03;
    elseif(length(thr2)>1), fprintf(2,'ERROR (Threshold Floor must be a 1x1 number)\n'); return;
    elseif(~isnumeric(thr2)), fprintf(2,'ERROR (Threshold Floor must must be a 1x1 number)\n'); return;
    elseif(thr2<=0), fprintf(2,'ERROR (Threshold Floor must must be a 1x1 number)\n'); return;
    end    
end
    
%% main
tfd = tfd./max(max(tfd));
[M, N] = size(tfd);
[all_comp, ~] = CompSep(tfd, N, M, W, thr1, thr2);
all_comp = all_comp(:,:,1:np);
omega_set = zeros(N, np);
temp = zeros(N, np);
lim = zeros(2, np);
for i = 1:np
    [~,I] = max(all_comp(:,:,i));
    omega_set(:,i) = I;
    [ll,~,~] = find(omega_set(:,i) > 1);
    lim(:,i) = [min(ll) ; max(ll)];   
    omega_set(lim(1,i):lim(2,i),i)=(omega_set(lim(1,i):lim(2,i),i))./(2*M);
    temp(:,i) = omega_set(:,i) > 1;
end
timesupport = sum(temp,2);
end

%% Auxiliary Functions
% This function separates the signal components from its time-frequency distribution (TFD)
function [tfd_all_kom,tfd_est_ostalo]=CompSep(tfd,tlength,flength,W,stop1,stop2)
N=tlength;
E_ukupno=sum(sum(tfd)); E_ostatak=E_ukupno;
    tfd_est_ostalo=tfd; 
    br=1;
    tfd_all_kom=zeros(flength,N,3);
while E_ostatak>=stop2*E_ukupno,
        tfd_est_kom=zeros(flength,N);
        [aaa,bbb]=max(tfd);
        [tfd_max_0,t_max_0]=max(aaa);
        f_max_br_0=bbb(t_max_0);
        tfd_slice=tfd(:,t_max_0); 
        D_f_minus=W; D_f_plus=W;
        if (f_max_br_0-D_f_minus)<1,    D_f_minus=f_max_br_0-1; end; 
        if (f_max_br_0+D_f_plus)>flength,     D_f_plus=flength-f_max_br_0;  end;
        window_brise_ostalo=zeros(flength,1); window_brise_comp=ones(flength,1); 
        window_brise_ostalo(f_max_br_0-D_f_minus:f_max_br_0+D_f_plus)=ones(max(size(f_max_br_0-D_f_minus:f_max_br_0+D_f_plus)),1); 
        window_brise_comp(f_max_br_0-D_f_minus:f_max_br_0+D_f_plus)=zeros(max(size(f_max_br_0-D_f_minus:f_max_br_0+D_f_plus)),1); 
        tfd_est_kom(:,t_max_0)=tfd_slice.*window_brise_ostalo ; 
        tfd_est_ostalo(:,t_max_0)=tfd_slice.*window_brise_comp; 
        t_max=t_max_0+1; f_max_br=f_max_br_0; tfd_max=tfd_max_0; 
        while t_max<=N && tfd_max>stop1*tfd_max_0  
            tfd_slice=tfd(:,t_max); 
            window_brise_ostalo=zeros(flength,1); window_brise_comp=ones(flength,1); 
            D_f_minus=W; D_f_plus=W;
            if (f_max_br-D_f_minus)<1, D_f_minus=f_max_br-1; end; 
            if (f_max_br+D_f_plus)>flength, D_f_plus=flength-f_max_br;  end;    
            window_brise_ostalo(f_max_br-D_f_minus:f_max_br+D_f_plus)=ones(max(size(f_max_br-D_f_minus:f_max_br+D_f_plus)),1); 
            window_brise_comp(f_max_br-D_f_minus:f_max_br+D_f_plus)=zeros(max(size(f_max_br-D_f_minus:f_max_br+D_f_plus)),1); 
            tfd_est_kom(:,t_max)=tfd_slice.*window_brise_ostalo ; %PPPPP
            tfd_est_ostalo(:,t_max)=tfd_slice.*window_brise_comp;  

            [tfd_max,f_max_br]=max(tfd_slice.*window_brise_ostalo);   
            t_max=t_max+1;
        end;
        t_max=t_max_0-1; f_max_br=f_max_br_0; tfd_max=tfd_max_0;
        while t_max>=1 && tfd_max>stop1*tfd_max_0  
            tfd_slice=tfd(:,t_max); 
            window_brise_ostalo=zeros(flength,1); window_brise_comp=ones(flength,1); 
            D_f_minus=W; D_f_plus=W;
            if (f_max_br-D_f_minus)<1, D_f_minus=f_max_br-1; end; 
            if (f_max_br+D_f_plus)>flength, D_f_plus=flength-f_max_br;  end;
            window_brise_ostalo(f_max_br-D_f_minus:f_max_br+D_f_plus)=ones(max(size(f_max_br-D_f_minus:f_max_br+D_f_plus)),1);
            window_brise_comp(f_max_br-D_f_minus:f_max_br+D_f_plus)=zeros(max(size(f_max_br-D_f_minus:f_max_br+D_f_plus)),1);
            tfd_est_kom(:,t_max)=tfd_slice.*window_brise_ostalo ;
            tfd_est_ostalo(:,t_max)=tfd_slice.*window_brise_comp;
            [tfd_max,f_max_br]=max(tfd_slice.*window_brise_ostalo);
            t_max=t_max-1;
        end;
    
    tfd=tfd_est_ostalo;    
    tfd_all_kom(:,:,br)=tfd_est_kom;
    br=br+1;    
    E_ostatak=sum(sum(tfd_est_ostalo));
end; 
end