function [pp,gv_alpha,gv_beta,gv_delta] = Initialization_vagwo(np,nx,varmax,varmin,velmax,velmin)
pp=zeros(np,nx); 
gv_alpha=zeros(np,nx);
gv_beta=zeros(np,nx);
gv_delta=zeros(np,nx);
for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
    gv_alpha(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_beta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_delta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
end