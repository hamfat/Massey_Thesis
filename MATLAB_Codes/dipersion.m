syms v n v1 w d gl gk gca vl vk vca v2 v3 v4 psi c
%model parameters

%global d gl gk gca vl vk vca v2 v3 v4 psi c

%auxiliary functions

minf=0.5*(1+tanh((v-v1)/v2));
lamdan=cosh((v-v3)/(2*v4));
ninf=0.5*(1+tanh((v-v3)/v4));

dV=-gl*(v-vl)-gk*n*(v-vk)-gca*minf*(v-vca);

dN=psi*lamdan*(ninf-n);

k=jacobian([dV,dN],[v,n]);
 initial=[-0.28982,0.15767];% spiral source
%initial=[-0.65183,0.0034313]; % saddle
%initial=[-0.7141,0.00172]; %nodal sink
eq=initial;
K=[];
figure;
hold on;

%for v1=[ -0.25]
%for c=[0.005 0.00576 0.00618]
for w=linspace(-200,200,3000)
    d=0.0001;
    v1=-0.2466;
    v2=0.3125;
    v3=-0.1375;
    v4=.1812;
    vk=-1.125;
    vca=1.0;
    vl=-0.875;
    gk= 1.0; 
    gl=0.25;
    gca=0.4997;
    psi=0.3;
    %c=0.00611628125857;
    c=0.0043683356085;
    J=eval(subs(k,[v,n],[eq(1),eq(2)]));
    T=-d*w^2+2*w*1j*c+J(1,1)+J(2,2);
    Del=(-d*w^2+w*1j*c+J(1,1))*(w*1j*c+J(2,2))-(J(1,2)*J(2,1));
    p=[1 T Del];
    RTs=roots(p);
    realRTs = real(RTs);
    imRTs = imag(RTs);
   %plot(w,realRTs(1),'k.',w,realRTs(2),'k.',w,imRTs(1),'r.',w,imRTs(2),'r.')
   plot(-realRTs(1),imRTs(1),'b.',-realRTs(2),imRTs(2),'b.')
   %plot([0 0], [-2 2], 'k','Linewidth',1)
    
end

%end
