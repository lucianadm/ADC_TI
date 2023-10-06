function  k=logistico(varargin) 
m=1;
n=1; 
persistent semilla1;
if (isempty(semilla1))
    semilla1=rand;
end
if (nargin==1)
 n  = varargin{1};
m  = varargin{2};
end
if (nargin==2)
 n  = varargin{1};
 m  = varargin{2};
end
if (nargin==3)
 n  = varargin{1};
 m  = varargin{2};
 semilla1=varargin{3};
end
x=semilla1;  
clear k;
for (col=1:m)
for (oo=1:n)
k(oo,col)=(-4*(x-1/2).^2+1);
 x=k(oo,col);
 end
end
semilla1=k(n,m); 
