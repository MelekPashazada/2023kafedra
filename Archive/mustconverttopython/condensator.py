mm=[0,40]
dp=[0,40]
r=[0,40]
kk=2
ak=2
p0=1

p0=400
m0=0.2
patm=1
bet=0.78
bc=2.5*0.0001
tau=1
hh=20
rk=1000
rc=0.1
N=10
ak=0
kk=0.05*0.000000000001
q=1000000
qq=q*(1/(2*3.14*hh))
taum=tau
taut=0.001

print(mm)

def zx(x):
    zx=1-0.1162*0.01*x+0.3744*0.00001*x*x-0.2965*0.000000001*x*x*x-0.1975*0.00000000001*x*x*x*x


def zpp(x):
    zpp=-0.1162*0.01+2*0.3744*0.00001*x-3*0.2965*0.000000001*x*x-4*0.1975*0.00000000001*x*x*x

def kp(x):
    kp=kk*ak**(x-p0)

def sp(x):
    sp=0.814286*x

def spp(x):
    spp=0.814286

def ap(x):
    ap=1.116+0.1157*0.01*x+0.23674*0.000001*x*x

def app(x):
    app=0.1157*0.01+2*0.23674*0.000001*x

def cp(x):
    cp=(0.637*0.0001-0.5057*0.000001*x+0.6265*0.00000001*x*x-0.1595*0.0000000001*x*x*x+0.13*0.0000000000001*x*x*x*x)

def cpp(x):
    cpp=(-0.5057*0.000001+2*0.6265*0.00000001*x-3*0.1595*0.0000000001*x*x+4*0.13*0.0000000000001*x*x*x)

def gmap(x):
    gmap=194.899-0.42974*0.1*x+0.1335*0.0001*x*x-0.6053*0.000001*x*x*x+0.622*0.000000001*x*x*x*x

def gampp(x):
    gampp=-0.42974*0.1+2*0.1335*0.0001*x-3*0.6053*0.000001*x*x+4*0.622*0.000000001*x*x*x

def mugp(x):
    mugp=((0.0126+0.257*0.0001*x+0.1633*0.0000001*x*x)*0.00000001/(86400))

def mukp(x):
    mukp=((0.6+0.3295*0.01*x+0.1044*0.0001*x*x-0.1558*0.0000001*x*x*x+0.85*0.00000000001*x*x*x*x)*0.00000001/(86400))

def fq(x):
    fq=2.0833*x*x*x*x+4.9167*x*x*x-5.5708*x*x-0.277*x+0.882

def fk(x):
    fk=1.8864*x*x+0.1889*x+0.0005

def nn1(x):
    nn1=((x*bet*(1-cp(x)*gamp(x)))/(mugp(x)*zx(x)))*kp(x)

def nn2(x):
    nn2=((x*bet*cp(x))/(mugp(x)*zx(x)))*kp(x) 

def mm1(x):
    mm1=(sp(x)/(mukp(x)*ap(x)))*kp(x) 

def mm2(x):
    mm2=(1/(mukp(x)*ap(x)))*kp(x)
    
def ll1(x):
    ll1=(x*bet*(1-cp(x)*gamp(x))/zx(x))

def ll11(x):
    ll11=(1/sqr(zx(x)))*(bet*zx(x)*((1-cp(x)*gamp(x))-x*(cpp(x)*gamp(x)+cp(x)*gampp(x)))-x*bet*(1-cp(x)*gamp(x))*zpp(x))

def ll2(x):
    ll2=(x*bet*cp(x))/zx(x)

def ll22(x):
    ll22=(1/sqr(zx(x)))*(bet*zx(x)*(cp(x)+x*cpp(x))-bet*x*cp(x)*zpp(x))

def q1(x):
    q1=(sp(x)/ap(x))

def q2(x):
    q2=(1/ap(x))

def q11(x):
    q11=(1/sqr(ap(x)))*(spp(x)*ap(x)-sp(x)*app(x))

def q22(x):
    q22=-(1/sqr(ap(x)))*app(x)

def ff0(x):
    ff0=(q1(x)-ll1(x))/(q2(x)-ll2(x))   

def ff1(x,y,z):
    ff1=ll11(x)*z*(1-y)+q11(x)*z*y

def ff2(x,y,z):
    ff2=ll22(x)*z*(1-y)+q22(x)*z*y

def ff3(x,y,z):
    ff3=ll1(x)*(1-y)+q1(x)*y

def ff4(x,y,z):
    ff4=ll2(x)*(1-y)+q2(x)*y

def ff5(x,y,z):
    ff5=z*q2(x)-z*ll2(x)

def f1g(x,y,z):
    f1g=0.5*(fq(y)*nn1(x)+fq(y1)*nn1(x1))

def f2g(x,x1,y,y1):
    f2g=0.5*(fq(y)*nn2(x)+fq(y1)*nn2(x1))

def f1k(x,x1,y,y1):
    f1k=0.5*(fk(y)*mm1(x)+fk(y1)*mm1(x1))

def f2k(x,x1,y,y1):
    f2k=0.5*(fk(y)*mm2(x)+fk(y1)*mm2(x1))

for i in range(0,n):
    r[i]=i*ln(rk/rc)/n
for i in range(0,n):
    p[i,o]=po
    sig[i,o]=0
    m[i,o]=mo

h=ln(rk/rc)/n
j=1
k1=1
k=0
k2=1

while j>0:
    as1=(f1g(p[0,j-1],p[0,j-1],sig[0,j-1],sig[0,j-1])+f1k(p[0,j-1],p[0,j-1],sig[0,j-1],sig[0,j-1]))
    a[0,j]=1
    c[0,j]=-qq*h/as1

    for i in range(1,n-1):
        aa=(exp(-2*r[i])*taut/((rc*rc)*h*h))
        m[i,j]=m0*exp(bc*(p[i,j-1]-p0))

        df1=aa*(f1g(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1])+f1k(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1]))-aa*ff0(p[i,j-1])*(f2g(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1])+f2k(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1]))

        df2=aa*(f1g(p[i,j-1],p[i-1,j-1],sig[i,j-1],sig[i,j-1])+f1k(p[i,j-1],p[i-1,j-1],sig[i,j-1],sig[i,j-1]))-aa*ff0(p[i,j-1])*(f2g(p[i,j-1],p[i-1,j-1],sig[i,j-1],sig[i,j-1])+f2k(p[i,j-1],p[i-1,j-1],sig[i,j-1],sig[i,j-1]))

        a[i,j]=df1/(df1+df2*(1-a[i-1,j])+(ff1(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff2(p[i,j-1],sig[i,j-1],m[i,j-1])))

        c[i,j]=(df2*c[i-1,j]+p[i,j-1]*(ff1(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff2(p[i,j-1],sig[i,j-1],m[i,j-1]))-(ff3(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff4(p[i,j-1],sig[i,j-1],m[i,j-1]))*(m[i,j]-m[i,j-1]))/(df1+df2*(1-a[i-1,j])+(ff1(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff2(p[i,j-1],sig[i,j-1],m[i,j-1])));


    p[n-1,j]=c[n-1,j]/(1-a[n-1,j])
    p[n,j]=p[n-1,j]
    i=n-1


while i!=0:
    p[i-1,j]=a[i-1,j]*p[i,j]+c[i-1,j]
    i=i-1


m[0,j]=m0*exp(bc*(p[0,j-1]-p0))

aa=(exp(-2*r[0])*taut/((rc*rc)*h*h))

st2=aa*(f2g(p[1,j-1],p[2,j-1],sig[1,j-1],sig[2,j-1])+f2k(p[1,j-1],p[2,j-1],sig[1,j-1],sig[2,j-1]))*(p[2,j]-p[1,j])-aa*(f2g(p[0,j-1],p[1,j-1],sig[0,j-1],sig[1,j-1])+f2k(p[0,j-1],p[1,j-1],sig[0,j-1],sig[1,j-1]))*(p[1,j]-p[0,j])

st1=aa*(f1g(p[1,j-1],p[1,j-1],sig[1,j-1],sig[1,j-1])+f1k(p[1,j-1],p[1,j-1],sig[1,j-1],sig[1,j-1]))*(p[2,j-1]-p[1,j-1])-aa*(f1g(p[0,j-1],p[0,j-1],sig[0,j-1],sig[0,j-1])+f1k(p[0,j-1],p[0,j-1],sig[0,j-1],sig[0,j-1]))*(p[1,j-1]-p[0,j-1])

i=0

stt=st1-ff0(p[0,j-1])*st2-(ff3(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff4(p[i,j-1],sig[i,j-1],m[i,j-1]))*(m[i,j]-m[i,j-1])


sig[i,j]=sig[i,j-1]+(st2/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))-(ff2(p[i,j-1],sig[i,j-1],m[i,j-1])/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))*(1/(ff1(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff2(p[i,j-1],sig[i,j-1],m[i,j-1])))*stt-(ff4(p[i,j-1],sig[i,j-1],m[i,j-1])/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))*(m[i,j]-m[i,j-1])

i=1
while i==n:
    if i==1:
        sdd=sig[0,j]
    else :
        sdd=sig[i,j-1]
        
    m[i,j]=m0*exp(bc*(p[i,j-1]-p0))
    aa=(exp(-2*r[i])*taut/((rc*rc)*h*h))


    st2=aa*(f2g(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1])+f2k(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1]))*(p[i+1,j]-p[i,j])-aa*(f2g(p[i,j-1],p[i-1,j-1],sig[i,j-1],sdd)+f2k(p[i,j-1],p[i-1,j-1],sig[i,j-1],sdd))*(p[i,j]-p[i-1,j])

    st1=aa*(f1g(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1])+f1k(p[i,j-1],p[i+1,j-1],sig[i+1,j-1],sig[i+1,j-1]))*(p[i+1,j]-p[i,j])-aa*(f1g(p[i,j-1],p[i-1,j-1],sig[i,j-1],sdd)+f1k(p[i,j-1],p[i-1,j-1],sig[i,j-1],sdd))*(p[i,j]-p[i-1,j])

    stt=st1-ff0(p[i,j-1])*st2-(ff3(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff4(p[i,j-1],sig[i,j-1],m[i,j-1]))*(m[i,j]-m[i,j-1])

    sig[i,j]=sig[i,j-1]+(st2/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))-(ff2(p[i,j-1],sig[i,j-1],m[i,j-1])/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))*(1/(ff1(p[i,j-1],sig[i,j-1],m[i,j-1])-ff0(p[i,j-1])*ff2(p[i,j-1],sig[i,j-1],m[i,j-1])))*stt-(ff4(p[i,j-1],sig[i,j-1],m[i,j-1])/ff5(p[i,j-1],sig[i,j-1],m[i,j-1]))*(m[i,j]-m[i,j-1])
    
    i=i+1

sig[n,j]=sig[n-1,j]











if trunc((taut*j)/10)==k2:
    gk=(sp(p[0,j])+(mukp(p[0,j])/mugp(p[0,j]))*(fk(sig[0,j])/fq(sig[0,j]))*ap(p[0,j])*(p[0,j]/zx(p[0,j]))*bet*(1-cp(p[0,j])*gamp(p[0,j])))/(1+cp(p[0,j])*(mukp(p[0,j])/mugp(p[0,j]))*(fk(sig[0,j])/fq(sig[0,j]))*ap(p[0,j])*(p[0,j]/zx(p[0,j])))

    print((taut*j):2:2,' ',(m[0,j]):4:4,' ',(p[0,j]):2:2,' ',(p[1,j]):2:2,' ',(p[2,j]):2:2,' ',(p[3,j]):2:2,' ',(p[4,j]):2:2,' ',(p[5,j]):2:2,' ',(p[6,j]):2:2,' ',(p[7,j]):2:2,' ',(p[8,j]):2:2,' ',(p[9,j]):2:2,' ',(p[10,j]):2:2,' ',(p[11,j]):2:2,' ',(p[12,j]):2:2,' ',(p[13,j]):2:2,' ',(p[14,j]):2:2,' ',(p[15,j]):2:2,' ',(sig[0,j]):5:5,' ',(sig[1,j]):5:5,' ',(sig[2,j]):5:5,' ',(sig[3,j]):5:5,' ',sig[4,j]:5:5,' ',sig[5,j]:5:5,' ',sig[6,j]:5:5,' ',sig[7,j]:5:5,' ',sig[8,j]:5:5,' ',sig[9,j]:5:5,' ',sig[10,j]:5:5,' ',(sig[11,j]):5:5,' ',(sig[12,j]):5:5,' ',(sig[13,j]):5:5,' ',(sig[14,j]):5:5,' ',sig[15,j]:5:5,' ',(q/gk):2:2,' ',gk:2:2)

    print("\n\n")


    print(f,(taut*j):2:2,' ',' ',(m[0,j]):4:4,' ',(p[0,j]):2:2,' ',(p[1,j]):2:2,' ',(p[2,j]):2:2,' ',(p[3,j]):2:2,' ',(p[4,j]):2:2,' ',(p[5,j]):2:2,' ',(p[6,j]):2:2,' ',(p[7,j]):2:2,' ',(p[8,j]):2:2,' ',(p[9,j]):2:2,' ',(p[10,j]):2:2,' ',(p[11,j]):2:2,' ',(p[12,j]):2:2,' ',(p[13,j]):2:2,' ',(p[14,j]):2:2,' ',(p[15,j]):2:2,' ',(sig[0,j]):5:5,' ',(sig[1,j]):5:5,' ',(sig[2,j]):5:5,' ',(sig[3,j]):5:5,' ',sig[4,j]:5:5,' ',sig[5,j]:5:5,' ',sig[6,j]:5:5,' ',sig[7,j]:5:5,' ',sig[8,j]:5:5,' ',sig[9,j]:5:5,' ',sig[10,j]:5:5,' ',(sig[11,j]):5:5,' ',(sig[12,j]):5:5,' ',(sig[13,j]):5:5,' ',(sig[14,j]):5:5,' ',sig[15,j]:5:5,' ',(m[0,j]):5:5,' ',m[n,j]:5:5,' ',(q/gk):2:2,' ',gk:2:2)

    print(f)
    k2=k2+1

j=j+1 
close(f)
assign(f,'d:\rell.txt')
rewrite(f)

    



