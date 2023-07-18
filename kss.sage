#
# Sage Code to search for KSS curves with D=1 or D=3
# M.Scott February 2018

# set K value here
K=12
# set CM discriminant, must be 1 or 3
D=3

# Set search parameters. We will be searching through polynomials of degree nb
# ... with nz non-zero coefficients of absolute size less than or equal to lim  

nz=3  # must be > 1
lim=1 # must be positive

# Show All solutions for each rho, or just first solution, which may not be the prettiest.
showall=False

# some successful searches
# D=1, K=8,  nz=2, lim=1
# D=3, K=12, nz=3, lim=1   # Finds BLS and BN curves!
# D=1, K=16, nz=2, lim=2
# D=3, K=18, nz=2, lim=2
# D=1, K=32, nz=2, lim=3
# D=3, K=36, nz=2, lim=2
# D=1, K=40, nz=2, lim=2
# D=3, K=54, nz=2, lim=1

def igcd(x,y) :
# integer GCD, returns GCD of x and y 
    if y==0 :
        return x
    while True :
        r=x%y
        if r==0 :
            break
        x=y
        y=r
    return y

def mylcm(a,b) :
    return (a*b)/gcd(a,b)

def flat(p) :
    lcm=1
    c=p.coefficients()
    #print c
    #print p.degree()
    for i in range(len(c)) :
        d=c[i].denom()
        lcm=mylcm(lcm,d)
    return lcm
    
def content(p) :
    con=1
    c=p.coefficients()
    con=c[0][0]
    #print c
    for i in range(1,len(c)) :
        d=c[i][0]
        con=gcd(con,d)
    return Integer(con)

def iter_bits(x,n) :
    gotone=False 
    for i in range(0,n-1) :
        if x[i]==1 and x[i+1]==0 :
            gotone=True
            x[i+1]=1
            x[i]=0
            if x[0]==1 :
                break
            k=1
            while True :
                if x[k] != 0 :
                    break
                k=k+1
            for j in range(0,i-k) : 
                x[j]=x[j+k]
                x[j+k]=0
            break            
    return gotone

def iter_nums(v,m,lim) :
    for k in range(0,m) :
        if v[k]==-1 : 
            v[k]=1
            break
        if v[k]<lim :
            v[k]=v[k]+1
            break
        v[k]=-lim
    for k in range(0,m) :
        if v[k]!=lim :
            return True
    return False

u=[]
q=[]
s=[]
v=[]
rhobestn=8  # only interested in rho<1.6
rhobestd=5

EL=0
if D==1 :
    EL=mylcm(4,K)
if D==3 :
    EL=mylcm(3,K)

nb=euler_phi(EL)
#print("nb= ",nb)

C.<z> = CyclotomicField(EL)
#R.<x> = PolynomialRing(ZZ)

for i in range(0,nb) :
    q.append(0)
    u.append(0)

for i in range(0,nz) :
    s.append(0)
    v.append(0)

more_bits=True
k=0
while True :
    if k==0 :
        for j in range(nz) :
            u[j]=1
    else :
        more_bits=iter_bits(u,nb)
    j=0
    for i in range(nb) :
        if u[i]!=0 :
            s[j]=i
            j=j+1
    more_nums=True
    n=0
    while True :
        if n==0 :
            for j in range(nz) :
                v[j]=-lim
        else:
            more_nums=iter_nums(v,nz,lim)
        
        for j in range(nb) :
            q[j]=0
        for j in range(nz) :
            q[s[j]]=v[j]
        
# all above here manages the search loop            

        pb=0    # create next polyomial in QQ
        for j in range(nb) :
            pb=pb+q[j]*z^j

        print(q,end="\r")

        r=pb.minpoly()
        #print("\nr= ",r)
        M.<w> = NumberField(r)
        rz=M.roots_of_unity()
        nrz=len(rz)
        #print("\nRoots = ",nrz)

        # get CM discriminant D as a polynomial

        isrmd=1/(0*w-D)
        #print("isrmd= ",isrmd)

        if isrmd.is_square() :
            sd=isrmd.sqrt()
            #print("\nd= ",sd*sd)
            for i in range(nrz) :

            # search though K-th roots of unity

                if igcd(i+1,K) != 1 :
                    continue;
                pru=rz[i]^(EL/K)    # k-th root of unity
                #print("\nunity= ",pru^K)
                ft=pru+1
                fy=sd*(pru-1)
    
                t=ft.polynomial()
                y=fy.polynomial()

                p=(t*t+D*y*y)/4

                rhon=p.degree()
                rhod=r.degree()
                ig=igcd(rhon,rhod)
                rhon/=ig
                rhod/=ig

                if rhon<rhod :
                    continue

                if showall :
                    if rhobestd*rhon>rhobestn*rhod :
                        continue  # rho is not interesting
                else :
                    if rhobestd*rhon>=rhobestn*rhod :
                        continue  # rho is not interesting

                if not p.is_irreducible() :
                    continue

            #print ("rho= ",rhon,"/",rhod)

            # solution looks interesting...
            # convert polynomials over QQ to ZZ (with one common integer divisor m)

                plcm=flat(p)
                tlcm=flat(t)

                m=mylcm(plcm,tlcm)
                p=p*m
                t=t*m

                b=0
                tries1=0
                tries2=0
                fail=False;
                #print ("True")
                while True :
                    # try to find any residue class that works

                    tries1=tries1+1
                    if tries1>200000 :  # give up..
                        fail=True
                        break
                    if p(b)%m != 0 :
                        b=b+1
                        continue
                    if t(b)%m !=0 :
                        b=b+1
                        continue

                    tries1=0
                    sp=p(x=m*x+b).expand()/m

                    # try 100 times to find p that doesn't have an integer factor
                
                    if content(sp)!=1 :
                        tries2=tries2+1
                        if tries2>100 :
                            fail=True
                            break
                        b=b+1
                        continue
                    tries2=0
                    if not sp.polynomial(ZZ).is_irreducible() :
                        b=b+1
                        continue

                    if fail :
                        break
                    st=t(x=m*x+b).expand()/m
                    sr=r(x=m*x+b).expand()
 
                    c=content(sr)
                    sr=sr/c
                    if not sr.polynomial(ZZ).is_irreducible() :
                        b=b+1
                        continue
                    break

                if fail :
                    continue

                ct=gcd(content(t),m)
                mt=m/ct
                t=t/ct
                sp=sp.polynomial(ZZ)
                st=st.polynomial(ZZ)
                sr=sr.polynomial(ZZ)

                np=sp+1-st

# check its the right embedding degree!
                isitreal=True
                for j in range(2,K) :
                    if K%j!=0 :
                        continue
                    if ((sp^j)-1)%sr==0 :
                        isitreal=False
                        #print("That one is a dud, actual embedding degree is ",j)
                        #print ("rho= ",rhon,"/",rhod)
                        break
                if not isitreal :
                    continue

                #if ((sp^K)-1)%sr!=0 :
                #    continue;

                if np%sr == 0 :
                    cf=np/sr
                    print ("\nSolution found, rho= ",rhon,"/",rhod)
                    print ("p= (",p,")/",m)
                    print ("t= (",t,")/",mt)
                    print ("r= ",r)
                    print ("For sample residue class ",m,"*x +",b)
                    print ("p= ",sp)
                    print ("t= ",st)
                    print ("r= ",sr)
                    print ("c= ",cf)
                    print ("i= ",i)
                    print ("0?= ",((sp^K)-1)%sr); # check embedding degree
                    print ("")

                    rhobestn=rhon 
                    rhobestd=rhod

# all below here manages the search loop    

        if not more_nums :
            break
        n=n+1
    if not more_bits :
        break
    k=k+1

    