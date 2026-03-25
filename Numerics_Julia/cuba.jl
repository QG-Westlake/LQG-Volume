using Cubature
using LinearAlgebra


function xl(g)
    (g[1,1]+g[2,2])/2
end

function lambda1(x)
    x+sqrt(x^2-1)
end

function it(x)
    return sin(x[1])*x[2]^2
end

function it2(x)
    return (sin(x[1])^2)*sin(x[2])
end

function sn(nl,nu,f1)
    rs=0
    for i in nl:nu
        rs=rs+f1
    end
end

function pr(n)
    if n==0
        return 1.0
    elseif n>0
        x1=1.0
        for x in 1:n
            x1=x1*x
        end
        return x1
    end
end

function prlg(n)
    if n==0
        return 1.0
    elseif n>0
        x1=0.0
        for x in 1:n
            x1=x1+log(x)
        end
        return x1
    end
end



#j-representation of SL(2,C) element#
function glc(j,g)
    a=floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).+floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).*im
    for x1 in -j:j
        for x2 in -j:j
            for l in 0:2j
                if j+x1-l>=0&&j-x2-l>=0&&-x1+x2+l>=0&&l>=0
                    a[Int32(j+x1+1),Int32(j+x2+1)]=a[Int32(j+x1+1),Int32(j+x2+1)]+(sqrt(pr(j+x1))*sqrt(pr(j-x1))*sqrt(pr(j+x2))*sqrt(pr(j-x2)))/(pr(j+x1-l)*pr(j-x2-l)*pr(-x1+x2+l)*pr(l))*(g[1,1]^(j-x2-l))*(g[2,2]^(j+x1-l))*(g[1,2]^(-x1+x2+l))*(g[2,1]^(l))
                end
            end
        end
    end
    return a
end

function glc2(j,g)
    a=floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).+floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).*im
    for x1 in -j:j
        for x2 in -j:j
            for l in 0:2j
                if j+x1-l>=0&&j-x2-l>=0&&-x1+x2+l>=0&&l>=0
                    a[Int32(j+x1+1),Int32(j+x2+1)]=a[Int32(j+x1+1),Int32(j+x2+1)]+1/(pr(j+x1-l)*pr(j-x2-l)*pr(-x1+x2+l)*pr(l))*(g[1,1]^(j-x2-l))*(g[2,2]^(j+x1-l))*(g[1,2]^(-x1+x2+l))*(g[2,1]^(l))
                end
            end
            a[Int32(j+x1+1),Int32(j+x2+1)]=log(a[Int32(j+x1+1),Int32(j+x2+1)])+(prlg(j+x1)+prlg(j-x1)+prlg(j+x2)+prlg(j-x2))/2
        end
    end
    return a
end

#j-representation of SU(2) element#
function hp(j,s)
    sp=sqrt(s[1]^2+s[2]^2+s[3]^2)
    if sp!=0
        rep1=[cos(sp)+im*s[3]*sin(sp)/sp (im*s[1]+s[2])*sin(sp)/sp;im*(s[1]+im*s[2])*sin(sp)/sp cos(sp)-im*s[3]*sin(sp)/sp]
    else
        rep1=[1 0;0 1]
    end 
    return glc(j,rep1)
end

#SU(2) element#
function hpp123(s1,s2,s3,rep1)
    rep1.=[cos(s1)+im*cos(s2)*sin(s1) (im*cos(s3)+sin(s3))*sin(s1)*sin(s2);(im*cos(s3)-sin(s3))*sin(s1)*sin(s2) cos(s1)-im*cos(s2)*sin(s1)]
    #return 0
    #else
    #    rep1=[1.0+im*0.0 0.0+im*0.0;0.0+im*0.0 1.0+im*0.0]
    #end 
end

#\tilde{\theta}
function ttk1(x)
    vec1 = x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    vecz=[0 z1*sin(theta) z1*cos(theta)]
    vecw=[0 0 w1]
    h1=hp(1,vec1)*transpose(vecz)
    #println("h1=",h1)
    mw=sqrt(conj(vecw)*transpose(vecw))
    #println("mw=",mw)
    mz=sqrt(transpose(conj(h1))*h1)
    #println("mz=",mz)
    a1=conj(vecw)*h1/(mw*mz)
    return a1[1]
end

#\tilde{\chi}
function ttk2(x)
    vec1 = x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    vecw=[0 w2*sin(chi) w2*cos(chi)]
    vecz=[0 0 z2]
    h1=hp(1,vec1)*transpose(vecz)
    mw=sqrt(conj(vecw)*transpose(vecw))
    mz=sqrt(transpose(conj(h1))*h1)
    a1=conj(vecw)*h1/(mw*mz)
    return a1[1]
end
    
function f1(z, w, x)
    acosh(cos(z)*cos(conj(w))+sin(z)*sin(conj(w))*ttk1(x))
end

function f2(z, w, x)
    acosh(cos(z)*cos(conj(w))+sin(z)*sin(conj(w))*conj(ttk2(x)))
end

#Integrand for open 4 vertex
function itfro(x,gd1,gd2,gd3,gd4,h1,h2,h3,h4,n,t,k)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k)
    #
    #println(k)
    fs1=acosh(tr(gd1*h1*adjoint(k))/2)
    #println(fs1)
    fs2=acosh(tr(gd2*h2*adjoint(k))/2)
    fs3=acosh(tr(gd3*h3*adjoint(k))/2)
    fs4=acosh(tr(gd4*h4*adjoint(k))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                    #println("pr1=",pr1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                    #println("pr2=",pr2)
                end
                if fs3-2*pi*im*n2 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n2)
                    #println("pr2=",pr2)
                end
                if fs4-2*pi*im*n2 == 0
                    pr4=1.0+im*0.0
                else
                    pr4=fn1(fs4,n2)
                    #println("pr2=",pr2)
                end
                result=result+pr1*pr2*pr3*pr4*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2+(fs3-2*pi*im*n1)^2+(fs4-2*pi*im*n1)^2)/t)
            #println(exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t))
        end
    end
    
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

#Integrand for 1st order Volume operator on open 4 vertex
function itfropv(x,gd1,gd2,gd3,gd4,h1,h2,h3,h4,n,t,k)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k)
    #
    #println(k)
    fs1=acosh(tr(gd1*h1*adjoint(k))/2)
    #println(fs1)
    fs2=acosh(tr(gd2*h2*adjoint(k))/2)
    fs3=acosh(tr(gd3*h3*adjoint(k))/2)
    fs4=acosh(tr(gd4*h4*adjoint(k))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                    #println("pr1=",pr1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                    #println("pr2=",pr2)
                end
                if fs3-2*pi*im*n2 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n2)
                    #println("pr2=",pr2)
                end
                if fs4-2*pi*im*n2 == 0
                    pr4=1.0+im*0.0
                else
                    pr4=fn1(fs4,n2)
                    #println("pr2=",pr2)
                end
                result=result+pr1*pr2*pr3*pr4*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2+(fs3-2*pi*im*n1)^2+(fs4-2*pi*im*n1)^2)/t)
            #println(exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t))
        end
    end
    
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

function itfro6(x,gd1,gd2,gd3,gd4,gd5,gd6,h1,h2,h3,h4,h5,h6,n,t,k)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k)
    #
    #println(k)
    fs1=acosh(tr(gd1*h1*adjoint(k))/2)
    #println(fs1)
    fs2=acosh(tr(gd2*h2*adjoint(k))/2)
    fs3=acosh(tr(gd3*h3*adjoint(k))/2)
    fs4=acosh(tr(gd4*h4*adjoint(k))/2)
    fs5=acosh(tr(gd5*h5*adjoint(k))/2)
    fs6=acosh(tr(gd6*h6*adjoint(k))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                    #println("pr1=",pr1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                    #println("pr2=",pr2)
                end
                if fs3-2*pi*im*n2 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n2)
                    #println("pr2=",pr2)
                end
                if fs4-2*pi*im*n2 == 0
                    pr4=1.0+im*0.0
                else
                    pr4=fn1(fs4,n2)
                    #println("pr2=",pr2)
                end
                if fs5-2*pi*im*n2 == 0
                    pr5=1.0+im*0.0
                else
                    pr5=fn1(fs5,n2)
                    #println("pr2=",pr2)
                end
                if fs6-2*pi*im*n2 == 0
                    pr6=1.0+im*0.0
                else
                    pr6=fn1(fs6,n2)
                    #println("pr2=",pr2)
                end
                result=result+pr1*pr2*pr3*pr4*pr5*pr6*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2+(fs3-2*pi*im*n1)^2+(fs4-2*pi*im*n1)^2+(fs5-2*pi*im*n1)^2+(fs6-2*pi*im*n1)^2)/t)
            #println(exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t))
        end
    end
    
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

#Integrand for 2 flower graph
function itfr(x,gd1,gd2,h1,h2,n,t,k)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k)
    #
    #println(k)
    fs1=acosh(tr(gd1*k*h1*adjoint(k))/2)
    #println(fs1)
    fs2=acosh(tr(gd2*k*h2*adjoint(k))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                    #println("pr1=",pr1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                    #println("pr2=",pr2)
                end
                result=result+pr1*pr2*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            #println(exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t))
        end
    end
    
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

#Integrand for the 3-Bridge graph
function it3r(x,gd1,gd2,gd3,h1,h2,h3,n,t,k1,k2)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k1)
    hpp123(x[4],x[5],x[6],k2)
    #
    #println(k)
    fs1=acosh(tr(gd1*k1*h1*adjoint(k2))/2)
    #println(fs1)
    fs2=acosh(tr(gd2*k1*h2*adjoint(k2))/2)
    fs3=acosh(tr(gd3*k1*h3*adjoint(k2))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
            for n3 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                    #println("pr1=",pr1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                    #println("pr2=",pr2)
                end
                if fs3-2*pi*im*n2 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n3)
                    #println("pr2=",pr2)
                end
                result=result+pr1*pr2*pr3*exp(((fs3-2*pi*im*n3)^2+(fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            end
            #println(exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t))
        end
    end
    
    return real((sin(x[1])^2)*sin(x[2])*(sin(x[4])^2)*sin(x[5])*result/(2*pi^2)/(2*pi^2))
end



function itfi(x)
    fs1=f1(z1,w1,x)
    fs2=f2(z2,w2,x)
    res1=0
    for n1 in 0:n
        for n2 in 0:n
            if fs1 != 0
                if fs2 != 0
                    res1=res1+(fs1-2*pi*im*n1)/sinh(fs1-2*pi*im*n1)*(fs2-2*pi*im*n2)/sinh(fs2-2*pi*im*n2)*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
                else
                    res1=res1+(fs1-2*pi*im*n1)/sinh(fs1-2*pi*im*n1)*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
                end
            else
                if fs2 != 0
                    res1=res1+(fs2-2*pi*im*n2)/sinh(fs2-2*pi*im*n2)*exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
                else
                    res1=res1+exp(((fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
                end
            end
        end
    end
    return imag(res1)
end

function fn1(f,nf)
    return (f-2*pi*im*nf)/sinh(f-2*pi*im*nf)
end

# real part of the integrand for 3 flower graph
function iter(x,gd1,gd2,gd3,h1,h2,h3,n,t,k)
    #vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    #s1=vp1[1]
    #s2=vp1[2]
    #s3=vp1[3]
    hpp123(x[1],x[2],x[3],k)
    #
    fs1=acosh(tr(gd1*k*h1*adjoint(k))/2)
    fs2=acosh(tr(gd2*k*h2*adjoint(k))/2)
    fs3=acosh(tr(gd3*k*h3*adjoint(k))/2)
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
            for n3 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                end
                if fs3-2*pi*im*n3 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n3)
                end
                result=result+pr1*pr2*pr3*exp(((fs3-2*pi*im*n3)^2+(fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            end
        end
    end
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

function iter_test(x,gd1,gd2,gd3,h1,h2,h3,n,t,vp1,k)
    vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    s1=vp1[1]
    s2=vp1[2]
    s3=vp1[3]
    hpp123(s1,s2,s3,k)
    k1=k[1,1]
    k2=k[1,2]
    k3=k[2,1]
    k4=k[2,2]
    #rs1=gd1*k*h1*inv(k)
    #rs2=gd2*k*h2*inv(k)
    #rs3=gd3*k*h3*inv(k)
    fs1=acosh((k4*(h1[1,1]*(gd1[1,1]*k1+gd1[1,2]*k3)+h1[2,1]*(gd1[1,1]k2+gd1[1,2]*k4))-k3*(h1[1,2]*(gd1[1,1]*k1+gd1[1,2]*k3)+h1[2,2]*(gd1[1,1]*k2+gd1[1,2]*k4))-k2*(h1[1,1]*(gd1[2,1]*k1+gd1[2,2]*k3)+h1[2,1]*(gd1[2,1]*k2+gd1[2,2]*k4))+k1*(h1[1,2]*(gd1[2,1]*k1+gd1[2,2]*k3)+h1[2,2]*(gd1[2,1]*k2+gd1[2,2]*k4)))/(k1*k4-k2*k3)/2)
    #
    fs2=acosh((k4*(h2[1,1]*(gd2[1,1]*k1+gd2[1,2]*k3)+h2[2,1]*(gd2[1,1]k2+gd2[1,2]*k4))-k3*(h2[1,2]*(gd2[1,1]*k1+gd2[1,2]*k3)+h2[2,2]*(gd2[1,1]*k2+gd2[1,2]*k4))-k2*(h2[1,1]*(gd2[2,1]*k1+gd2[2,2]*k3)+h2[2,1]*(gd2[2,1]*k2+gd2[2,2]*k4))+k1*(h2[1,2]*(gd2[2,1]*k1+gd2[2,2]*k3)+h2[2,2]*(gd2[2,1]*k2+gd2[2,2]*k4)))/(k1*k4-k2*k3)/2)
    #
    #
    fs3=acosh((k4*(h3[1,1]*(gd3[1,1]*k1+gd3[1,2]*k3)+h3[2,1]*(gd3[1,1]k2+gd3[1,2]*k4))-k3*(h3[1,2]*(gd3[1,1]*k1+gd3[1,2]*k3)+h3[2,2]*(gd3[1,1]*k2+gd3[1,2]*k4))-k2*(h3[1,1]*(gd3[2,1]*k1+gd3[2,2]*k3)+h3[2,1]*(gd3[2,1]*k2+gd3[2,2]*k4))+k1*(h3[1,2]*(gd3[2,1]*k1+gd3[2,2]*k3)+h3[2,2]*(gd3[2,1]*k2+gd3[2,2]*k4)))/(k1*k4-k2*k3)/2)
    #
    #
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
            for n3 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                end
                if fs3-2*pi*im*n3 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n3)
                end
                result=result+pr1*pr2*pr3*exp(((fs3-2*pi*im*n3)^2+(fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            end
        end
    end
    return real((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

# imagine part of the integrand for 3 flower graph
function itei(x,gd1,gd2,gd3,h1,h2,h3,n,t,vp1,k)
    vp1.=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    s1=vp1[1]
    s2=vp1[2]
    s3=vp1[3]
    hpp123(s1,s2,s3,k)
    k1=k[1,1]
    k2=k[1,2]
    k3=k[2,1]
    k4=k[2,2]
    #rs1=gd1*k*h1*inv(k)
    #rs2=gd2*k*h2*inv(k)
    #rs3=gd3*k*h3*inv(k)
    fs1=acosh((gd1[2]*(-k3^2*h1[2]+k4^2*h1[3]+k3*k4*(h1[1]-h1[4]))+gd1[1]*(k1*k4*h1[1]-k1*k3*h1[2]+k2*k4*h1[3]-k2*k3*h1[4]))/(-k2*k3+k1*k4)/2+(gd1[3]*k1^2*h1[2]-gd1[3]*k2^2*h1[3]-gd1[4]*k2*(k3*h1[1]+k4*h1[3])+gd1[3]*k1*k2*(-h1[1]+h1[4])+gd1[4]*k1*(k3*h1[2]+k4*h1[4]))/(-k2*k3+k1*k4)/2)
    #
    fs2=acosh((gd2[2]*(-k3^2*h2[2]+k4^2*h2[3]+k3*k4*(h2[1]-h2[4]))+gd2[1]*(k1*k4*h2[1]-k1*k3*h2[2]+k2*k4*h2[3]-k2*k3*h2[4]))/(-k2*k3+k1*k4)/2+(gd2[3]*k1^2*h2[2]-gd2[3]*k2^2*h2[3]-gd2[4]*k2*(k3*h2[1]+k4*h2[3])+gd2[3]*k1*k2*(-h2[1]+h2[4])+gd2[4]*k1*(k3*h2[2]+k4*h2[4]))/(-k2*k3+k1*k4)/2)
    #
    fs3=acosh((gd3[2]*(-k3^2*h3[2]+k4^2*h3[3]+k3*k4*(h3[1]-h3[4]))+gd3[1]*(k1*k4*h3[1]-k1*k3*h3[2]+k2*k4*h3[3]-k2*k3*h3[4]))/(-k2*k3+k1*k4)/2+(gd3[3]*k1^2*h3[2]-gd3[3]*k2^2*h3[3]-gd3[4]*k2*(k3*h3[1]+k4*h3[3])+gd3[3]*k1*k2*(-h3[1]+h3[4])+gd3[4]*k1*(k3*h3[2]+k4*h3[4]))/(-k2*k3+k1*k4)/2)
    #
    result=0.0+im*0.0
    for n1 in 0:n
        for n2 in 0:n
            for n3 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1.0+im*0.0
                else
                    pr1=fn1(fs1,n1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1.0+im*0.0
                else
                    pr2=fn1(fs2,n2)
                end
                if fs3-2*pi*im*n3 == 0
                    pr3=1.0+im*0.0
                else
                    pr3=fn1(fs3,n3)
                end
                result=result+pr1*pr2*pr3*exp(((fs3-2*pi*im*n3)^2+(fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            end
        end
    end
    return imag((sin(x[1])^2)*sin(x[2])*result/(2*pi^2))
end

#the integrand for 3 flower graph
function ite(x,gd1,gd2,gd3,h1,h2,h3,n,t)
    vp1=x[1]*[cos(x[3])*sin(x[2]) sin(x[3])*sin(x[2]) cos(x[2])]
    k=hpp(vp1)
    rs1=gd1*k*h1*inv(k)
    rs2=gd2*k*h2*inv(k)
    rs3=gd3*k*h3*inv(k)
    fs1=acosh(xl(rs1))
    fs2=acosh(xl(rs2))
    fs3=acosh(xl(rs3))
    result=0
    for n1 in 0:n
        for n2 in 0:n
            for n3 in 0:n
                if fs1-2*pi*im*n1 == 0
                    pr1=1
                else
                    pr1=fn1(fs1,n1)
                end
                if fs2-2*pi*im*n2 == 0
                    pr2=1
                else
                    pr2=fn1(fs2,n2)
                end
                if fs3-2*pi*im*n3 == 0
                    pr3=1
                else
                    pr3=fn1(fs3,n3)
                end
                result=result+pr1*pr2*pr3*exp(((fs3-2*pi*im*n3)^2+(fs2-2*pi*im*n2)^2+(fs1-2*pi*im*n1)^2)/t)
            end
        end
    end
    return result
end

#SL(2,C) element#
function slce(z)
    p1=imag(z[1])
    z1=real(z[1])
    p2=imag(z[2])
    z2=real(z[2])
    p3=imag(z[3])
    z3=real(z[3])
    tt=sqrt(p1^2+2*im*p1*z1+(p2+im*z2)^2+(p3-z1+im*z3)*(p3+z1+im*z3))
    if tt!=0
        a=cosh(tt/2)-(p3+im*z3)*sinh(tt/2)/tt
        b=-(p1-im*p2+im*z1+z2)*sinh(tt/2)/tt
        c=-(p1+im*p2+im*z1-z2)*sinh(tt/2)/tt
        d=cosh(tt/2)+(p3+im*z3)*sinh(tt/2)/tt
        return [a b;c d]
    end
    return [1 0;0 1]
end

function slce2(z)
    p1=imag(z[1])
    r1=real(z[1])
    p2=imag(z[2])
    r2=real(z[2])
    p3=imag(z[3])
    r3=real(z[3])
    tt1=sqrt(r1^2+r2^2+r3^2)
    tt2=sqrt(p1^2+p2^2+p3^2)
    if tt1!=0
        a1=cos(tt1/2)-im*r3*sin(tt1/2)/tt1
        b1=(-im*r1-r2)*sin(tt1/2)/tt1
        c1=(-im*r1+r2)*sin(tt1/2)/tt1
        d1=cos(tt1/2)+im*r3*sin(tt1/2)/tt1
        m1=[a1 b1;c1 d1]
    else
        m1=[1.0+im*0.0 0.0+im*0.0;0.0+im*0.0 1.0+im*0.0]
    end
    if tt2!=0
        a2=cosh(tt2/2)-p3*sinh(tt2/2)/tt2
        b2=-(p1-im*p2)*sinh(tt2/2)/tt2
        c2=-(p1+im*p2)*sinh(tt2/2)/tt2
        d2=cosh(tt2/2)+p3*sinh(tt2/2)/tt2
        m2=[a2 b2;c2 d2]
    else
        m2=[1.0+im*0.0 0.0+im*0.0;0.0+im*0.0 1.0+im*0.0]
    end
    #println(m1)
    #println(m2)
    return m2*m1
end
#Peakedness
function pkj(j,g)
    glc(j,g)
end

#@time (val,err) = hcubature(it2, [0 0 0], [pi pi 2*pi], reltol=1e-12, abstol=1e-12, maxevals=0)

#[val/(2*pi^2) err]

#----------------------------------------------------------------------------------------------------
#2 flower graph overlap function
#=
z1=0.6
z2=0.75
w1=0.6
w2=0.8

theta=0
chi=0

n=0
t=0.2

(val1,err1) = hcubature(itf, [0 0 0], [pi pi 2*pi], reltol=1e-12, abstol=1e-12, maxevals=0)
(val2,err2) = hcubature(itfi, [0 0 0], [pi pi 2*pi], reltol=1e-12, abstol=1e-12, maxevals=0)

println(val1)
println(val2)
=#






#=
println(cos(ttk1([0 0 0])))
println(cos(conj(ttk2([0 0 0]))))
println(cos(z1)*cos(conj(w1)))
println(sin(z1)*sin(conj(w1)))
print(cos(z1)*cos(conj(w1))+sin(z1)*sin(conj(w1)))
=#
#cos(z1)*cos(conj(w1))+sin(z1)*sin(conj(w1))*ttk1([0 0 0])





#--------------------------------------------------------------------------------------------------------
#3 flower graph normalization
#=
jp=BigFloat(5)*[1 1 1 1]

z1v=[0.0+im*0.0 0.0+im*0.0 0.0+im*jp[1]] #SL2C element for edge 1 in |psi>
z2v=[pi/2+im*0.0 0.0+im*jp[1] 0.0+im*0.0] #SL2C element for edge 2 in |psi>
z3v=[0.0+im*jp[1] pi/2+im*0.0 0.0+im*0.0] #SL2C element for edge 3 in |psi>

w1v=[0.0+im*0.0 0.0+im*0.0 0.0+im*jp[1]] #SL2C element for edge 1 in <psi|
w2v=[pi/2+im*0.0 0.0+im*jp[1] 0.0+im*0.0] #SL2C element for edge 2 in <psi|
w3v=[0.0+im*jp[1] pi/2+im*0.0 0.0+im*0.0] #SL2C element for edge 3 in <psi|

g1=slce(w1v)
gd1=conj(transpose(g1))
g2=slce(w2v)
gd2=conj(transpose(g2))
g3=slce(w3v)
gd3=conj(transpose(g3))
h1=slce(z1v)
h2=slce(z2v)
h3=slce(z3v)

n=0
t=0.5

(val1,err1) = hcubature(iter, [0 0 0], [pi pi 2*pi], reltol=1e-12, abstol=1e-12, maxevals=0)
=#

