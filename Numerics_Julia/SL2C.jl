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

    function MatrixExp(l)
        a=BigFloat(l[1,1])
        b=BigFloat(l[1,2])
        c=BigFloat(l[2,1])
        d=BigFloat(l[2,2])
        dt=sqrt((a-d)^2+4*b*c)
        m11=exp((a+d)/2)*(dt*cosh(dt/2)+(a-d)*sinh(dt/2))
        m12=2*b*exp((a+d)/2)*sinh(dt/2)
        m21=2*c*exp((a+d)/2)*sinh(dt/2)
        m22=exp((a+d)/2)(dt*cosh(dt/2)+(d-a)*sinh(dt/2))
        if dt == 0
            return exp((a+d)/2)*[1+(a-d)/2 b;c 1-(a-d)/2]
        else
            return [m11 m12;m21 m22]/dt
        end
    end

    function pv(j,A)
        [j[1]*[-sin(A)/2,-sqrt(3)sin(A)/2,cos(A)] j[2]*[sin(A),0,cos(A)] j[3]*[0,0,1]; j[4]*[-sin(A)/2,sqrt(3)sin(A)/2,cos(A)]]
    end

    function rlist(p)
        [cos(p[2]/2)*exp(im*(p[1]+p[3])/2) im*sin(p[2]/2)*exp(-im*(p[1]-p[3])/2)]
    end

    function glc(j,g)
        a=floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).+floatabc.(zeros((Int32(2j+1),Int32(2j+1)))).*im
        for x1 in -j:j
            for x2 in -j:j
                for l in -2j:2j
                    if j+x1-l>=0&&j-x2-l>=0&&-x1+x2+l>=0&&l>=0
                        a[Int32(j+x1+1),Int32(j+x2+1)]=a[Int32(j+x1+1),Int32(j+x2+1)]+(sqrt(pr(j+x1))*sqrt(pr(j-x1))*sqrt(pr(j+x2))*sqrt(pr(j-x2)))/(pr(j+x1-l)*pr(j-x2-l)*pr(-x1+x2+l)*pr(l))*(g[1,1]^(j-x2-l))*(g[2,2]^(j+x1-l))*(g[1,2]^(-x1+x2+l))*(g[2,1]^(l))
                    end
                end
            end
        end
        return a
    end

    function g1(j,p,u)
        a1=(rlist(u))[1]
        b1=(rlist(u))[2]
        p1=p[1]
        p2=p[2]
        p3=p[3]
        u1=[[a1 -conj(b1)]; [b1 conj(a1)]]
        gslc=zeros((2,2))+im*zeros((2,2))
        if p1^2+p2^2+p3^2==0
            u2=[[1 0;0 1]]
        else
            u2=[[cosh(1/2*sqrt(p1^2+p2^2+p3^2))-p3*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/(sqrt(p1^2+p2^2+p3^2)) -(p1-im*p2)*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/sqrt(p1^2+p2^2+p3^2)]; [-(p1+im*p2)*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/sqrt(p1^2+p2^2+p3^2) cosh(1/2*sqrt(p1^2+p2^2+p3^2))+p3*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/(sqrt(p1^2+p2^2+p3^2))]]
        end
        a1=u2*u1
        for x1 in 1:2
            for x2 in 1:2
                t1=BigFloat(real(a1[x1,x2]))
                t2=BigFloat(imag(a1[x1,x2]))
                gslc[x1,x2]=t1+im*t2
            end
        end
        return glc(j,gslc)
    end


    #SL2C element with the definition: 
    function g1(j,z)
        a1=(rlist(u))[1]
        b1=(rlist(u))[2]
        p1=imag(z[1])
        p2=imag(z[2])
        p3=imag(z[3])
        u1=[[a1 -conj(b1)]; [b1 conj(a1)]]
        gslc=zeros((2,2))+im*zeros((2,2))
        if p1^2+p2^2+p3^2==0
            u2=[[1 0;0 1]]
        else
            u2=[[cosh(1/2*sqrt(p1^2+p2^2+p3^2))-p3*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/(sqrt(p1^2+p2^2+p3^2)) -(p1-im*p2)*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/sqrt(p1^2+p2^2+p3^2)]; [-(p1+im*p2)*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/sqrt(p1^2+p2^2+p3^2) cosh(1/2*sqrt(p1^2+p2^2+p3^2))+p3*sinh(1/2*sqrt(p1^2+p2^2+p3^2))/(sqrt(p1^2+p2^2+p3^2))]]
        end
        a1=u2*u1
        for x1 in 1:2
            for x2 in 1:2
                #t1=BigFloat(real(a1[x1,x2]))
                #t2=BigFloat(imag(a1[x1,x2]))
                #gslc[x1,x2]=t1+im*t
            end
        end
        return glc(j,a1)
    end


    