
using WignerSymbols
    function X(b,c)
        2*b*(2*b+1)*(2*b+2)*(2*c)*(2*c+1)*(2*c+2)
    end

    function Ax(x,y)
        sqrt((2*x+1)*(2*y+1))
    end

    function SJ1(a,b,c,d,e,f)
        if abs(b-a)<=c<=abs(b+a)&&abs(c-d)<=e<=abs(c+d)&&abs(b-d)<=f<=abs(b+d)&&abs(f-a)<=e<=abs(f+a)&&a+b+c==floor(a+b+c)&&c+d+e==floor(c+d+e)&&b+d+f==floor(b+d+f)&&a+e+f==floor(a+e+f)
            return wigner6j(a,b,c,d,e,f)
        end 
        return 0
    end

    function qa(j)
        a2=Fai4(j)[:,2]
        al2=zeros((size(a2,1),size(a2,1)))
        for x1 in a2
            for x2 in a2
                i1=x1-minimum(a2)+1
                i2=x2-minimum(a2)+1
                if i2==i1-1
                    al2[Int32(i1),Int32(i2)]=1/sqrt((2*x1-1)*(2*x1+1))*sqrt(((j[1]+j[2]+x1+1)*(-j[1]+j[2]+x1))*(j[1]-j[2]+x1)*(j[1]+j[2]-x1+1)*(j[3]+j[4]+x1+1)*(-j[3]+j[4]+x1)*(j[3]-j[4]+x1)*(j[3]+j[4]-x1+1))
                end
            end
        end
        return al2-transpose(al2)
    end
    
    function qa2(j)
        a2=Fai4(j)[:,2]
        al2=zeros((size(a2,1),size(a2,1)))
        for x1 in a2
            for x2 in a2
                i1=x1-minimum(a2)+1
                i2=x2-minimum(a2)+1
                if i2==i1-1
                    al2[Int32(i1),Int32(i2)]=1/(sqrt((2*x1-1))*sqrt((2*x1+1)))*((sqrt(j[1]+j[2]+x1+1)*sqrt(-j[1]+j[2]+x1))*sqrt(j[1]-j[2]+x1)*sqrt(j[1]+j[2]-x1+1)*sqrt(j[3]+j[4]+x1+1)*sqrt(-j[3]+j[4]+x1)*sqrt(j[3]-j[4]+x1)*sqrt(j[3]+j[4]-x1+1))
                end
            end
        end
        return al2-transpose(al2)
    end

    function qas(j)
        a2=Fai4(j)[:,2]
        resVs=spzeros(size(a2,1),size(a2,1))
        for x1 in a2
            for x2 in a2
                i1=x1-minimum(a2)+1
                i2=x2-minimum(a2)+1
                if i2==i1-1
                    resVs[Int32(i1),Int32(i2)]=1/(sqrt((2*x1-1))*sqrt((2*x1+1)))*((sqrt(j[1]+j[2]+x1+1)*sqrt(-j[1]+j[2]+x1))*sqrt(j[1]-j[2]+x1)*sqrt(j[1]+j[2]-x1+1)*sqrt(j[3]+j[4]+x1+1)*sqrt(-j[3]+j[4]+x1)*sqrt(j[3]-j[4]+x1)*sqrt(j[3]+j[4]-x1+1))
                end
            end
        end
        return resVs-transpose(resVs)
    end

    function qatt(j)
        a2=Fai4(j)[:,2]
        al2=zeros((size(a2,1),size(a2,1)))
        for i1 in 1:size(a2,1)
            for i2 in 1:size(a2,1)
                x1=a2[i1]
                x2=a2[i2]
                if i2==i1-1
                    al2[Int32(i1),Int32(i2)]=1/sqrt((2*x1-1)*(2*x1+1))*sqrt(((j[1]+j[2]+x1+1)*(-j[1]+j[2]+x1))*(j[1]-j[2]+x1)*(j[1]+j[2]-x1+1)*(j[3]+j[4]+x1+1)*(-j[3]+j[4]+x1)*(j[3]-j[4]+x1)*(j[3]+j[4]-x1+1))
                end
            end
        end
        return al2-transpose(al2)
    end
    
    function ro1(ls1)
        A1=0
        t=0
        for t1 in 1:3
            for ii1 in 2:3
                if ls1[ii1-1]>ls1[ii1]
                    A1=ls1[ii1]
                    ls1[ii1]=ls1[ii1-1]
                    ls1[ii1-1]=A1
                    t=t+1
                end
            end
        end
        if rem(t,2)==0
            return 1
        else
            return -1
        end
    end
    
    function qo(j,ml,mr,I,J,K)
        d1=0
        cB1=0
        cB2=0
        cA1=0
        cA2=0
        cC1=0
        csum=0
        for mi in 1:size(ml,1)
            if mi != I&&mi != J&&mi != K
                if ml[mi]==mr[mi]
                    d1=1
                end
            end
        end
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    if i!=j&&j!=k
                        sign1=ro1([i j k])
                        I1=[I J K][i]
                        J1=[I J K][j]
                        K1=[I J K][k]
                        if ml[I1]==mr[I1]-1 && ml[J1]==mr[J1]+1
                            cB1=sqrt(j[I1]*(j[I1]+1)-ml[I1]*(ml[I1]-1))
                            println(j,ml,mr)
                            cA1=sqrt(j[J1]*(j[J1]+1)-ml[J1]*(ml[J1]+1))
                        end
                        if ml[I1]==mr[I1]+1 && ml[J1]==mr[J1]-1
                            cB2=sqrt(j[J1]*(j[J1]+1)-ml[J1]*(ml[J1]-1))
                            cA2=sqrt(j[I1]*(j[I1]+1)-ml[I1]*(ml[I1]+1))
                        end
                        if ml[K1]==mr[K1]+1
                            cC1=ml[K1]
                        end
                        csum=csum+sign1*(cB1*cA1-cA2*cB2)*cC1
                    end
                end
            end
        end
        return csum
    end

    function qo123(ls123,I,J,K)
        vres=spzeros(Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)),Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)))
        for ma in 1:Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1))
            ma1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
            ma2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
            ma3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
            ma4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
            i1=Int(ls123[1]-ma1+1)
            i2=Int(ls123[2]-ma2+1)
            i3=Int(ls123[3]-ma3+1)
            i4=Int(ls123[4]-ma4+1)
            for mb in 1:Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1))
                mb1=ls123[1]-rem(mb-1,(2*ls123[1]+1))
                mb2=ls123[2]-rem(fld(mb-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                mb3=ls123[3]-rem(fld(fld(mb-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                mb4=ls123[4]-rem(fld(fld(fld(mb-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                j1=Int(ls123[1]-mb1+1)
                j2=Int(ls123[2]-mb2+1)
                j3=Int(ls123[3]-mb3+1)
                j4=Int(ls123[4]-mb4+1)
                a11=qo(ls123,[ma1 ma2 ma3 ma4],[mb1 mb2 mb3 mb4],I,J,K)
                if a11!=0
                    vres[ma,mb]=a11
                end
            end
        end
    end
                
            

    function qvt(j,a,ap,I,J,K,n,jm1,jm2,A1,B1)
        kd1=1
        kd2=1
        for s1 in 2:I-1
            if a[s1]==ap[s1]
                kd1=kd1*1
            else
                kd1=kd1*0
            end
        end
        for s2 in K:n
            if a[s2]==ap[s2]
                kd2=kd2*1
            else
                kd2=kd2*0
            end
        end       
        if kd1*kd2 != 0
            # function pr1(n1)
            #     rp1 = Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1)
            # end
            p1 = 1/4*(-1+0*im)^(j[K]+j[I]+a[I-1]+a[K])*(-1+0*im)^(a[I]-ap[I])*(-1+0*im)^(sum(j[I+1:J-1]))*(-1+0*im)^(-sum(j[J+1:K-1]))*X(j[I],j[J])^(1/2)*X(j[J],j[K])^(1/2)*Ax(a[I],ap[I])*Ax(a[J],ap[J])
            p2 = Wj6([a[I-1] j[I] a[I] ap[I] j[I]],jm2,B1)
            if p2 != 0
                p3 = Wj6([a[K] j[K] a[K-1] ap[K-1] j[K]],jm2,B1)
                if p3 != 0
                    p4 = ((-1+0*im)^(ap[J]+ap[J-1])*Wj6([a[J] j[J] ap[J-1] a[J-1] j[J]],jm2,B1)*Wj6([ap[J-1] j[J] ap[J] a[J] j[J]],jm2,B1)-(-1+0*im)^(a[J]+a[J-1])*Wj6([ap[J] j[J] ap[J-1] a[J-1] j[J]],jm2,B1)*Wj6([a[J-1] j[J] ap[J] a[J] j[J]],jm2,B1))
                    if p4 != 0
                        p5=1
                        p6=1
                        if I+1<=J-1
                            for n1 in I+1:J-1
                            p5 = p5*(Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1))
                            end
                        else
                            p5=1
                        end
                        if J+1<=K-1
                            for n1 in J+1:K-1
                                p6 = p6*(Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1))
                            end
                        else
                            p6=1
                        end
                        return p1*p2*p3*p4*p5*p6
                    end
                end
            end     
        end
        return 0
    end

    function q1t(j,a,ap,J,K,n,jm1,jm2,A1,B1)
        kd2=1
        for s2 in K:n
            if a[s2]==ap[s2]
                kd2=kd2*1
            else
                kd2=kd2*0
            end
        end       
        if kd2 != 0
            # function pr1(n1)
            #     rp1 = Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1)
            # end
            p1 = 1/4*(-1+0*im)^(j[K]-j[1]+a[K]+1)*(-1+0*im)^(sum(j[2:J-1]))*(-1+0*im)^(-sum(j[J+1:K-1]))*X(j[1],j[J])^(1/2)*X(j[J],j[K])^(1/2)*Ax(a[J],ap[J])
            p3 = Wj6([a[K] j[K] a[K-1] ap[K-1] j[K]],jm2,B1)
            if p3 != 0
                if J>=2
                    p4 = ((-1+0*im)^(ap[J]+ap[J-1])*Wj6([a[J] j[J] ap[J-1] a[J-1] j[J]],jm2,B1)*Wj6([ap[J-1] j[J] ap[J] a[J] j[J]],jm2,B1)-(-1+0*im)^(a[J]+a[J-1])*Wj6([ap[J] j[J] ap[J-1] a[J-1] j[J]],jm2,B1)*Wj6([a[J-1] j[J] ap[J] a[J] j[J]],jm2,B1))
                else
                    p4=1
                end
                if p4 != 0
                    p5=1
                    p6=1
                    if 2<=J-1
                        for n1 in 2:J-1
                            p5 = p5*(Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1))
                        end
                    else
                        p5=1
                    end
                    if J+1<=K-1      
                        for n1 in J+1:K-1       
                            p6 = p6*(Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1))
                        end
                    else
                        p6=1
                    end
                    return p1*p3*p4*p5*p6
                end
            end    
        end
        return 0
    end

    function q12t(j,a,ap,K,n,jm1,jm2,A1,B1)
        kd2=1
        for s2 in K:n
            if a[s2]==ap[s2]
                kd2=kd2*1
            else
                kd2=kd2*0
            end
        end       
        if kd2 != 0
            #function pr1(n1)
            #    rp1 = Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*SJ1(j[n1],ap[n1-1],ap[n1],1,a[n1],a[n1-1])
            #end
            p1 = 1/2*(-1+0*im)^(j[K]-j[1]-j[2]+a[K]+1)*(-1+0*im)^(-sum(j[3:K-1]))*X(j[2],j[K])^(1/2)*Ax(a[2],ap[2])
            p2 = Wj6([j[1] j[2] a[2] ap[2] j[2]],jm1,A1)
            p3 = Wj6([a[K] j[K] a[K-1] ap[K-1] j[K]],jm2,B1)
            if p3 != 0
                p4 = a[2]*(a[2]+1)-ap[2]*(ap[2]+1)
                if p4 != 0   
                    if 3<=K-1                 
                        #p6 = prod(n1->Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*SJ1(j[n1],ap[n1-1],ap[n1],1,a[n1],a[n1-1]),3:K-1)
                        p6=1
                        for n1 in 3:K-1
                            p6=p6*Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1)
                        end
                    else
                        p6=1
                    end
                    return p1*p2*p3*p4*p6
                end
            end    
        end
        return 0
    end

    function q12tNew(j,a,ap,K,n,jm1,jm2,A1,B1)
        kd2=1
        for s2 in K:n
            if a[s2]==ap[s2]
                kd2=kd2*1
            else
                kd2=kd2*0
            end
        end       
        if kd2 != 0
            #function pr1(n1)
            #    rp1 = Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*SJ1(j[n1],ap[n1-1],ap[n1],1,a[n1],a[n1-1])
            #end
            p1 = 1/4*(-1+0*im)^(j[K]-j[1]-j[2]+a[K]+1)*(-1+0*im)^(-sum(j[3:K-1]))*X(j[2],j[K])^(1/2)*Ax(a[2],ap[2])
            p2 = Wj6([j[1] j[2] a[2] ap[2] j[2]],jm1,A1)
            p3 = Wj6([a[K] j[K] a[K-1] ap[K-1] j[K]],jm2,B1)
            if p3 != 0
                p4 = a[2]*(a[2]+1)-ap[2]*(ap[2]+1)
                if p4 != 0   
                    if 3<=K-1                 
                        #p6 = prod(n1->Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*SJ1(j[n1],ap[n1-1],ap[n1],1,a[n1],a[n1-1]),3:K-1)
                        p6=1
                        for n1 in 3:K-1
                            p6=p6*Ax(ap[n1],a[n1])*(-1+0*im)^(ap[n1-1]+a[n1-1]+1)*Wj6([j[n1] ap[n1-1] ap[n1] a[n1] a[n1-1]],jm1,A1)
                        end
                    else
                        p6=1
                    end
                    return p1*p2*p3*p4*p6
                end
            end    
        end
        return 0
    end

    function q123t(j,a,ap,n,jm1,jm2,A1,B1)
        kd2=1
        for s2 in 3:n
            if a[s2]==ap[s2]
                kd2=kd2*1
            else
                kd2=kd2*0
            end
        end  
        if kd2 != 0
            p1 = 1/2*(-1+0*im)^(-j[2]-j[1]+j[3])*X(j[2],j[3])^(1/2)*Ax(a[2],ap[2])*(a[2]*(a[2]+1)*(-1+0*im)^(a[3])-ap[2]*(ap[2]+1)*(-1+0*im)^(ap[3]))
            p2 = Wj6([j[1] j[2] a[2] ap[2] j[2]],jm1,A1)
            p3 = Wj6([a[3] j[3] a[2] ap[2] j[3]],jm2,B1)
            return p1*p2*p3
        end
        return 0
    end

    function vq1(j,n1)
        a=Fai6(j)
        l=size(a)[1]
        vol = zeros(l,l)
        for i1 in 1:l
            for j1 in 1:l
                vol[i1,j1]=q123t(j,a[i1,:],a[j1,:],n1)+q1t(j,a[i1,:],a[j1,:],3,5,n1)+q1t(j,a[i1,:],a[j1,:],5,6,n1)-q1t(j,a[i1,:],a[j1,:],2,6,n1)-qvt(j,a[i1,:],a[j1,:],2,3,4,n1)+qvt(j,a[i1,:],a[j1,:],3,4,5,n1)-qvt(j,a[i1,:],a[j1,:],4,5,6,n1)-qvt(j,a[i1,:],a[j1,:],2,4,6,n1)
            end
        end
        return vol
    end

function v123(jv1,al1,I,J,K,n,jm1,jm2,A1,B1)
    rs1=spzeros(size(al1)[1],size(al1)[1])*im
    if I==1&&J==2&&K==3
        for i in 1:size(al1)[1]
            for j in 1:size(al1)[1]
                a=q123t(jv1,al1[i,:],al1[j,:],n,jm1,jm2,A1,B1)
                if a!=0
                    rs1[i,j]=a
                end
            end
        end
    elseif I==1&&J==2
        for i in 1:size(al1)[1]
            for j in 1:size(al1)[1]
                a=-q12t(jv1,al1[i,:],al1[j,:],K,n,jm1,jm2,A1,B1)
                if a!=0
                    rs1[i,j]=a
                end
            end
        end
    elseif I==1
        for i in 1:size(al1)[1]
            for j in 1:size(al1)[1]
                a=q1t(jv1,al1[i,:],al1[j,:],J,K,n,jm1,jm2,A1,B1)
                if a!=0
                    rs1[i,j]=a
                end
            end
        end
    else    
        for i in 1:size(al1)[1]
            for j in 1:size(al1)[1]
                a=qvt(jv1,al1[i,:],al1[j,:],I,J,K,n,jm1,jm2,A1,B1)
                if a!=0
                    rs1[i,j]=a
                end
            end
        end      
    end
    return rs1
end

# function v1232(jv1,al1,I,J,K,n,jm1,jm2,A1,B1)
#     rs1=zeros(size(al1)[1],size(al1)[1])*im
#     if I==1&&J==2&&K==3
#         for i in 1:size(al1)[1]
#             for j in 1:size(al1)[1]
#                 a=q123t(jv1,al1[i,:],al1[j,:],n,jm1,jm2,A1,B1)
#                 if a!=0
#                     rs1[i,j]=a
#                 end
#             end
#         end
#     elseif I==1&&J==2
#         for i in 1:size(al1)[1]
#             for j in 1:size(al1)[1]
#                 a=-q12t(jv1,al1[i,:],al1[j,:],K,n,jm1,jm2,A1,B1)
#                 if a!=0
#                     rs1[i,j]=a
#                 end
#             end
#         end
#     elseif I==1
#         for i in 1:size(al1)[1]
#             for j in 1:size(al1)[1]
#                 a=q1t(jv1,al1[i,:],al1[j,:],J,K,n,jm1,jm2,A1,B1)
#                 if a!=0
#                     rs1[i,j]=a
#                 end
#             end
#         end
#     else    
#         for i in 1:size(al1)[1]
#             for j in 1:size(al1)[1]
#                 a=qvt(jv1,al1[i,:],al1[j,:],I,J,K,n,jm1,jm2,A1,B1)
#                 if a!=0
#                     rs1[i,j]=a
#                 end
#             end
#         end      
#     end
#     return rs1
# end

function v1232(jv1,al1,I,J,K,n,jm1,jm2,A1,B1)
    rs1=spzeros(size(al1)[1],size(al1)[1])*im
        for i in 1:size(al1)[1]
            for j in 1:size(al1)[1]
                if I==1&&J==2&&K==3
                    a=q123t(jv1,al1[i,:],al1[j,:],n,jm1,jm2,A1,B1)
                elseif I==1&&J==2
                    a=-q12t(jv1,al1[i,:],al1[j,:],K,n,jm1,jm2,A1,B1)
                elseif I==1
                    a=q1t(jv1,al1[i,:],al1[j,:],J,K,n,jm1,jm2,A1,B1)
                else
                    a=qvt(jv1,al1[i,:],al1[j,:],I,J,K,n,jm1,jm2,A1,B1)
                end
                if a!=0
                    rs1[i,j]=a
                end
            end
        end
    # return 0
    return rs1
end