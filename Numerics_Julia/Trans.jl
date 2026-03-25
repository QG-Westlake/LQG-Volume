    function lb1(j)
        j*(j+1)/2
    end

    function pv(j,A)
        [j[1]*[-sin(A)/2,-sqrt(3)sin(A)/2,cos(A)] j[2]*[sin(A),0,cos(A)] j[3]*[0,0,1] j[4]*[-sin(A)/2,sqrt(3)sin(A)/2,cos(A)]]
    end

    function pv2(j,A1,A2)
        sa=(1+cos(A2))*(cos(A1)+cos(A2))/(-1+cos(A1))
        sb=(1-(1+cos(A1)+cos(A2))^2)*(-cot(A1/2)^2)/(cos(A1)^2+2*cos(A1)*(1+cos(A2))+cos(A2)*(2+cos(A2)))
        sc=2*(1+cos(A2))*(cos(A1)+cos(A2))/(-1+cos(A1))
        if sa>=0&&sb>=0&&sc>=0
            h1=j[1]*[[1,0,0] [cos(A1),abs(sin(A1)),0] [cos(A2),(abs(sin(A1))*(1+cos(A2)))/(-1+cos(A1)),-sqrt(2)*sqrt(sa)] [-1-cos(A1)-cos(A2),-abs(cos(A1)+cos(A2))*sqrt(sb),sqrt(sc)]]
            if abs(sum(h1[:,1].+h1[:,2].+h1[:,3].+h1[:,4]))<=1e-10
                return (h1,abs(sum(h1[:,1].+h1[:,2].+h1[:,3].+h1[:,4])))
            end
        end
        return "skip"
    end

    #Randomize the entries of a list
    function randl(ls1)
        l=size(ls1)[1]
        for i in 1:l
            li=Int(floor(rand()*(l-i+1)))+1
            a=ls1[li]
            splice!(ls1, li)
            append!(ls1, a)
        end
    end

    #generating j sequence for the 3 flower graph
    function seq1(js1,js2,js3,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    ls1[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1,:]=[i1/2 j1/2 k1/2 i1/2 j1/2 k1/2]
                    ls2[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1]=i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1
                end
            end
        end
        randl(ls2)
    end

    #generating j sequence for open 4 vertex
    function seq3o(js1,js2,js3,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                        ls1[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1,:]=[i1/2 j1/2 k1/2]
                        ls2[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1]=i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1
                end
            end
        end
        randl(ls2)
    end

    #generating j sequence for open 4 vertex
    function seq4o(js1,js2,js3,js4,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:js4
                        ls1[i1*(js2+1)*(js3+1)*(js4+1)+j1*(js3+1)*(js4+1)+k1*(js4+1)+l1+1,:]=[i1/2 j1/2 k1/2 l1/2]
                        ls2[i1*(js2+1)*(js3+1)*(js4+1)+j1*(js3+1)*(js4+1)+k1*(js4+1)+l1+1]=i1*(js2+1)*(js3+1)*(js4+1)+j1*(js3+1)*(js4+1)+k1*(js4+1)+l1+1
                    end
                end
            end
        end
        randl(ls2)
    end

    #generating j sequence for open 5 vertex
    function seq5o(js1,js2,js3,js4,js5,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:js4
                        for ml1 in 0:js5 
                            ls1[i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)+j1*(js3+1)*(js4+1)*(js5+1)+k1*(js4+1)*(js5+1)+l1*(js5+1)+ml1+1,:]=[i1/2 j1/2 k1/2 l1/2 ml1/2]
                            ls2[i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)+j1*(js3+1)*(js4+1)*(js5+1)+k1*(js4+1)*(js5+1)+l1*(js5+1)+ml1+1]=i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)+j1*(js3+1)*(js4+1)*(js5+1)+k1*(js4+1)*(js5+1)+l1*(js5+1)+ml1+1
                        end
                    end
                end
            end
        end
        randl(ls2)
    end

    #generating j sequence for open 6 vertex
    function seq6o(js1,js2,js3,js4,js5,js6,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:js4
                        for m1 in 0:js5
                            for n1 in 0:js6
                                ls1[i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1)+j1*(js3+1)*(js4+1)*(js5+1)*(js6+1)+k1*(js4+1)*(js5+1)*(js6+1)+l1*(js5+1)*(js6+1)+m1*(js6+1)+n1+1,:]=[i1/2 j1/2 k1/2 l1/2 m1/2 n1/2]
                                ls2[i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1)+j1*(js3+1)*(js4+1)*(js5+1)*(js6+1)+k1*(js4+1)*(js5+1)*(js6+1)+l1*(js5+1)*(js6+1)+m1*(js6+1)+n1+1]=i1*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1)+j1*(js3+1)*(js4+1)*(js5+1)*(js6+1)+k1*(js4+1)*(js5+1)*(js6+1)+l1*(js5+1)*(js6+1)+m1*(js6+1)+n1+1
                            end
                        end
                    end
                end
            end
        end
        randl(ls2)
    end

    function seq6oMod(js1,js2,js3,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    ls1[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1,:]=[i1/2 j1/2 k1/2 i1/2 j1/2 k1/2]
                    ls2[i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1]=i1*(js2+1)*(js3+1)+j1*(js3+1)+k1+1
                end
            end
        end
        randl(ls2)
    end

    function seq7oMod(js1,js2,js3,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:(js1+js2+js3+js1+js2+js3)
                        ls1[i1*(js2+1)*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+j1*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+k1*(js1+js2+js3+js1+js2+js3+1)+l1+1,:]=[i1/2 j1/2 k1/2 i1/2 j1/2 k1/2 l1/2]
                        ls2[i1*(js2+1)*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+j1*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+k1*(js1+js2+js3+js1+js2+js3+1)+l1+1]=i1*(js2+1)*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+j1*(js3+1)*(js1+js2+js3+js1+js2+js3+1)+k1*(js1+js2+js3+js1+js2+js3+1)+l1+1
                    end
                end
            end
        end
        randl(ls2)
    end

    function seq7oModN(js1,js2,js3,ls1,ls2)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:(js1+js2+js3+js1+js2+js3)
                        if l1 <=(i1+j1+k1+i1+j1+k1)
                            append!(ls1,[[i1/2 j1/2 k1/2 i1/2 j1/2 k1/2 l1/2]])
                        end
                    end
                end
            end
        end
        for ii in 1:size(ls1,1)
            append!(ls2,ii)
        end
        randl(ls2)
    end


    function seq4o2(js1,js2,js3,js4,ls1)
        for i1 in 0:js1
            for j1 in 0:js2
                for k1 in 0:js3
                    for l1 in 0:js4
                        ls1[i1*(js2+1)*(js3+1)*(js4+1)+j1*(js3+1)*(js4+1)+k1*(js4+1)+l1+1,:]=[i1/2 j1/2 k1/2 l1/2]
                    end
                end
            end
        end
        ls1=seqm(ls1)
        ls2=[i for i in 1:size(ls1,1)]
        randl(ls2)
        return (ls1,ls2)
    end

    function seqm(list1)
        s=1
        while s<=size(list1,1)
            if size(Fai4(list1[s,:]),1)==0
                list1=rm2!(list1,s)
            s=s-1
            end
        s=s+1
        end
        return list1
    end

    function rm2!(list1,i)
        a=zeros(size(list1,1)-1,size(list1,2))
        for mi in 1:i-1
            a[mi,:].=list1[mi,:]
        end
        for mi in i+1:size(list1,1)
            a[mi-1,:]=list1[mi,:]
        end
        list1=zeros(size(list1,1)-1,size(list1,2))
        list1=a
    end

    function selector(M)
        Mr=real.(M)
        Mi=imag.(M)
        Mra=zeros(size(Mr,1),size(Mr,2))
        Mrb=zeros(size(Mr,1),size(Mr,2))
        Mia=zeros(size(Mr,1),size(Mr,2))
        Mib=zeros(size(Mr,1),size(Mr,2))
        a1=0
        a2=0
        for i in 1:size(Mr,1)
            for j in 1:size(Mr,2)
                if Mr[i,j] != 0
                    Mra[i,j]=Mr[i,j]/abs(Mr[i,j])
                    Mrb[i,j]=log10(abs(Mr[i,j]))
                    a1=Mrb[i,j]
                else
                    Mra[i,j]=0
                    Mrb[i,j]=a1
                end
                if Mi[i,j] != 0
                    Mia[i,j]=Mi[i,j]/abs(Mi[i,j])
                    Mib[i,j]=log10(abs(Mi[i,j]))
                    a2=Mib[i,j]
                else
                    Mia[i,j]=0
                    Mib[i,j]=a2
                end
            end
        end
        return (Mra,Mrb,Mia,Mib)
    end

    function svdm123(g)
        (g1,g2,g3,g4)=selector(g)
        (u1,s1,v1)=svd(g2)
        (u2,s2,v2)=svd(g4)
        for i in 1:size(g2,1)
            ar=exp10.(u1[:,1:i]*Diagonal(s1[1:i])*transpose(v1)[1:i,:]).*g1
            as1=abs(sum(real.(g).-ar))
            print((as1/size(g2,1)/size(g2,1),abs(sum(real.(g)))/size(g2,1)/size(g2,1)))
            if as1/size(g2,1)/size(g2,1) <= abs(sum(real.(g)))/size(g2,1)/size(g2,1)*10e-5
                for j in 1:size(g2,1)
                    ai=exp10.(u2[:,1:i]*Diagonal(s2[1:i])*transpose(v2)[1:i,:]).*g3
                    as2=abs(sum(imag.(g).-ai))
                    if as2/size(g2,1)/size(g2,1) <= abs(sum(imag.(g)))/size(g2,1)/size(g2,1)*10e-5
                        return (u1[:,1:i],s1[1:i],transpose(v1)[1:i,:],u2[:,1:j],s2[1:j],transpose(v2)[1:j,:],g1,g3)
                    end
                end
            end
        end
        return (u1,s1,v1,u2,s2,v2,g1,g3)
    end

    function selector1(M1)
        Ma=abs.(M1)
        Maa=zeros(size(M1,1),size(M1,2))
        for i in 1:size(Ma,1)
            for j in 1:size(Ma,2)
                if Ma[i,j] != 0
                    Maa[i,j]=1
                else
                    Maa[i,j]=0
                end
            end
        end
        return Maa
    end

    function svdm1(g,n)
        (u1,s1,v1)=svd(log.(g))
        g0=selector1(g)
        for i in 1:size(s1,1)
            ar=exp.(u1[:,1:i]*Diagonal(s1[1:i])*adjoint(v1)[1:i,:])
            #as1=sum((abs.(g.-ar)./abs.(g)).*g0)/size(g,1)/size(g,2)
            as2=maximum((abs.(g.-ar)./abs.(g)).*g0)
            #print((as1,as2)," ")
            if as2 <= 10^(-n)
                return (u1[:,1:i],s1[1:i],transpose(v1)[1:i,:])
            end
        end
        return (u1,s1,transpose(v1))
    end
            
        
        

#     function seqd(list1,core)
#         list2=list1
#         for i in 1:core
#             if i <= mod(size(list1,1),core)
#                 for j in 1:Int(ceil(size(list1,1)/core))
#                     list2[(i-1)*ceil(size(list1,1)/core)+j,:].=list1[ceil(size(list1,1)/core)*j+1]
#                 end
#             if mod(i,2)==1
#             end
#             else
#                 numcore=mod(size(list1,1),core)
#                 num1=Int(ceil(size(list1,1)/core))*numcore
#                 realcore=i-numcore
#                 for j in 1:Int(floor(size(list1,1)/core))
                    
#                 end
#             end
#         end
#     end
                

    # 4-D intertwiner i_{m1 m2 m3 m4}
    function Inter4(ls1, ls2, Atl, jml, a, s1, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi)
        Wj1_mod(ls1,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        @inbounds ai.=broadcast(^,-1,Int.(ls1[1:a,3].-s1)).*rs2
        Wj1_mod(ls2,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        @inbounds bi.=rs2
        @inbounds ai.=ai.*bi
        return 0.0
    end

    function Inter4_Compact(ls1, ls2, Atl, jml, a, s1, ij1, ai)
#         println(ls1,ls2,a,s1)
        for all1 in 1:a
            ai[all1:all1].=Wj1(ls1[all1,:],jml,Atl,ij1)*(-1)^(Int(ls1[all1,3]-s1))*Wj1(ls2[all1,:],jml,Atl,ij1)
        end
        return 0.0
    end

    function Inter5_Compact(ls1, ls2, ls3, Atl, jml, a, s1, s2, ij1, ai)
#         println(ls1,ls2,a,s1)
        for all1 in 1:a
            ai[all1:all1].=Wj1(ls1[all1,:],jml,Atl,ij1)*(-1)^(Int(ls1[all1,3]-s1))*Wj1(ls2[all1,:],jml,Atl,ij1)*(-1)^(Int(ls2[all1,3]-s2))*Wj1(ls3[all1,:],jml,Atl,ij1)
        end
        return 0.0
    end

    function Inter6_Compact(ls1, ls2, ls3, ls4, Atl, jml, a, s1, s2, s3, ij1, ai)
        #         println(ls1,ls2,a,s1)
                for all1 in 1:a
                    ai[all1:all1].=Wj1(ls1[all1,:],jml,Atl,ij1)*(-1)^(Int(ls1[all1,3]-s1))*Wj1(ls2[all1,:],jml,Atl,ij1)*(-1)^(Int(ls2[all1,3]-s2))*Wj1(ls3[all1,:],jml,Atl,ij1)*(-1)^(Int(ls3[all1,3]-s3))*Wj1(ls4[all1,:],jml,Atl,ij1)
                end
                return 0.0
            end

    # 6-D intertwiner
    function cond1(j1,j2,j3,m1,m2,m3)
        if abs(m1)<=j1&&abs(m2)<=j2&&abs(m3)<=j3&&rem(2*(j1+j1+j1),2)==0&&m1+m2==-m3
            return true
        end
        return false
    end

    function Inter6_mod(ja, ia, m1, m2, m3, m4, m5, threejlist, jml, wlist)
        s1 = -m1-m2
        s2 = s1-m3
        s3 = s2-m4
        m6 = s3-m5
        if abs(m6)<=ja[6]
            ai=(-1)^(ia[2]-s1)*Wj1([ja[1] ja[2] ia[2] m1 m2 s1],jml,threejlist,wlist)
            if ai != 0
                bi=(-1)^(ia[3]-s2)*Wj1([ia[2] ja[3] ia[3] -s1 m3 s2],jml,threejlist,wlist)
                if bi != 0
                    ci=(-1)^(ia[4]-s3)*Wj1([ia[3] ja[4] ia[4] -s2 m4 s3],jml,threejlist,wlist)
                    if ci != 0
                        di=Wj1([ia[4] ja[5] ja[6] -s3 m5 m6],jml,threejlist,wlist)
                        if di !=0
                            return ai*bi*ci*di
                        end
                    end
                end
            end
        end
        return 0.0
    end

    function Inter6_mod2(ls1, ls2, ls3, ls4, Atl, jml, a, s1, s2, s3, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi, ci, di)
        Wj1_mod(ls1,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ai.=broadcast(^,-1,(@view ls1[1:a,3]).-s1).*rs2
        Wj1_mod(ls2,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        bi.=broadcast(^,-1,(@view ls2[1:a,3]).-s2).*rs2
        Wj1_mod(ls3,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ci.=broadcast(^,-1,(@view ls3[1:a,1]).-s3).*broadcast(^,-1,((@view ls3[1:a,1]).+(@view ls3[1:a,2]).+(@view ls3[1:a,3]))).*rs2
        Wj1_mod(ls4,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        di.=rs2
        ai.=(((ai.*bi).*ci).*di)
        return 0.0
    end

    function Inter7(ls1, ls2, ls3, ls4, ls5, Atl, jml, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi, ci, di, ei)
        Wj1_mod(ls1,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ai.=broadcast(^,-1,(@view ls1[1:a,3]).-s1).*rs2
        Wj1_mod(ls2,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        bi.=broadcast(^,-1,(@view ls2[1:a,3]).-s2).*rs2
        Wj1_mod(ls3,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ci.=broadcast(^,-1,(@view ls3[1:a,1]).-s3).*broadcast(^,-1,((@view ls3[1:a,1]).+(@view ls3[1:a,2]).+(@view ls3[1:a,3]))).*rs2
        Wj1_mod(ls4,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        di.=broadcast(^,-1,(@view ls4[1:a,1]).-s4).*broadcast(^,-1,((@view ls4[1:a,1]).+(@view ls4[1:a,2]).+(@view ls4[1:a,3]))).*rs2
        Wj1_mod(ls5,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ei.=rs2
        ai.=((((ai.*bi).*ci).*di).*ei)
        return 0.0
    end

    function Inter7NewBackup(ls1, ls2, ls3, ls4, ls5, Atl, jml, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi, ci, di, ei)
        Wj1_mod(ls1,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ai.=broadcast(^,-1,(@view ls1[1:a,1]).-(@view ls1[1:a,2]).+s1).*rs2
        Wj1_mod(ls2,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        bi.=broadcast(^,-1,(@view ls2[1:a,1]).-(@view ls2[1:a,2]).+s2).*rs2
        Wj1_mod(ls3,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ci.=broadcast(^,-1,(@view ls3[1:a,3]).-(@view ls3[1:a,2]).+s3).*broadcast(^,-1,((@view ls3[1:a,1]).+(@view ls3[1:a,2]).+(@view ls3[1:a,3]))).*rs2
        Wj1_mod(ls4,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        di.=broadcast(^,-1,(@view ls4[1:a,3]).-(@view ls4[1:a,2]).+s4).*broadcast(^,-1,((@view ls4[1:a,1]).+(@view ls4[1:a,2]).+(@view ls4[1:a,3]))).*rs2
        Wj1_mod(ls5,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ei.=broadcast(^,-1,(@view ls5[1:a,3]).-(@view ls5[1:a,1]).-(@view ls5[1:a,5])).*rs2
        ai.=((((ai.*bi).*ci).*di).*ei)
        return 0.0
    end

    function Inter7New(ls1, ls2, ls3, ls4, ls5, Atl, jml, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi, ci, di, ei)
        Wj1_mod(ls1,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ai.=broadcast(^,-1,(@view ls1[1:a,3]).-s1).*rs2
        Wj1_mod(ls2,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        bi.=broadcast(^,-1,(@view ls2[1:a,3]).-s2).*rs2
        Wj1_mod(ls3,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ci.=broadcast(^,-1,(@view ls3[1:a,1]).-s3).*broadcast(^,-1,((@view ls3[1:a,1]).+(@view ls3[1:a,2]).+(@view ls3[1:a,3]))).*rs2
        Wj1_mod(ls4,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        di.=broadcast(^,-1,(@view ls4[1:a,1]).-s4).*broadcast(^,-1,((@view ls4[1:a,1]).+(@view ls4[1:a,2]).+(@view ls4[1:a,3]))).*rs2
        Wj1_mod(ls5,jml,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        ei.=rs2
        ai.=((((ai.*bi).*ci).*di).*ei)
        return 0.0
    end

    function Inter6_mod3(ls1, ls2, ls3, ls4, Atl, jml, a, s1, s2, s3, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, ai, bi, ci, di)
        ai.=broadcast(^,-1,(@view ls1[1:a,3]).-s1).*broadcast((x1,x2,x3,y1,y2,y3)->tjst(x1,x2,x3,y1,y2,y3),ls1[1:a,1],ls1[1:a,2],ls1[1:a,3],ls1[1:a,4],ls1[1:a,5],ls1[1:a,6])
        bi.=broadcast(^,-1,(@view ls2[1:a,3]).-s2).*broadcast((x1,x2,x3,y1,y2,y3)->tjst(x1,x2,x3,y1,y2,y3),ls2[1:a,1],ls2[1:a,2],ls2[1:a,3],ls2[1:a,4],ls2[1:a,5],ls2[1:a,6])
        ci.=broadcast(^,-1,(@view ls3[1:a,1]).-s3).*broadcast((x1,x2,x3,y1,y2,y3)->tjst(x1,x2,x3,y1,y2,y3),ls3[1:a,1],ls3[1:a,2],ls3[1:a,3],ls3[1:a,4],ls3[1:a,5],ls3[1:a,6])
        di.=broadcast((x1,x2,x3,y1,y2,y3)->tjst(x1,x2,x3,y1,y2,y3),ls4[1:a,1],ls4[1:a,2],ls4[1:a,3],ls4[1:a,4],ls4[1:a,5],ls4[1:a,6])
        ai.=(((ai.*bi).*ci).*di)
        return 0.0
    end

    # 4-D recoupling sequence (j1=a1, a2, a3=j4, a4=0)
    function Fai4(j)
        l1 = [0 0 0 0];
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            if abs(i1-j[3])<=j[4]<=abs(i1+j[3])&&rem(2*(j[2]+j[1]),2)==rem(2*(j[3]+j[4]),2)
                l1=vcat(l1, [j[1] i1 j[4] 0])
            end
        end
        if size(l1)[1]>1
            return l1[2:end,:]
        else
            return []
        end
    end

    function Fai5(j)
        l1 = [0 0 0 0 0];
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            for j1 in abs(j[3]-i1):abs(j[3]+i1)
                if abs(j1-j[4])<=j[5]<=abs(j1+j[4]) && mod((j[5]+j[4]+j1)*2,2)==0
                    l1=vcat(l1, [j[1] i1 j1 j[5] 0])
                end
            end
        end
        return l1[2:end,:]
    end

    function Fai4a(j)
        l1 = [0 0 0 0];
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            for j1 in abs(j[3]-i1):abs(j[3]+i1)
                for k1 in abs(j[4]-j1):abs(j[4]+j1)
                    l1=vcat(l1, [j[1] i1 j1 k1])
                end
            end
        end
        return l1[2:end,:]
    end

    # 6-D recoupling sequence (j1=a1, a2, a3, a4, a5=j6, a6=0)
    function Fai6(j)
        l1 = [0 0 0 0 0 0];
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            for j1 in abs(j[3]-i1):abs(j[3]+i1)
                for k1 in abs(j[4]-j1):abs(j[4]+j1)
                    if abs(k1-j[5])<=j[6]<=abs(k1+j[5]) && mod((j[6]+j[5]+k1)*2,2)==0
                        l1=vcat(l1, [j[1] i1 j1 k1 j[6] 0])
                    end
                end
            end
        end
        return l1[2:end,:]
    end

    function Fai7(j)
        l1 = [0 0 0 0 0 0 0];
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            for j1 in abs(j[3]-i1):abs(j[3]+i1)
                for k1 in abs(j[4]-j1):abs(j[4]+j1)
                    for kk1 in abs(j[5]-k1):abs(j[5]+k1)
                        if abs(kk1-j[6])<=j[7]<=abs(kk1+j[6]) && mod((j[7]+j[6]+kk1)*2,2)==0
                            l1=vcat(l1, [j[1] i1 j1 k1 kk1 j[7] 0])
                        end
                    end
                end
            end
        end
        return l1[2:end,:]
    end

    function Fain(j,l,n)
        if size(j,1)>=3
        end
        for i1 in abs(j[2]-j[1]):abs(j[2]+j[1])
            for j1 in abs(j[3]-i1):abs(j[3]+i1)
                for k1 in abs(j[4]-j1):abs(j[4]+j1)
                    if abs(k1-j[5])<=j[6]<=abs(k1+j[5]) && mod((j[6]+j[5]+k1)*2,2)==0
                        l1=vcat(l1, [j[1] i1 j1 k1 j[6] 0])
                    end
                end
            end
        end
        return l1[2:end,:]
    end
    
    function TransformM_arbitrary(nl,j,t,g,jml,aj1,ij1)
        1
    end
    
    # 4-dimensional transformation matrix
function TransformM4(j,t,g01,g02,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai4(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        as=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g02))
        #println(g2c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        a=size(A1,1)
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[3]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        for ma in 1:(2*j1+1)*(2*j2+1)*(2*j3+1)
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = m1-m2
            m4 = s1+m3
            if abs(m4)<=j4
            #Block Two    
                ls1[1:a,4].=-m1
                ls1[1:a,5].=m2
                ls1[1:a,6].=s1
                ls2[1:a,4].=-s1
                ls2[1:a,5].=-m3  
                ls2[1:a,6].=m4
            #Block Three
                #println(g2c)
                i3=(-1)^(j1-m1)*(-1)^(j3-m3)*g1c[Int32(j1+m1+1),Int32(j2+m2+1)]*g2c[Int32(j3+m3+1),Int32(j4+m4+1)]
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as.=(i2.*i3).*i1
                        #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                        #print([float(i1) float(i2) float(i3)],"\n")
                        #println()
                #println([m1 m2 m3])
                #println(i1)
                #println(i2)
                #println(i3)
                Alist=(Alist.+as)
                count1=count1+1
                if m3 != m3p
                    #print(m5p)
                end
                m3p=m3
            end
        end
        #println(count1)
        return Alist
    end

function TransformM1(j,t,g01)
        g1c=sqrt(2*j+1)*exp(-t*(lb1(j)))*glc(j,g01)
    end

function TransformM3(j,t,g01,g02,g03,jml,aj1,ij1)
        Alist=Array{ComplexF64,1}(undef,Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1)))
        for i in 1:size(Alist)[1]
            Alist[i]=0
        end
        g1c=conj(sqrt(2*j[1]+1)*exp(-t*(lb1(j[1])))*glc(j[1],g01))
        g2c=conj(sqrt(2*j[2]+1)*exp(-t*(lb1(j[2])))*glc(j[2],g02))
        g3c=conj(sqrt(2*j[3]+1)*exp(-t*(lb1(j[3])))*glc(j[3],g03))
        j12=Int((2*j[1]+1)*(2*j[2]+1))
        j123=Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        for m1 in -j[1]:j[1]
            for  m2 in -j[2]:j[2]
                if abs(m1+m2)<=j[3]
                    coef1=(-1)^(Int(j[1]+j[2]+j[3]-m1-m2-(-m1-m2)))
                    @inbounds kron!(kronC1,@view(g1c[Int(j[1]+m1+1),:]),@view(g2c[Int32(j[2]+m2+1),:]))
                    @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j[3]-m1-m2+1),:]))
                    Alist.=Alist.+kronC2*Wj1([j[1] j[2] j[3] m1 m2 -m1-m2],jml,aj1,ij1)*coef1
                end
            end
        end
        return Alist
    end

function TransformM34(j,t,g01,g02,g03,g04,jml,aj1,ij1)
        al1=Fai4(j)
        Alist=zeros(size(al1,1),Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1)*(2*j[4]+1)))*im
        for i in 1:size(Alist)[1]
            Alist[i]=0
        end
        g1c=adjoint(sqrt(2*j[1]+1)*exp(-t*(lb1(j[1])))*glc(j[1],g01))
        g2c=adjoint(sqrt(2*j[2]+1)*exp(-t*(lb1(j[2])))*glc(j[2],g02))
        g3c=adjoint(sqrt(2*j[3]+1)*exp(-t*(lb1(j[3])))*glc(j[3],g03))
        g4c=adjoint(sqrt(2*j[4]+1)*exp(-t*(lb1(j[4])))*glc(j[4],g03))
        j12=Int((2*j[1]+1)*(2*j[2]+1))
        j123=Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1))
        j1234=Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1)*(2*j[4]+1))
        kronC1=zeros(j12)*im
        kronC2=zeros(j123)*im
        kronC3=zeros(j1234)*im
        for ma in 1:Int((2*j[1]+1)*(2*j[2]+1)*(2*j[3]+1))
            #Block One
            m1=j[1]-rem(ma-1,(2*j[1]+1))
            m2=j[2]-rem(fld(ma-1,(2*j[1]+1)),(2*j[2]+1))
            m3=j[3]-rem(fld(fld(ma-1,(2*j[1]+1)),(2*j[2]+1)),(2*j[3]+1))
                    if abs(m1+m2+m3)<=j[4]
                        coef1=(-1)^(Int(j[1]+j[2]+j[3]+j[4]-m1-m2-m3-(-m1-m2-m3)))
                        @inbounds kron!(kronC1,@view(g1c[Int(j[1]+m1+1),:]),@view(g2c[Int32(j[2]+m2+1),:]))
                        @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j[3]+m3+1),:]))
                        @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j[4]-m1-m2-m3+1),:]))
                        for i in 1:size(al1,1)
                            Alist[i,:].=Alist[i,:].+kronC3*Wj1([j[1] j[2] al1[i,2] m1 m2 -m1-m2],jml,aj1,ij1)*Wj1([al1[i,2] j[3] j[4] m1+m2 m3 -m1-m2-m3],jml,aj1,ij1)*coef1*(-1)^(Int(al1[i,2]+m1+m2))*sqrt(2*al1[i,2]+1)
                        end
#                         println("f1=",[j[1] j[2] j[3] j[4] m1 m2 m3 -(m1+m2+m3)],Alist[1,1])
#                         println("g=",g1c[Int(j[1]+m1+1),1],g2c[Int(j[2]+m2+1),1],g3c[Int(j[3]+m3+1),1],g4c[Int(j[4]-m1-m2-m3+1),1])
                    end
        end
        return Alist
    end


function TransformM3tt(j,t,g01,g02,g03,jml,aj1,ij1,m3,m4,m5)
        g11=conj(glc(j[1],g01))
        g12=conj(glc(j[2],g02))
        g13=conj(glc(j[3],g03))
        g1c=conj(sqrt(2*j[1]+1)*exp(-t*(lb1(j[1])))*glc(j[1],g01))
        g2c=conj(sqrt(2*j[2]+1)*exp(-t*(lb1(j[2])))*glc(j[2],g02))
        g3c=conj(sqrt(2*j[3]+1)*exp(-t*(lb1(j[3])))*glc(j[3],g03))
        a1=sqrt(2*j[1]+1)*exp(-t*(lb1(j[1])))
        a2=sqrt(2*j[2]+1)*exp(-t*(lb1(j[2])))
        a3=sqrt(2*j[3]+1)*exp(-t*(lb1(j[3])))
        as=0
        for m1 in -j[1]:j[1]
            for m2 in -j[2]:j[2]
                if abs(m1+m2)<=j[3]
                    resg=g11[Int(j[1]+m1+1),Int(j[1]-m3+1)]*g12[Int(j[2]+m2+1),Int(j[2]-m4+1)]*g13[Int(j[3]-m1-m2+1),Int(j[3]-m5+1)]
                    coef1=(-1)^(Int(j[1]+j[2]+j[3]-m1-m2-(-m1-m2)))
                    wcj=Wj1([j[1] j[2] j[3] m1 m2 -m1-m2],jml,aj1,ij1)
                    as=as+resg*coef1*wcj*a1*a2*a3
                end
            end
        end
#         println(as)
        return as
    end

function TransformM4m(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as1.=sA.*i2*i1
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
                @inbounds broadcast!(+,ia,ia,kronC4)
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM4m(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as1.=sA.*i2*i1
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
                @inbounds broadcast!(+,ia,ia,kronC4)
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM3mtt(jj,t,g01,g02,g03,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=jj[1]
        j2=jj[2]
        j3=jj[3]
        j4=jj[4]
        ia=zeros(Int(2*j4+1),a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)))*im
        as1=zeros(size(A1,1))*im
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g1c=exp(-t*(lb1(jj[1])))*adjoint(glc(jj[1],g01))*10^(-20)
        #println(g1c)
        g2c=exp(-t*(lb1(jj[2])))*adjoint(glc(jj[2],g02))*10^(-20)
        g3c=exp(-t*(lb1(jj[3])))*adjoint(glc(jj[3],g03))*10^(-20)
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j1+1)*sqrt(2*j2+1)*sqrt(2*j3+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1).*sqrt(2*j4+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64,2}(undef, Int(a),j123)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4=s1-m3
            mi5 = Int(j4+m4+1)
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
#                 print(ls1,ls1[:,3].-s1)
                Inter4_Compact(ls1, ls2, Atl, jml, a, s1, ij1, i2)
                as1.=sA.*i2*i1*(-1)^(j1+j2+j3+j4-m1-m2-m3-m4)
#                 print(as1)
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,as1,transpose(kronC2))
#                 println(size(ia),size(kronC4))
#                 print(kronC4)
                ia[mi5,:,:].=@view(ia[mi5,:,:]).+kronC3
#                 println("f2=",[j1 j2 j3 j4 m1 m2 m3 m4],ia[1,1])
#                 println("g2=",sqrt(2*j1+1)*g1c[Int32(j1+m1+1),1],sqrt(2*j2+1)*g2c[Int32(j2+m2+1),1],sqrt(2*j3+1)*g3c[Int32(j3+m3+1),1],sqrt(2*j4+1)*g4c[Int32(j4+m4+1),1])
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM4mtt(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g1c=exp(-t*(lb1(j[1])))*adjoint(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*adjoint(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*adjoint(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*adjoint(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j1+1)*sqrt(2*j2+1)*sqrt(2*j3+1)*sqrt(2*j4+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
#                 print(ls1,ls1[:,3].-s1)
                Inter4_Compact(ls1, ls2, Atl, jml, a, s1, ij1, i2)
                as1.=sA.*i2*i1*(-1)^(j1+j2+j3+j4-m1-m2-m3-m4)
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
                @inbounds broadcast!(+,ia,ia,kronC4)
#                 println("f2=",[j1 j2 j3 j4 m1 m2 m3 m4],ia[1,1])
#                 println("g2=",sqrt(2*j1+1)*g1c[Int32(j1+m1+1),1],sqrt(2*j2+1)*g2c[Int32(j2+m2+1),1],sqrt(2*j3+1)*g3c[Int32(j3+m3+1),1],sqrt(2*j4+1)*g4c[Int32(j4+m4+1),1])
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM4mtt_s(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        # g1c=exp(-t*(lb1(j[1])))*adjoint(glc(j[1],g01))*10^(-25)
        # #println(g1c)
        # g2c=exp(-t*(lb1(j[2])))*adjoint(glc(j[2],g02))*10^(-25)
        # g3c=exp(-t*(lb1(j[3])))*adjoint(glc(j[3],g03))*10^(-25)
        # g4c=exp(-t*(lb1(j[4])))*adjoint(glc(j[4],g04))*10^(-25)
        g1c=exp(-t*(lb1(j[1])))*adjoint(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*adjoint(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*adjoint(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*adjoint(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j1+1)*sqrt(2*j2+1)*sqrt(2*j3+1)*sqrt(2*j4+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
#                 print(ls1,ls1[:,3].-s1)
                Inter4_Compact(ls1, ls2, Atl, jml, a, s1, ij1, i2)
                as1.=sA.*i2*i1*(-1)^(j1+j2+j3+j4-m1-m2-m3-m4)
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
                @inbounds broadcast!(+,ia,ia,kronC4)
#                 println("f2=",[j1 j2 j3 j4 m1 m2 m3 m4],ia[1,1])
#                 println("g2=",sqrt(2*j1+1)*g1c[Int32(j1+m1+1),1],sqrt(2*j2+1)*g2c[Int32(j2+m2+1),1],sqrt(2*j3+1)*g3c[Int32(j3+m3+1),1],sqrt(2*j4+1)*g4c[Int32(j4+m4+1),1])
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

    function TransformM4mtt_sNew(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g1c=exp(-t*(lb1(j[1])))*adjoint(glc(j[1],g01))/10^25
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*adjoint(glc(j[2],g02))/10^25
        g3c=exp(-t*(lb1(j[3])))*adjoint(glc(j[3],g03))/10^25
        g4c=exp(-t*(lb1(j[4])))*adjoint(glc(j[4],g04))/10^25
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.05
        #Initialzation of ls
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j1+1)*sqrt(2*j2+1)*sqrt(2*j3+1)*sqrt(2*j4+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
#                 print(ls1,ls1[:,3].-s1)
                Inter4_Compact(ls1, ls2, Atl, jml, a, s1, ij1, i2)
                as1.=sA.*i2*i1*(-1)^(j1+j2+j3+j4-m1-m2-m3-m4)
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
                @inbounds broadcast!(+,ia,ia,kronC4)
#                 println("f2=",[j1 j2 j3 j4 m1 m2 m3 m4],ia[1,1])
#                 println("g2=",sqrt(2*j1+1)*g1c[Int32(j1+m1+1),1],sqrt(2*j2+1)*g2c[Int32(j2+m2+1),1],sqrt(2*j3+1)*g3c[Int32(j3+m3+1),1],sqrt(2*j4+1)*g4c[Int32(j4+m4+1),1])
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM5mtt(jj,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=jj[1]
        j2=jj[2]
        j3=jj[3]
        j4=jj[4]
        j5=jj[5]
        ia=zeros(Int(2*j5+1),a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g1c=exp(-t*(lb1(jj[1])))*adjoint(glc(jj[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(jj[2])))*adjoint(glc(jj[2],g02))
        g3c=exp(-t*(lb1(jj[3])))*adjoint(glc(jj[3],g03))
        g4c=exp(-t*(lb1(jj[4])))*adjoint(glc(jj[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        ls3=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=view(A1,1:a,3)
        ls3[1:a,1].=view(A1,1:a,3)
        ls3[1:a,2].=j4
        ls3[1:a,3].=j5
        i1=sqrt(2*j1+1)*sqrt(2*j2+1)*sqrt(2*j3+1)*sqrt(2*j4+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1).*sqrt.((A1[:,3].*2).+1).*sqrt(2*j5+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            m4=j4-rem(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1))
            s1 = -m1-m2
            s2=s1-m3
            m5 = s2-m4
            mi5 = Int(j5+m5+1)
            if abs(m5)<=j5
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=s2
                @view(ls3[1:a,4]).=-s2
                @view(ls3[1:a,5]).=m4  
                @view(ls3[1:a,6]).=m5
            #Block Three
                #println(g2c)
#                 print(ls1,ls1[:,3].-s1)
                Inter5_Compact(ls1, ls2, ls3, Atl, jml, a, s1, s2, ij1, i2)
                as1.=sA.*i2*i1*(-1)^(j1+j2+j3+j4+j5-m1-m2-m3-m4-m5)
#                 print(as1)
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,as1,transpose(kronC3))
#                 println(size(ia),size(kronC4))
#                 print(kronC4)
                ia[mi5,:,:].=@view(ia[mi5,:,:]).+kronC4
#                 println("f2=",[j1 j2 j3 j4 m1 m2 m3 m4],ia[1,1])
#                 println("g2=",sqrt(2*j1+1)*g1c[Int32(j1+m1+1),1],sqrt(2*j2+1)*g2c[Int32(j2+m2+1),1],sqrt(2*j3+1)*g3c[Int32(j3+m3+1),1],sqrt(2*j4+1)*g4c[Int32(j4+m4+1),1])
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

    function TransformM6mNew(j,t,g01,g02,g03,g04,g05,g06,jml,mv123,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        j5=j[5]
        j6=j[6]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1)*(2*j6+1)))*im
        as1=zeros(size(A1,1))*im
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g5c=exp(-t*(lb1(j[5])))*conj(glc(j[5],g05))
        g6c=exp(-t*(lb1(j[6])))*conj(glc(j[6],g06))
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=j[1]
        m2=0.0
        m3=0.0
        m4=0.0
        m5=0.0
        m5p=0.0
        m6=0.0
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        ls3=zeros(a,6)
        ls4=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=view(A1,1:a,3)
        ls3[1:a,1].=view(A1,1:a,4)
        ls3[1:a,2].=j4
        ls3[1:a,3].=view(A1,1:a,3)
        ls4[1:a,1].=j5
        ls4[1:a,2].=j6
        ls4[1:a,3].=view(A1,1:a,4)
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1).*sqrt.((A1[:,3].*2).+1).*sqrt.((A1[:,4].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        j12345=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1))
        j123456=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1)*(2*j6+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64}(undef, j12345)
        kronC5=Array{ComplexF64}(undef, j123456)
        kronC6=Array{ComplexF64,2}(undef, Int(a),j123456)
        for ma in 1:(2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1)
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            m4=j4-rem(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1))
            m5=j5-rem(fld(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1)),(2*j5+1))
            s1 = -m1-m2
            s2 = s1-m3
            s3 = s2-m4
            m6 = s3-m5
            if abs(m6)<=j6
            #Block Two    
                ls1[1:a,4].=m1
                ls1[1:a,5].=m2
                ls1[1:a,6].=s1
                ls2[1:a,4].=-s1
                ls2[1:a,5].=m3  
                ls2[1:a,6].=s2
                ls3[1:a,4].=s3
                ls3[1:a,5].=m4  
                ls3[1:a,6].=-s2
                ls4[1:a,4].=m5
                ls4[1:a,5].=m6
                ls4[1:a,6].=-s3
#                 if j6+m6+1==1.5
#                     print([j1 j2 j3 j4 j5 j6; m1 m2 m3 m4 m5 m6])
#                 end
            #Block Three
                Inter6_mod2(ls1,ls2,ls3,ls4,Atl,jml,a,s1,s2,s3,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi, ci, di) 
                as1.=sA.*i2*i1
                @inbounds kron!(kronC1,@view(g1c[Int32(j1+m1+1),:]),@view(g2c[Int32(j2+m2+1),:]))
                @inbounds kron!(kronC2,kronC1,@view(g3c[Int32(j3+m3+1),:]))
                @inbounds kron!(kronC3,kronC2,@view(g4c[Int32(j4+m4+1),:]))
                @inbounds kron!(kronC4,kronC3,@view(g5c[Int32(j5+m5+1),:]))
                @inbounds kron!(kronC5,kronC4,@view(g6c[Int32(j6+m6+1),:]))
                @inbounds kron!(kronC6,as1,transpose(kronC5))
                @inbounds broadcast!(+,ia,ia,kronC6)
            end
        end
        #println(count1)
        return (Alist,ia)
    end

    function TransformM6Final(j,g1c,g2c,g3c,mv123)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        j5=j[5]
        j6=j[6]
        #count2=0
        #count3=0
            #Block One
        m1=mv123[1]
        m2=mv123[2]
        m3=mv123[3]
        m4=mv123[4]
        m5=mv123[5]
        m6=mv123[6]
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)
            #Block Three
        #println(count1)
        return i1*(-1)^(j1-m1)*(-1)^(j2-m2)*(-1)^(j3-m3)*g1c[Int32(j1-m1+1),Int32(j4+m4+1)]*g2c[Int32(j2-m2+1),Int32(j5+m5+1)]*g3c[Int32(j3-m3+1),Int32(j6+m6+1)]
        # return g1c[Int32(j1+m1+1),Int32(j4+m4+1)]*g2c[Int32(j2+m2+1),Int32(j5+m5+1)]*g3c[Int32(j3+m3+1),Int32(j6+m6+1)]
    end

    function TransformM6Final2(j,g1c,g2c,g3c,mv123)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        j5=j[5]
        j6=j[6]
        #count2=0
        #count3=0
            #Block One
        m1=mv123[1]
        m2=mv123[2]
        m3=mv123[3]
        m4=mv123[4]
        m5=mv123[5]
        m6=mv123[6]
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)
            #Block Three
        #println(count1)
        return i1*g1c[Int32(j1+m1+1),Int32(j4+m4+1)]*g2c[Int32(j2+m2+1),Int32(j5+m5+1)]*g3c[Int32(j3+m3+1),Int32(j6+m6+1)]
    end

function TransformM4m_SVD(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros(size(A1,1))*im
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        gk1=kron(g1c,g2c)
        gk2=kron(g3c,g4c)
        (g1u,g1s,g1v)=svdm1(gk1)
        (g2u,g2s,g2v)=svdm1(gk2)
#         print(size(s1r))
#         print(sum((gk1.-(u1r*Diagonal(s1r)*v1r.+u1i*Diagonal(s1i)*v1i.*im))./gk1)/size(gk1,1)/size(gk1,1))
        sizes1=size(g1s,1)
        sizes2=size(g2s,1)
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        km1=zeros(sizes1*sizes2)
        km2=zeros(a,sizes1*sizes2)
        skk1=zeros(Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))
        skk2=zeros(Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))
        ar=zeros(a,sizes1*sizes2)
        resk=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))
        sv1=Diagonal(g1s)*adjoint(g1v)
        sv2=Diagonal(g2s)*adjoint(g2v)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as1.=sA.*i2*i1
                is1=Int32(j1+m1)*Int32(2*j2+1)+Int32(j2+m2+1)
                is2=Int32(j3+m3)*Int32(2*j4+1)+Int32(j4+m4+1)
                kron!(km1,exp.(u1r[is1,:]),exp.(u2r[is2,:]))
                kron!(km2,as1,transpose(km1))
                ar.=ar.+log.(km2)
            end
        end

        for i in 1:sizes1*sizes2
            si1=rem(i-1,sizes1)+1
            si2=fld(i-1,sizes1)+1
#             print((j1,j2,j3,j4),(si1,si2))
#             print(size(gk1),size(gk2))
#             print(size(skkr1),size(sv1r),size(sv2r))
            kron!(skk1,exp.(sv1[si1,:]),exp.(sv2[si2,:]))
            kron!(resk,ar[:,i],transpose(skk1))
            ia.=ia.+resk
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end


function TransformM4m_single(j,ms,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2,A1,a,t3)
        Alist=Array{ComplexF64,1}(undef,a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a)*im
        as1=zeros(size(A1,1))*im
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=Array{Float64,2}(undef,a,6)
        ls2=Array{Float64,2}(undef,a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as1.=sA.*i2*i1*g1c[Int32(j1+m1+1),Int32(j1+ms[1]+1)]*g2c[Int32(j2+m2+1),Int32(j2+ms[2]+1)]*g3c[Int32(j3+m3+1),Int32(j3+ms[3]+1)]*g4c[Int32(j4+m4+1),Int32(j4+ms[4]+1)]
                ia.=ia.+as1
                t3[1]=t3[1]+1
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

function TransformM4m_matmul(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai4(j)
        a=size(A1,1)
        Alist=zeros(a)+im*zeros(a)
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        as1=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        sA=sqrt.((A1[:,2].*2).+1)
        j12=Int((2*j1+1)*(2*j2+1))
        j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
        j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
        kronC1=Array{ComplexF64,2}(undef, Int(2*j1+1),Int(2*j2+1))
        kronC2=Array{ComplexF64,2}(undef, j12, Int(2*j3+1))
        kronC3=Array{ComplexF64,2}(undef, j123, Int(2*j4+1))
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                as1.=sA.*i2*i1
                @inbounds matmul!(kronC1,Matrix(reshape(@view(g1c[Int32(j1+m1+1),:]),(Int(2*j1+1),1))),Matrix(reshape(@view(g2c[Int32(j2+m2+1),:]),(1,Int(2*j2+1)))))
                @inbounds matmul!(kronC2,Matrix(reshape(transpose(kronC1),(j12,1))),Matrix(reshape(g3c[Int32(j3+m3+1),:],(1,Int(2*j3+1)))))
                @inbounds matmul!(kronC3,Matrix(reshape(transpose(kronC2),(j123,1))),Matrix(reshape(g4c[Int32(j4+m4+1),:],(1,Int(2*j4+1)))))
                @inbounds matmul!(kronC4,Matrix(reshape(as1,(a,1))),Matrix(reshape(transpose(kronC3),(1,j1234))))
                @inbounds ia.=ia.+kronC4
            end
        end
        as1=nothing
        ls1=nothing
        ls2=nothing
        sA=nothing
        ij1=nothing
        i2=nothing
        bi=nothing
        ci=nothing
        di=nothing
        ja1=nothing
        ja2=nothing
        ja3=nothing
        ja=nothing
        i=nothing
        j=nothing
        #println(count1)
        return (Alist,ia)
    end

# function TransformM4m_tulio(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2)
#         A1=Fai4(j)
#         a=size(A1,1)
#         Alist=zeros(a)+im*zeros(a)
#         j1=j[1]
#         j2=j[2]
#         j3=j[3]
#         j4=j[4]
#         ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
#         as1=zeros((size(A1,1)))+im*zeros((size(A1,1)))
#         g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
#         #println(g1c)
#         g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
#         g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
#         g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
#         #println(g2c)
#         #println(g3c)
#         #println(g4c)
#         count1=0
#         #count2=0
#         #count3=0
#         i1=0.0
#         i3=0.0
#         m1=0.0
#         m2=0.0
#         m3=0.0
#         m3p=0.0
#         ls1=zeros(a,6)
#         ls2=zeros(a,6)
#         #Initialzation of ls
#         ls1[1:a,1].=j1
#         ls1[1:a,2].=j2
#         ls1[1:a,3].=view(A1,1:a,2)
#         ls2[1:a,1].=view(A1,1:a,2)
#         ls2[1:a,2].=j3
#         ls2[1:a,3].=j4
#         i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
#         ij1=Array{Int32}(undef, a, 2)
#         i2=Array{Float64,1}(undef,a)
#         bi=Array{Float64,1}(undef,a)
#         ci=Array{Float64,1}(undef,a)
#         di=Array{Float64,1}(undef,a)
#         ja1=Array{Float64,1}(undef,a)
#         ja2=Array{Float64,1}(undef,a)
#         ja3=Array{Float64,1}(undef,a)
#         ja=Array{Float64,1}(undef,a)
#         i=Array{Float64,1}(undef,a)
#         j=Array{Float64,1}(undef,a)
#         sA=sqrt.((A1[:,2].*2).+1)
#         j12=Int((2*j1+1)*(2*j2+1))
#         j123=Int((2*j1+1)*(2*j2+1)*(2*j3+1))
#         j1234=Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
#         Q1=Array{ComplexF64,1}(undef,Int(2*j1+1))
#         Q2=Array{ComplexF64,1}(undef,Int(2*j2+1))
#         Q3=Array{ComplexF64,1}(undef,Int(2*j3+1))
#         Q4=Array{ComplexF64,1}(undef,Int(2*j4+1))
#         kronC1=Array{ComplexF64,2}(undef, Int(2*j1+1),Int(2*j2+1))
#         kronC2=Array{ComplexF64,2}(undef, j12, Int(2*j3+1))
#         kronC3=Array{ComplexF64,2}(undef, j123, Int(2*j4+1))
#         kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
#         for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
#             #Block One
#             m1=j1-rem(ma-1,(2*j1+1))
#             m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
#             m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
#             s1 = -m1-m2
#             m4 = s1-m3
#             if abs(m4)<=j4
#             #Block Two    
#                 @view(ls1[1:a,4]).=m1
#                 @view(ls1[1:a,5]).=m2
#                 @view(ls1[1:a,6]).=s1
#                 @view(ls2[1:a,4]).=-s1
#                 @view(ls2[1:a,5]).=m3  
#                 @view(ls2[1:a,6]).=m4
#             #Block Three
#                 #println(g2c)
#                 Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
#                 as1.=sA.*i2*i1
#                 Q1.=@view(g1c[Int32(j1+m1+1),:])
#                 Q2.=@view(g2c[Int32(j2+m2+1),:])
#                 Q3.=@view(g3c[Int32(j3+m3+1),:])
#                 Q4.=@view(g4c[Int32(j4+m4+1),:])
#                 @tullio kronC1[x,y] := Q1[x] * Q2[y]
#                 @tullio kronC2[x,y] := transpose(kronC1)[x] * Q3[y]
#                 @tullio kronC3[x,y] := transpose(kronC2)[x] * Q4[y]
#                 @tullio kronC4[x,y] := as1[x] * transpose(kronC3)[y]
#                 @inbounds ia.=ia.+kronC4
#             end
#         end
#         as1=nothing
#         ls1=nothing
#         ls2=nothing
#         sA=nothing
#         ij1=nothing
#         i2=nothing
#         bi=nothing
#         ci=nothing
#         di=nothing
#         ja1=nothing
#         ja2=nothing
#         ja3=nothing
#         ja=nothing
#         i=nothing
#         j=nothing
#         #println(count1)
#         return (Alist,ia)
#     end

function TransformM4m_CUDA(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai4(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        as1=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        a=size(A1,1)
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        ia=zeros(a,Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)))*im
        sA=sqrt.((A1[:,2].*2).+1)
        for ma in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1))
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                @view(ls1[1:a,4]).=m1
                @view(ls1[1:a,5]).=m2
                @view(ls1[1:a,6]).=s1
                @view(ls2[1:a,4]).=-s1
                @view(ls2[1:a,5]).=m3  
                @view(ls2[1:a,6]).=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                for mb in 1:Int((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))
                    m21=j1-rem(mb-1,(2*j1+1))
                    m22=j2-rem(fld(mb-1,(2*j1+1)),(2*j2+1))
                    m23=j3-rem(fld(fld(mb-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
                    m24=j4-rem(fld(fld(fld(mb-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1))
                    @inbounds i3=g1c[Int32(j1+m1+1),Int32(j1+m21+1)]*g2c[Int32(j2+m2+1),Int32(j2+m22+1)]*g3c[Int32(j3+m3+1),Int32(j3+m23+1)]*g4c[Int32(j4+m4+1),Int32(j4+m24+1)]
                    @inbounds Alist=@view ia[:,Int(mb)]
                    @inbounds as1.=(i2*i3*i1).*sA
                    #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                    #print([float(i1) float(i2) float(i3)],"\n")
                    #println()
                    #println([m1 m2 m3])
                    #println(i1)
                    #println(i2)
                    #println(i3)
                    @inbounds Alist.=(Alist.+as1)
                    #count1=count1+1
                    if m3 != m3p
                        #print(m5p)
                    end
                    m3p=m3
                end
            end
        end
        #println(count1)
        return (Alist,ia)
    end

function TransformM4m2(j,t,g01,g02,g03,g04,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai4(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        as=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        #println(g1c)
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        #println(g2c)
        #println(g3c)
        #println(g4c)
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m3p=0.0
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        a=size(A1,1)
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=j4
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        for ma in 1:(2*j1+1)*(2*j2+1)*(2*j3+1)
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            s1 = -m1-m2
            m4 = s1-m3
            if abs(m4)<=j4
            #Block Two    
                ls1[1:a,4].=m1
                ls1[1:a,5].=m2
                ls1[1:a,6].=s1
                ls2[1:a,4].=-s1
                ls2[1:a,5].=m3  
                ls2[1:a,6].=m4
            #Block Three
                #println(g2c)
                Inter4(ls1,ls2,Atl,jml,a,s1,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi)
                    i3=g1c[Int32(j1+m21+1),Int32(j1+m1+1)]*g2c[:,Int32(j2+m2+1)]*sum(g3c[:,Int32(j3+m3+1)])*sum(g4c[:,Int32(j4+m4+1)])
                    as.=(i2.*i3).*i1
                    #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                    #print([float(i1) float(i2) float(i3)],"\n")
                    #println()
                    #println([m1 m2 m3])
                    #println(i1)
                    #println(i2)
                    #println(i3)
                    Alist=(Alist.+as)
                    #count1=count1+1
                    if m3 != m3p
                        #print(m5p)
                    end
                    m3p=m3
            end
        end
        #println(count1)
        return Alist
    end


    # 6-dimensional transformation matrix
    function TransformM6(j,t,g01,g02,g03,g04,g05,g06)
        A1=Fai6(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1=glc(j[1],g01)
        g2=glc(j[2],g02)
        g3=glc(j[3],g03)
        g4=glc(j[4],g04)
        g5=glc(j[5],g05)
        g6=glc(j[6],g06)
        for a in 1:size(A1,1)
            sum=BigFloat(0.0)+im*BigFloat(0.0)
            for m1 in -j[1]:j[1]
                for m2 in -j[2]:j[2]
                    for m3 in -j[3]:j[3]
                        for m4 in -j[4]:j[4]
                            for m5 in -j[5]:j[5]
                                m6 = -m1-m2-m3-m4-m5
                                if -j[6]<=m6<=j[6]
                                    i2=Inter6_mod(j,A1[a,:],[m1 m2 m3 m4 m5])
                                    if i2 != 0
                                        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)*sqrt(2*j[4]+1)*sqrt(2*j[5]+1)*sqrt(2*j[6]+1)*exp(-t*(lb1(j[1])+lb1(j[2])+lb1(j[3])+lb1(j[4])+lb1(j[5])+lb1(j[6])))
                                        i3=conj(g1[1,Int32(j[1]-m1+1)])*conj(g2[1,Int32(j[2]-m2+1)])*conj(g3[1,Int32(j[3]-m3+1)])*conj(g4[1,Int32(j[4]-m4+1)])*conj(g5[1,Int32(j[5]-m5+1)])*conj(g6[1,Int32(j[6]-m6+1)])
                                        as=i1*i2*i3
                                        #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                                        #print([float(i1) float(i2) float(i3)],"\n")
                                        #println()
                                        sum=sum+as
                                    end
                                end
                            end
                        end
                    end
                end
            end
            Alist[a]=sum
        end
        return Alist
    end

    # 6-dimensional transformation matrix (3 flower graph)
    function TransformM6f_mod3(j,t,g01,g02,g03,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai6(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        as=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m4=0.0
        m5=0.0
        m5p=0.0
        m6=0.0
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        j5=j[5]
        j6=j[6]
        a=size(A1,1)
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        ls3=zeros(a,6)
        ls4=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=view(A1,1:a,3)
        ls3[1:a,1].=view(A1,1:a,4)
        ls3[1:a,2].=j4
        ls3[1:a,3].=view(A1,1:a,3)
        ls4[1:a,1].=j5
        ls4[1:a,2].=j6
        ls4[1:a,3].=view(A1,1:a,4)
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        for ma in 1:(2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1)
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            m4=j4-rem(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1))
            m5=j5-rem(fld(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1)),(2*j5+1))
            s1 = -m1-m2
            s2 = s1-m3
            s3 = s2-m4
            m6 = s3-m5
            if abs(m6)<=j6
            #Block Two    
                ls1[1:a,4].=m1
                ls1[1:a,5].=m2
                ls1[1:a,6].=s1
                ls2[1:a,4].=-s1
                ls2[1:a,5].=m3  
                ls2[1:a,6].=s2
                ls3[1:a,4].=s3
                ls3[1:a,5].=m4  
                ls3[1:a,6].=-s2
                ls4[1:a,4].=m5
                ls4[1:a,5].=m6
                ls4[1:a,6].=-s3
            #Block Three
                i3=(-1)^(j1-m1)*(-1)^(j2-m2)*(-1)^(j3-m3)*g1c[Int32(j1-m1+1),Int32(j4+m4+1)]*g2c[Int32(j2-m2+1),Int32(j5+m5+1)]*g3c[Int32(j3-m3+1),Int32(j6+m6+1)]
                Inter6_mod2(ls1,ls2,ls3,ls4,Atl,jml,a,s1,s2,s3,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi, ci, di) 
                as.=(i2.*i3).*i1
                        #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                        #print([float(i1) float(i2) float(i3)],"\n")
                        #println()
                Alist=(Alist.+as)
                count1=count1+1
                if m5 != m5p
                    #print(m5p)
                end
                m5p=m5
            end
        end
        #println(count1)
        return Alist
    end

    function TransformM6m(j,t,g01,g02,g03,g04,g05,g06,jml,Atl,rs1,rs11,rs12,rs2)
        A1=Fai6(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        as=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1c=exp(-t*(lb1(j[1])))*conj(glc(j[1],g01))
        g2c=exp(-t*(lb1(j[2])))*conj(glc(j[2],g02))
        g3c=exp(-t*(lb1(j[3])))*conj(glc(j[3],g03))
        g4c=exp(-t*(lb1(j[4])))*conj(glc(j[4],g04))
        g5c=exp(-t*(lb1(j[5])))*conj(glc(j[5],g05))
        g6c=exp(-t*(lb1(j[6])))*conj(glc(j[6],g06))
        count1=0
        #count2=0
        #count3=0
        i1=0.0
        i3=0.0
        m1=0.0
        m2=0.0
        m3=0.0
        m4=0.0
        m5=0.0
        m5p=0.0
        m6=0.0
        j1=j[1]
        j2=j[2]
        j3=j[3]
        j4=j[4]
        j5=j[5]
        j6=j[6]
        a=size(A1,1)
        ls1=zeros(a,6)
        ls2=zeros(a,6)
        ls3=zeros(a,6)
        ls4=zeros(a,6)
        #Initialzation of ls
        ls1[1:a,1].=j1
        ls1[1:a,2].=j2
        ls1[1:a,3].=view(A1,1:a,2)
        ls2[1:a,1].=view(A1,1:a,2)
        ls2[1:a,2].=j3
        ls2[1:a,3].=view(A1,1:a,3)
        ls3[1:a,1].=view(A1,1:a,4)
        ls3[1:a,2].=j4
        ls3[1:a,3].=view(A1,1:a,3)
        ls4[1:a,1].=j5
        ls4[1:a,2].=j6
        ls4[1:a,3].=view(A1,1:a,4)
        i1=sqrt(2*j[1]+1)*sqrt(2*j[2]+1)*sqrt(2*j[3]+1)
        ij1=Array{Int32}(undef, a, 2)
        i2=Array{Float64,1}(undef,a)
        bi=Array{Float64,1}(undef,a)
        ci=Array{Float64,1}(undef,a)
        di=Array{Float64,1}(undef,a)
        ja1=Array{Float64,1}(undef,a)
        ja2=Array{Float64,1}(undef,a)
        ja3=Array{Float64,1}(undef,a)
        ja=Array{Float64,1}(undef,a)
        i=Array{Float64,1}(undef,a)
        j=Array{Float64,1}(undef,a)
        kronC1=Array{ComplexF64}(undef, j12)
        kronC2=Array{ComplexF64}(undef, j123)
        kronC3=Array{ComplexF64}(undef, j1234)
        kronC4=Array{ComplexF64,2}(undef, Int(a),j1234)
        for ma in 1:(2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1)*(2*j5+1)
            #Block One
            m1=j1-rem(ma-1,(2*j1+1))
            m2=j2-rem(fld(ma-1,(2*j1+1)),(2*j2+1))
            m3=j3-rem(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1))
            m4=j4-rem(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1))
            m5=j5-rem(fld(fld(fld(fld(ma-1,(2*j1+1)),(2*j2+1)),(2*j3+1)),(2*j4+1)),(2*j5+1))
            s1 = -m1-m2
            s2 = s1-m3
            s3 = s2-m4
            m6 = s3-m5
            if abs(m6)<=j6
            #Block Two    
                ls1[1:a,4].=m1
                ls1[1:a,5].=m2
                ls1[1:a,6].=s1
                ls2[1:a,4].=-s1
                ls2[1:a,5].=m3  
                ls2[1:a,6].=s2
                ls3[1:a,4].=s3
                ls3[1:a,5].=m4  
                ls3[1:a,6].=-s2
                ls4[1:a,4].=m5
                ls4[1:a,5].=m6
                ls4[1:a,6].=-s3
            #Block Three
                i3=(-1)^(j1-m1)*(-1)^(j2-m2)*(-1)^(j3-m3)*g1c[Int32(j1-m1+1),Int32(j4+m4+1)]*g2c[Int32(j2-m2+1),Int32(j5+m5+1)]*g3c[Int32(j3-m3+1),Int32(j6+m6+1)]
                Inter6_mod2(ls1,ls2,ls3,ls4,Atl,jml,a,s1,s2,s3,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi, ci, di) 
                as.=(i2.*i3).*i1
                        #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                        #print([float(i1) float(i2) float(i3)],"\n")
                        #println()
                Alist=(Alist.+as)
                count1=count1+1
                if m5 != m5p
                    #print(m5p)
                end
                m5p=m5
            end
        end
        #println(count1)
        return Alist
    end

    # 4-dimensional transformation matrix (alternate version)
    function TransformM4p(j,t,g01,g02,g03,g04)
        A1=Fai4(j)
        Alist=zeros((size(A1,1)))+im*zeros((size(A1,1)))
        g1=glc(j[1],g01)
        g2=glc(j[1],g02)
        g3=glc(j[1],g03)
        g4=glc(j[1],g04)
        i1=sqrt(2j[1]+1)*sqrt(2j[2]+1)*sqrt(2j[3]+1)*sqrt(2j[4]+1)*exp(-t*(lb1(j[1])+lb1(j[2])+lb1(j[3])+lb1(j[4])))
        for a in 1:size(A1,1)
            sum=BigFloat(0.0)+im*BigFloat(0.0)
            for m1 in -j[1]:j[1]
                for m2 in -j[2]:j[2]
                    for m3 in -j[3]:j[3]
                        m4=-m1-m2-m3
                        if -j[4]<=m4<=j[4]

                            i2=Inter4(j,A1[a,:],[m1 m2 m3 m4])
                            i3=conj(g1[1,Int32(j[1]-m1+1)])*conj(g2[1,Int32(j[2]-m2+1)])*conj(g3[1,Int32(j[3]-m3+1)])*conj(g4[1,Int32(j[4]-m4+1)])
                            as=i2*i3
                                #print([a j[1]-m1+1 j[2]-m2+1 j[3]-m3+1 j[4]-m4+1],float(as),"\n")
                                #print([float(i1) float(i2) float(i3)],"\n")
                                #println()
                                sum=sum+as
                        end
                    end
                end
            end
            Alist[a]=sum
        end
        return Alist*i1
    end

    #Calculate the normalization factor for 3-flower graph
    function Amp1(j,t,g1,g2,g3,jm,aj1)
        A1=Fai6(j)
        a=size(A1,1)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        rs123=TransformM6f_mod3(j,t1,g1,g2,g3,jm,aj1,rs1,rs11,rs12,rs2)
        return transpose(conj(rs123))*(((A1[:,2].*2).+1).*((A1[:,3].*2).+1).*((A1[:,4].*2).+1).*rs123)
    end
    
    #calculate the volume operator for 3-flower graph
    function TransM1(j,t,g1,g2,g3,jm,aj1)
        A1=Fai6(j)
        a=size(A1,1)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        rs123=TransformM6f_mod3(j,t1,g1,g2,g3,jm,aj1,rs1,rs11,rs12,rs2)
        vol1q=v123(j,A1,1,2,3,6)+v123(j,A1,1,3,5,6)+v123(j,A1,1,5,6,6)+v123(j,A1,1,2,6,6)-v123(j,A1,3,4,5,6)+v123(j,A1,4,5,6,6)+v123(j,A1,2,4,6,6)+v123(j,A1,2,3,4,6)
        return transpose(conj(rs123))*vol1q*(((A1[:,2].*2).+1).*((A1[:,3].*2).+1).*((A1[:,4].*2).+1).*rs123)
    end
    
    #Calculate the normalization factor for 2-flower graph
    function Amp4(j,t,g1,g2,jm,aj1)
        A1=Fai4(j)
        a=size(A1,1)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        rs123=TransformM4(j,t1,g1,g2,jm,aj1,rs1,rs11,rs12,rs2)
        return transpose(conj(rs123))*(((A1[:,2].*2).+1).*rs123)
    end

    #calculate the volume operator for 2-flower graph
    function TransM4(j,t,g1,g2,g3,jm,aj1)
        A1=Fai6(j)
        a=size(A1,1)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        rs123=TransformM6f_mod3(j,t1,g1,g2,g3,jm,aj1,rs1,rs11,rs12,rs2)
        vol1q=v123(j,A1,1,2,3,6)+v123(j,A1,1,3,5,6)+v123(j,A1,1,5,6,6)+v123(j,A1,1,2,6,6)-v123(j,A1,3,4,5,6)+v123(j,A1,4,5,6,6)+v123(j,A1,2,4,6,6)+v123(j,A1,2,3,4,6)
        return transpose(conj(rs123))*vol1q*(((A1[:,2].*2).+1).*((A1[:,3].*2).+1).*((A1[:,4].*2).+1).*rs123)
    end
    
    #Calculate the normalization factor for 4 vertex
    function Amp4o(j,t1,g1,g2,g3,g4,jm,aj1)
        A1=Fai4(j)
        if size(A1)[1] != 0
            a=size(A1,1)
            rs1=Array{Int32}(undef, a, 2)
            rs11=Array{Int32}(undef, a)
            rs12=Array{Int32}(undef, a)
            rs2=zeros(a)
            rs123=TransformM4m2(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2)
            return transpose(conj(rs123))*(rs123)
        else
            return 0.0
        end
    end

    function Amp4o2(j,t1,g1,g2,g3,g4,jm,aj1)
        A1=Fai4(j)
        if size(A1)[1] != 0
            a=size(A1,1)
            rs1=Array{Int32}(undef, a, 2)
            rs11=Array{Int32}(undef, a)
            rs12=Array{Int32}(undef, a)
            rs2=zeros(a)
            (rs123,abc)=TransformM4m(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2)
            #print(abc)
            sum(sum(abs.(abc).^2))
        else
            return 0.0
        end
    end

    #Calculate the transition factors for 4 vertex
    function Trans4o2(j,t,g1,g2,g3,g4,jm,aj1)
        A1=Fai4(j)
        a=size(A1,1)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4m_tulio(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2)
        return abc
    end

    function Trans4o(j,t,g1,g2,g3,g4,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4m(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end

    function Trans3ott(j,t1,g1,g2,g3,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM3mtt(j,t1,g1,g2,g3,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end

    function Trans4ott(j,t1,g1,g2,g3,g4,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4mtt(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end

    function Trans4ott_s(j,t1,g1,g2,g3,g4,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4mtt_s(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end

    function Trans4ott_sNew(j,t1,g1,g2,g3,g4,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4mtt_sNew(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end

    function Trans5ott(j,t1,g1,g2,g3,g4,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM5mtt(j,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a)
        return abc
    end


    function Trans6o(j,t,g1,g2,g3,g4,g5,g6,jm,aj1,A1,a)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        rs123=TransformM6m(j,t,g1,g2,g3,g4,g5,g6,jm,aj1,rs1,rs11,rs12,rs2)
        return abc
    end
    
    function Trans4o_single(j,ms,t,g1,g2,g3,g4,jm,aj1,A1,a,t3)
        rs1=Array{Int32}(undef, a, 2)
        rs11=Array{Int32}(undef, a)
        rs12=Array{Int32}(undef, a)
        rs2=zeros(a)
        (rs123,abc)=TransformM4m_single(j,ms,t1,g1,g2,g3,g4,jm,aj1,rs1,rs11,rs12,rs2,A1,a,t3)
        return abc
    end
    
    #Calculate expectation value

function Ampt4o(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1))
    ntest=[]
    t=0
    for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        t=t+1
        if rem(t,1000)==1
            print(t," ")
        end
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
            at1[i1]=Amp4o(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1)
        else
            at1[i1]=0.0+im*0.0
        end
    end
    at1
end