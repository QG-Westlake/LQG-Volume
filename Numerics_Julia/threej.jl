using LinearAlgebra

    function tjst(j1, j2, j3, m1, m2, m3)
        if abs(j1-j2) <= j3 <= abs(j1+j2) && abs(m1)<=j1&&abs(m2)<=j2&&abs(m3)<=j3&&m3==-m1-m2
            return wigner3j(j1, j2, j3, m1, m2, m3)
        else
            return 0
        end
    end

    function tjsl(jl)
        j=Int32(jl)
        a1 = zeros((2j+1, 2j+1, 4j+1, 2j+1, 2j+1, 4j+1))

        for j1 in 1:2j+1
            for j2 in 1:2j+1
                for j3 in 1:4j+1
                    for m1 in 1:2j+1
                        for m2 in 1:2j+1 
                            for m3 in 1:4j+1
                                a1[j1,j2,j3,m1,m2,m3]=tjst((j1-1)/2,(j2-1)/2,(j3-1)/2,(j1-1)/2-(m1-1),(j2-1)/2-(m2-1),(j3-1)/2-(m3-1))
                            end
                        end
                    end
                end
            end
        end
        return float(a1)
    end

    function jtoi(jmax,j,m,wlist)
        i=(2*j[1]+1)+(2*jmax[1]+1)*(2*j[2])+(2*jmax[1]+1)*(2*jmax[2]+1)*(2*j[3])
        j=(j[1]-m[1]+1)+(2*j[1]+1)*(j[2]-m[2])
        wlist[1]=Int32(i)
        wlist[2]=Int32(j)
        return 0
    end

    function jtoi_new(jmax,j,a,wlist)
        i=((view(j,1:a,1).*2).+1).+(view(j,1:a,2).*2).*(2*jmax[1]+1).+(view(j,a,3).*2).*((2*jmax[1]+1)*(2*jmax[2]+1))
        j=(view(j,1:a,1).-view(j,1:a,4).+1).+((view(j,1:a,1).*2).+1).*(view(j,1:a,2).-view(j,1:a,5))
        wlist[1:a,1]=broadcast(Int32,i)
        wlist[1:a,2]=broadcast(Int32,j)
        return 0
    end

    function itoj(ai,jmax,jlist)
        j1 = (rem(ai[1]-1,(2*jmax[1]+1)))/2
        j2 = (rem(fld(ai[1]-1,(2*jmax[1]+1)),(2*jmax[2]+1)))/2
        j3 = (fld(fld(ai[1]-1,(2*jmax[1]+1)),(2*jmax[2]+1)))/2
        m1 = j1-rem(ai[2]-1,Int32(2*j1+1)) 
        m2 = j2-fld(ai[2]-1,Int32(2*j1+1))
        jlist[1]=j1
        jlist[2]=j2 
        jlist[3]=j3
        jlist[4]=m1
        jlist[5]=m2
        jlist[6]=-m1-m2
        return 0
    end

    function itoj6(ai,jmax,jlist)
        j1 = (rem(ai-1,(2*jmax[1]+1)))/2
        j2 = (rem(fld(ai-1,(2*jmax[1]+1)),(2*jmax[2]+1)))/2
        j3 = (rem(fld(fld(ai-1,(2*jmax[1]+1)),(2*jmax[2]+1)),(2*jmax[3]+1)))/2
        j5 = (rem(fld(fld(fld(ai-1,(2*jmax[1]+1)),(2*jmax[2]+1)),(2*jmax[3]+1)),(2*jmax[4]+1)))/2
        j6 = (fld(fld(fld(fld(ai-1,(2*jmax[1]+1)),(2*jmax[2]+1)),(2*jmax[3]+1)),(2*jmax[4]+1)))/2
        jlist[1]=j1
        jlist[2]=j2 
        jlist[3]=j3
        jlist[4]=1
        jlist[5]=j5
        jlist[6]=j6
        return 0
    end

    function jtoi6(jmax,j)
        i=(2*j[1]+1)+(2*jmax[1]+1)*(2*j[2])+(2*jmax[1]+1)*(2*jmax[2]+1)*(2*j[3])+(2*jmax[1]+1)*(2*jmax[2]+1)*(2*jmax[3]+1)*(2*j[4])+(2*jmax[1]+1)*(2*jmax[2]+1)*(2*jmax[3]+1)*(2*jmax[4]+1)*(2*j[5])
        return Int(i)
    end
    
    function tjsl2(jl,a1)
        jlist = [0.0 0.0 0.0 0.0 0.0 0.0]
        j1=0.0
        j2=0.0
        j3=0.0
        m1=0.0
        m2=0.0
        for j1 in 1:size(a1)[1]
            for m1 in 1:size(a1)[2]
                itoj([j1 m1],jl,jlist)
                if abs(jlist[4])<=jlist[1]&&abs(jlist[5])<=jlist[2]&&abs(jlist[6])<=jlist[3]&&rem(2*(jlist[1]+jlist[2]+jlist[3]),2)==0
                    a1[j1,m1]=tjst(jlist[1],jlist[2],jlist[3],jlist[4],jlist[5],jlist[6])
                end
            end
        end
        return 0
    end

    function sixj1(jl,a1)
        jlist = [0.0 0.0 0.0 0.0 0.0 0.0]
        for i1 in 1:size(a1)[1]
            itoj6(i1,jl,jlist)
            j1=jlist[1]
            j2=jlist[2]
            j3=jlist[3]
            j5=jlist[5]
            j6=jlist[6]
            if abs(j2-j1)<=j3<=abs(j2+j1)&&abs(j3-1)<=j5<=abs(j3+1)&&abs(j2-1)<=j6<=abs(j2+1)&&abs(j6-j1)<=j5<=abs(j6+j1)&&j1+j2+j3==floor(j1+j2+j3)&&j3+1+j5==floor(j3+1+j5)&&j2+1+j6==floor(j2+1+j6)&&j1+j5+j6==floor(j1+j5+j6)
                a1[i1]=wigner6j(j1,j2,j3,1,j5,j6)
            end
        end
    end
    
    function Wj6(jl,jm,A)
        j1=jl[1]
        j2=jl[2]
        j3=jl[3]
        j5=jl[4]
        j6=jl[5]
        if abs(j2-j1)<=j3<=abs(j2+j1)&&abs(j3-1)<=j5<=abs(j3+1)&&abs(j2-1)<=j6<=abs(j2+1)&&abs(j6-j1)<=j5<=abs(j6+j1)&&j1+j2+j3==floor(j1+j2+j3)&&j3+1+j5==floor(j3+1+j5)&&j2+1+j6==floor(j2+1+j6)&&j1+j5+j6==floor(j1+j5+j6)
            i=jtoi6(jm,jl)
            if 1<= i <= size(A)[1]
                return Float64(A[i])
            end
        end
        return 0.0
    end

    function Wj1(jl,jm,A,ij1)
        if abs(jl[4])<=jl[1]&&abs(jl[5])<=jl[2]&&abs(jl[6])<=jl[3]
            jtoi(jm,[jl[1] jl[2] jl[3]],[jl[4] jl[5] jl[6]],ij1)
            s1=size(A)
            if 1<=ij1[1]<=s1[1]&&1<=ij1[2]<=s1[2]
                return Float64(A[ij1[1],ij1[2]])
            end
        end
        return 0.0
    end

    function Wj2(jl,jm,A)
        ij2=[0 0]
        if abs(jl[4])<=jl[1]&&abs(jl[5])<=jl[2]&&abs(jl[6])<=jl[3]
            jtoi(jm,[jl[1] jl[2] jl[3]],[jl[4] jl[5] jl[6]],ij2)
            s1=size(A)
            if 1<=ij2[1]<=s1[1]&&1<=ij2[2]<=s1[2]
                return Float64(A[ij2[1],ij2[2]])
            end
        end
        return 0.0
    end
    
    function Wj1_mod(jl,jm,Atl,a,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j)
        #rs1=Array{Int32,2}(undef, a, 2)
        #rs2=zeros(a)
        @inbounds i.=(((@view jl[1:a,1]).*2).+1).+((@view jl[1:a,2]).*2).*(2*(jm[1])+1).+((@view jl[1:a,3]).*2).*((2*(jm[1])+1)*(2*(jm[2])+1))
        @inbounds ja1.=(@view jl[1:a,1]).-abs.(@view jl[1:a,4])
        @inbounds ja2.=(@view jl[1:a,2]).-abs.(@view jl[1:a,5])
        @inbounds ja3.=(@view jl[1:a,3]).-abs.(@view jl[1:a,6])
        replace!(x -> x>=0 ? 1 : 0, ja1)
        replace!(x -> x>=0 ? 1 : 0, ja2)
        replace!(x -> x>=0 ? 1 : 0, ja3)
        @inbounds ja.=ja1.*ja2.*ja3
        @inbounds j.=((@view jl[1:a,1]).-(@view jl[1:a,4]).+1).+(((@view jl[1:a,1]).*2).+1).*((@view jl[1:a,2]).-(@view jl[1:a,5]))
        @inbounds ij1[1:a,1].=broadcast(Int32,i)
        @inbounds ij1[1:a,2].=broadcast(Int32,j)
        @inbounds rs1.=broadcast(Int32,(ij1.*ja).-(ja.-1))
        @inbounds rs11.=@view rs1[1:a,1]
        @inbounds rs12.=@view rs1[1:a,2]
        @inbounds rs2.=(broadcast((x,y)->Atl[x,y],rs11,rs12))
        @inbounds rs2.=rs2.*ja
        return 0
    end

    function Wj1_mod_backup(jl,jm,A,a,ij1)
        rs1=Array{Int32,2}(undef, a, 2)
        rs2=zeros(a)
        i=((view(jl,1:a,1).*2).+1).+(view(jl,1:a,2).*2).*(2*jm[1]+1).+(view(jl,a,3).*2).*((2*jm[1]+1)*(2*jm[2]+1))
        ja1=(view(jl,1:a,1).-broadcast(abs,view(jl,1:a,4)))
        ja2=(view(jl,1:a,2).-broadcast(abs,view(jl,1:a,5)))
        ja3=(view(jl,1:a,3).-broadcast(abs,view(jl,1:a,6)))
        replace!(x -> x>=0 ? 1 : 0, ja1)
        replace!(x -> x>=0 ? 1 : 0, ja2)
        replace!(x -> x>=0 ? 1 : 0, ja3)
        ja=ja1.*ja2.*ja3
        j=(view(jl,1:a,1).-view(jl,1:a,4).+1).+((view(jl,a,1).*2).+1).*(view(jl,a,2).-view(jl,a,5))
        ij1[1:a,1]=broadcast(Int32,i)
        ij1[1:a,2]=broadcast(Int32,j)
        rs1=broadcast(Int32,(ij1.*ja).-(ja.-1))
        #return A[view(rs1,1:a,1),view(rs1,1:a,2)].*ja
        rs2.=view(A,view(rs1,1:a,1:2))
        rs2.=rs2.*ja
        return rs2
    end