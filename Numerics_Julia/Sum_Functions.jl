function kronmc1(j,sj)
    m4=(sj[4]-1)/2-rem(j-1,sj[4])
    m3=(sj[3]-1)/2-rem(fld(j-1,sj[4]),sj[3])
    m2=(sj[2]-1)/2-rem(fld(fld(j-1,sj[4]),sj[3]),sj[2])
    m1=(sj[1]-1)/2-rem(fld(fld(fld(j-1,sj[4]),sj[3]),sj[2]),sj[1])
    return [m1 m2 m3 m4]
end

function kronmc3(j,sj)
    m3=(sj[3]-1)/2-rem(j-1,sj[3])
    m2=(sj[2]-1)/2-rem(fld(j-1,sj[3]),sj[2])
    m1=(sj[1]-1)/2-rem(fld(fld(j-1,sj[3]),sj[2]),sj[1])
    return [m1 m2 m3]
end

function transum(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        g1c=exp(-t1*(lb1(ls123[1])))*adjoint(glc(ls123[1],gf1))
        #println(g1c)
        g2c=exp(-t1*(lb1(ls123[2])))*adjoint(glc(ls123[2],gf2))
        g3c=exp(-t1*(lb1(ls123[3])))*adjoint(glc(ls123[3],gf3))
        g4c=exp(-t1*(lb1(ls123[4])))*adjoint(glc(ls123[4],gf4))
        ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)*sqrt(2*ls123[4]+1)
        kC1=zeros(Int(2*ls123[1]+1)*Int(2*ls123[2]+1))*im
        kC2=zeros(Int(2*ls123[1]+1)*Int(2*ls123[2]+1)*Int(2*ls123[3]+1))*im
        kC2=zeros(Int(2*ls123[1]+1)*Int(2*ls123[2]+1)*Int(2*ls123[3]+1)*Int(2*ls123[4]+1))*im
        state=zeros(Int(2*ls123[1]+1)*Int(2*ls123[2]+1)*Int(2*ls123[3]+1)*Int(2*ls123[4]+1),Int(2*ls123[1]+1)*Int(2*ls123[2]+1)*Int(2*ls123[3]+1)*Int(2*ls123[4]+1))*im
        for ma in 1:Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1))
            ma1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
            ma2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
            ma3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
            ma4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
            i1=Int(ls123[1]-ma1+1)
            i2=Int(ls123[2]-ma2+1)
            i3=Int(ls123[3]-ma3+1)
            i4=Int(ls123[4]-ma4+1)
            kron!(kC1,g1c[i1,:],g2c[i2,:])
            kron!(kC2,kC1,g3c[i3,:])
            kron!(kC3,kC2,g4c[i4,:])
            state[ma,:].=state[ma,:].+(kC3[:].*ic1)
        end
        for ma in 1:Int((2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1))
            ma1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
            ma2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
            ma3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
            ma4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
            i1=Int(ls123[1]-ma1+1)
            i2=Int(ls123[2]-ma2+1)
            i3=Int(ls123[3]-ma3+1)
            i4=Int(ls123[4]-ma4+1)
            v1=qo(ls123,[ma1 ma2 ma3 ma4],[mb1 mb2 mb3 mb4],1,2,3)-qo(ls123,[ma1 ma2 ma3 ma4],[mb1 mb2 mb3 mb4],1,2,4)+qo(ls123,[ma1 ma2 ma3 ma4],[mb1 mb2 mb3 mb4],1,3,4)-qo(ls123,[ma1 ma2 ma3 ma4],[mb1 mb2 mb3 mb4],2,3,4)
            at1[i1,2]=at1[i1,2]+adjoint(state[:,ma1])*v1*state[:,ma1]
        end
        state=zeros(1,1)
    end
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
    return 0
end

function tran123(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(a1[:,ja])*av1*(t1)^(3/2)*conj(a1[:,ja])
                    at1[i1,2]=at1[i1,2]+transpose(a1[:,ja])*conj(a1[:,ja])
                    for ii in 3:nn
                        at1[i1,ii]=at1[i1,ii]+transpose(a1[:,ja])*((vol1q*t1^3/8*im)^(ii-2))*conj(a1[:,ja])
                    end
                    at1[i1,nn+1]=at1[i1,nn+1]+transpose(a1[:,ja])*jl1*conj(a1[:,ja])*t1^2
                    at1[i1,nn+2]=at1[i1,nn+2]+transpose(a1[:,ja])*jl2*conj(a1[:,ja])*t1^2
                    at1[i1,nn+3]=at1[i1,nn+3]+transpose(a1[:,ja])*jl3*conj(a1[:,ja])*t1^2
                    at1[i1,nn+4]=at1[i1,nn+4]+transpose(a1[:,ja])*jl4*conj(a1[:,ja])*t1^2
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

@everywhere function tran1231(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    for ii in 3:nn
                        at1[i1,ii]=at1[i1,ii]+transpose(conj(a1[:,ja]))*(vol1q^(ii-2))*a1[:,ja]
                    end
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123test(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                for is in 1:a
                    at1[i1,1]=at1[i1,1]+adjoint(state[:,is])*av1*(t1)^(3/2)*state[:,is]
                    at1[i1,2]=at1[i1,2]+adjoint(state[:,is])*state[:,is]
                    at1[i1,3]=at1[i1,3]+adjoint(state[:,is])*(vol1q*t1^3/8*im)*state[:,is]
                    for ii in 4:nn
                        at1[i1,ii]=at1[i1,ii]+adjoint(state[:,is])*((vol1q*t1^3/8*im)^(2*(ii-3)))*state[:,is]
                    end
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123test_s(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott_s(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                for is in 1:a
                    at1[i1,1]=at1[i1,1]+adjoint(state[:,is])*av1*(t1)^(3/2)*state[:,is]
                    at1[i1,2]=at1[i1,2]+adjoint(state[:,is])*state[:,is]
                    at1[i1,3]=at1[i1,3]+adjoint(state[:,is])*(vol1q*t1^3/8*im)*state[:,is]
                    for ii in 4:nn
                        at1[i1,ii]=at1[i1,ii]+adjoint(state[:,is])*((vol1q*t1^3/8*im)^(2*(ii-3)))*state[:,is]
                    end
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123test_sj(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=[js1/2 js2/2 js3/2 js4/2]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott_s(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                for is in 1:a
                    at1[js1,1]=at1[js1,1]+adjoint(state[:,is])*av1*(t1)^(3/2)*state[:,is]
                    at1[js1,2]=at1[js1,2]+adjoint(state[:,is])*state[:,is]
                    at1[js1,3]=at1[js1,3]+adjoint(state[:,is])*(vol1q*t1^3/8*im)*state[:,is]
                    for ii in 4:nn
                        at1[js1,ii]=at1[js1,ii]+adjoint(state[:,is])*((vol1q*t1^3/8*im)^(2*(ii-3)))*state[:,is]
                    end
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[js1,:].=im*0.0
            end
        else
            at1[js1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    return 0
end

function tran123test_soverlaptest(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=zeros(4)
        ls123[1]=js1
        ls123[2]=js2
        ls123[3]=js3
        ls123[4]=js4
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott_s(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)*im/4
                (e1,ev)=eigen(vol1q)
                maxe1=maximum(real.(e1))
                maxe1loc=findall(x->x==maxe1, real.(e1))[1]
                resoverlap=zeros(size(e1,1)+1)
                av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                ist=1
                # print(transpose(inv(ev)[maxe1loc,:])*state[:,ist])
                for is in 1:a
                    for evi in 1:size(e1,1)
                        resoverlap[evi]=resoverlap[evi]+abs(transpose(inv(ev)[evi,:])*state[:,is])^2
                    end
                    resoverlap[size(e1,1)+1]=resoverlap[size(e1,1)+1]+adjoint(state[:,is])*state[:,is]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
        print(resoverlap)
        print(e1)
        print(ev)
    return (resoverlap,e1,ev)
end

function tran123test_soverlap(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=zeros(4)
        ls123[1]=js1
        ls123[2]=js2
        ls123[3]=js3
        ls123[4]=js4
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott_sNew(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)*im/4
                (e1,ev)=eigen(vol1q)
                maxe1=maximum(real.(e1))
                maxe1loc=findall(x->x==maxe1, real.(e1))[1]
                resoverlap=zeros(size(e1,1)+1)
                av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                ist=1
                # print(transpose(inv(ev)[maxe1loc,:])*state[:,ist])
                for is in 1:a
                    resoverlap[1:size(e1,1)].=resoverlap[1:size(e1,1)].+abs.(inv(ev)[:,:]*state[:,is]).^2
                    resoverlap[size(e1,1)+1]=resoverlap[size(e1,1)+1]+adjoint(state[:,is])*state[:,is]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0res1
        end
        print(ev)
    return (resoverlap,e1,ev)
end

function tranVtest_s666(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm61,jm62,B1,B2)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1),6)
    ls2=zeros((js1+1)*(js2+1)*(js3+1))
    seq6oMod(js1,js2,js3,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
        g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
        g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai6(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                vol1q=im/4*(v1232(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v1232(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v1232(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v1232(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v1232(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,2,3,4,6,jm61,jm62,B1,B2))
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                operator_v=av1*(t1)^(3/2)
                operator_q1=(vol1q*t1^3)
                operator_q2=(vol1q*t1^3)^(2)
                operator_q4=(vol1q*t1^3)^(2*2)
                operator_q6=(vol1q*t1^3)^(2*3)
                operator_q8=(vol1q*t1^3)^(2*4)
                operator_q10=(vol1q*t1^3)^(2*5)
                operator_q12=(vol1q*t1^3)^(2*6)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                lls3=Array{Float64,2}(undef,a,6)
                lls4=Array{Float64,2}(undef,a,6)
                # print(size(lls1))
                #Initialzation of ls
                lls1[:,1].=ls123[1]
                lls1[:,2].=ls123[2]
                lls1[:,3].=view(Al,:,2)
                lls2[:,1].=view(Al,:,2)
                lls2[:,2].=ls123[3]
                lls2[:,3].=view(Al,:,3)
                lls3[:,1].=view(Al,:,4)
                lls3[:,2].=ls123[4]
                lls3[:,3].=view(Al,:,3)
                lls4[:,1].=ls123[5]
                lls4[:,2].=ls123[6]
                lls4[:,3].=view(Al,1:a,4)
                i2=Array{Float64,1}(undef,a)
                ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                ij1=Array{Int32}(undef, a, 2)
                bi=Array{Float64,1}(undef,a)
                ci=Array{Float64,1}(undef,a)
                di=Array{Float64,1}(undef,a)
                ja1=Array{Float64,1}(undef,a)
                ja2=Array{Float64,1}(undef,a)
                ja3=Array{Float64,1}(undef,a)
                ja=Array{Float64,1}(undef,a)
                i=Array{Float64,1}(undef,a)
                j=Array{Float64,1}(undef,a)
                rs1=Array{Int32}(undef, a, 2)
                rs11=Array{Int32}(undef, a)
                rs12=Array{Int32}(undef, a)
                rs2=zeros(a)
                sA=sqrt.((Al[:,2].*2).+1).*sqrt.((Al[:,3].*2).+1).*sqrt.((Al[:,4].*2).+1)
                state=zeros(a)*im
                for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)
                    #Block One
                    m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                    m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                    m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                    m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                    m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                    s1 = -m1-m2
                    s2 = s1-m3
                    s3 = s2-m4
                    m6 = s3-m5
                    if abs(m6)<=ls123[6]
                    #Block Two    
                        lls1[1:a,4].=m1
                        lls1[1:a,5].=m2
                        lls1[1:a,6].=s1
                        lls2[1:a,4].=-s1
                        lls2[1:a,5].=m3  
                        lls2[1:a,6].=s2
                        lls3[1:a,4].=s3
                        lls3[1:a,5].=m4  
                        lls3[1:a,6].=-s2
                        lls4[1:a,4].=m5
                        lls4[1:a,5].=m6
                        lls4[1:a,6].=-s3
        #                 if j6+m6+1==1.5
        #                     print([j1 j2 j3 j4 j5 j6; m1 m2 m3 m4 m5 m6])
        #                 end
                    #Block Three
                        # print(ls123)
                        Inter6_mod2(lls1,lls2,lls3,lls4,aj1,jm,a,s1,s2,s3,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi, ci, di) 
                        com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                        state[:].=state[:].+(sA.*i2*com)
                    end
                end
                at1[i1,1]=at1[i1,1]+adjoint(state[:])*operator_v*state[:]
                at1[i1,2]=at1[i1,2]+adjoint(state[:])*state[:]
                at1[i1,3]=at1[i1,3]+adjoint(state[:])*operator_q1*state[:]
                at1[i1,4]=at1[i1,4]+adjoint(state[:])*operator_q2*state[:]
                at1[i1,5]=at1[i1,5]+adjoint(state[:])*operator_q4*state[:]
                at1[i1,6]=at1[i1,6]+adjoint(state[:])*operator_q6*state[:]
                at1[i1,7]=at1[i1,7]+adjoint(state[:])*operator_q8*state[:]
                at1[i1,8]=at1[i1,8]+adjoint(state[:])*operator_q10*state[:]
                at1[i1,9]=at1[i1,9]+adjoint(state[:])*operator_q12*state[:]
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123gv7New(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls11=Vector{Matrix{Float64}}()
    ls21=[]
    seq7oModN(js1,js2,js3,ls11,ls21)
    @everywhere ls1=$ls11
    @everywhere ls2=$ls21
    @sync @distributed for i1nn in 1:size(ls1,1)
            ls123=ls1[Int(ls2[i1nn])][:]        
            if ls123[7]<=ls123[1]+ls123[2]+ls123[3]+ls123[1]+ls123[2]+ls123[3]
                g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
                g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
                g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
                Al=Fai7(ls123)
                a=size(Al,1)
                if size(Al)[1]!=0
                        ov123=v1232(ls123,Al,1,2,3,7,jm61,jm62,B1,B2)*t1^3
                        ov135=v1232(ls123,Al,1,3,5,7,jm61,jm62,B1,B2)*t1^3
                        ov156=v1232(ls123,Al,1,5,6,7,jm61,jm62,B1,B2)*t1^3
                        ov126=v1232(ls123,Al,1,2,6,7,jm61,jm62,B1,B2)*t1^3
                        ov345=v1232(ls123,Al,3,4,5,7,jm61,jm62,B1,B2)*t1^3
                        ov456=v1232(ls123,Al,4,5,6,7,jm61,jm62,B1,B2)*t1^3
                        ov246=v1232(ls123,Al,2,4,6,7,jm61,jm62,B1,B2)*t1^3
                        ov234=v1232(ls123,Al,2,3,4,7,jm61,jm62,B1,B2)*t1^3
                        # vol1q=im/4*(ov123+ov135+ov156-ov126+ov345-ov456-ov246-ov234)
                        # (e1,ev)=eigen(vol1q)
                        # av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                        # operator_v=av1*(t1)^(3/2)
                        # operator_q1=(vol1q*t1^3)
                        # operator_q2=(vol1q*t1^3)^(2)
                        # operator_q4=(vol1q*t1^3)^(2*2)
                        # operator_q6=(vol1q*t1^3)^(2*3)
                        # operator_q8=(vol1q*t1^3)^(2*4)
                        # operator_q10=(vol1q*t1^3)^(2*5)
                        # operator_q12=(vol1q*t1^3)^(2*6)
                        # lls1=Array{Float64,2}(undef,a,6)
                        # lls2=Array{Float64,2}(undef,a,6)
                        # lls3=Array{Float64,2}(undef,a,6)
                        # lls4=Array{Float64,2}(undef,a,6)
                        # lls5=Array{Float64,2}(undef,a,6)
                        # ij1=Array{Int32}(undef, a, 2)
                        # lls1[1:a,1].=ls123[1]
                        # lls1[1:a,2].=ls123[2]
                        # lls1[1:a,3].=view(Al,1:a,2)
                        # lls2[1:a,1].=view(Al,1:a,2)
                        # lls2[1:a,2].=ls123[3]
                        # lls2[1:a,3].=view(Al,1:a,3)
                        # lls3[1:a,1].=view(Al,1:a,4)
                        # lls3[1:a,2].=ls123[4]
                        # lls3[1:a,3].=view(Al,1:a,3)
                        # lls4[1:a,1].=view(Al,1:a,5)
                        # lls4[1:a,2].=ls123[5]
                        # lls4[1:a,3].=view(Al,1:a,4)
                        # lls5[1:a,1].=ls123[6]
                        # lls5[1:a,2].=ls123[7]
                        # lls5[1:a,3].=view(Al,1:a,5)
                        # i2=Array{Float64,1}(undef,a)
                        # ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                        # ij1=Array{Int32}(undef, a, 2)
                        # bi=Array{Float64,1}(undef,a)
                        # ci=Array{Float64,1}(undef,a)
                        # di=Array{Float64,1}(undef,a)
                        # ei=Array{Float64,1}(undef,a)
                        # ja1=Array{Float64,1}(undef,a)
                        # ja2=Array{Float64,1}(undef,a)
                        # ja3=Array{Float64,1}(undef,a)
                        # ja=Array{Float64,1}(undef,a)
                        # i=Array{Float64,1}(undef,a)
                        # j=Array{Float64,1}(undef,a)
                        # rs1=Array{Int32}(undef, a, 2)
                        # rs11=Array{Int32}(undef, a)
                        # rs12=Array{Int32}(undef, a)
                        # rs2=zeros(a)
                        # ssum1=zeros(a,a)*im
                        # sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                        # state=zeros(a,Int32(2*ls123[7]+1))*im
                        # for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                        #     m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                        #     m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                        #     m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                        #     m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                        #     m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                        #     m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                        #     s1 = -m1-m2
                        #     s2 = s1-m3
                        #     s3 = s2-m4
                        #     s4 = s3-m5
                        #     m7 = s4-m6
                        #     if abs(m7)<=ls123[7]
                        #         @view(lls1[1:a,4]).=m1
                        #         @view(lls1[1:a,5]).=m2
                        #         @view(lls1[1:a,6]).=s1
                        #         @view(lls2[1:a,4]).=-s1
                        #         @view(lls2[1:a,5]).=m3
                        #         @view(lls2[1:a,6]).=s2
                        #         @view(lls3[1:a,4]).=s3
                        #         @view(lls3[1:a,5]).=m4  
                        #         @view(lls3[1:a,6]).=-s2
                        #         @view(lls4[1:a,4]).=s4
                        #         @view(lls4[1:a,5]).=m5  
                        #         @view(lls4[1:a,6]).=-s3
                        #         @view(lls5[1:a,4]).=m6
                        #         @view(lls5[1:a,5]).=m7
                        #         @view(lls5[1:a,6]).=-s4
                        #         # print(ls123,[m1,m2,m3,m4,m5,m6])
                        #         com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                        #         Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                        #         inter5n=i2[:].*sA*com
                        #         state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                        #     end
                        # end
                        # for lis1 in 1:Int(2*ls123[7]+1)
                        #     ssum1.=ssum1.+conj(state[:,lis1])*transpose(state[:,lis1])
                        # end
    #                     for lisi in 1:a
    #                         for lisj in 1:a
    #                             at1[i1nn,1]=at1[i1nn,1]+operator_v[lisi,lisj]*ssum1[lisi,lisj]
    #                             if lisi==lisj
    #                                 at1[i1nn,2]=at1[i1nn,2]+ssum2[lisi,lisj]
    #                             end
    #                             at1[i1nn,3]=at1[i1nn,3]+operator_q1[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,4]=at1[i1nn,4]+operator_q2[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,5]=at1[i1nn,5]+operator_q4[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,6]=at1[i1nn,6]+ov123[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,7]=at1[i1nn,7]+ov135[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,8]=at1[i1nn,8]+ov156[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,9]=at1[i1nn,9]+ov126[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,10]=at1[i1nn,10]+ov345[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,11]=at1[i1nn,11]+ov456[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,12]=at1[i1nn,12]+ov246[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1nn,13]=at1[i1nn,13]+ov234[lisi,lisj]*ssum1[lisi,lisj]
    # #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
    # #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
    # #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
    # #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
    #                         end
    #                     end
                        # for mj1 in 1:Int32(2*ls123[7]+1)
                        #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                        #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                        #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                        #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                        #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                        #     at1[i1nn,10]=at1[i1nn,10]+adjoint(state[:,mj1])*ov123*state[:,mj1]
                        #     at1[i1nn,11]=at1[i1nn,11]+adjoint(state[:,mj1])*ov135*state[:,mj1]
                        #     at1[i1nn,12]=at1[i1nn,12]+adjoint(state[:,mj1])*ov156*state[:,mj1]
                        #     at1[i1nn,13]=at1[i1nn,13]+adjoint(state[:,mj1])*ov126*state[:,mj1]
                        #     at1[i1nn,14]=at1[i1nn,14]+adjoint(state[:,mj1])*ov345*state[:,mj1]
                        #     at1[i1nn,15]=at1[i1nn,15]+adjoint(state[:,mj1])*ov456*state[:,mj1]
                        #     at1[i1nn,16]=at1[i1nn,16]+adjoint(state[:,mj1])*ov246*state[:,mj1]
                        #     at1[i1nn,17]=at1[i1nn,17]+adjoint(state[:,mj1])*ov234*state[:,mj1]
                        # end
                        state=nothing
                else
                    at1[i1nn,:].=at1[i1nn,:].+im*0.0
                end
                Al=nothing
                a=nothing
            end
    end
    return 0
end

function tran123gv7New2(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls11=Vector{Matrix{Float64}}()
    ls21=[]
    seq7oModN(js1,js2,js3,ls11,ls21)
    @everywhere ls1=$ls11
    @everywhere ls2=$ls21
    @sync @distributed for i1nn in 1:size(ls1,1)
            ls123=ls1[Int(ls2[i1nn])][:]      
            if ls123[7]<=ls123[1]+ls123[2]+ls123[3]+ls123[1]+ls123[2]+ls123[3]
                g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
                g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
                g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
                Al=Fai7(ls123)
                a=size(Al,1)
                if size(Al)[1]!=0
                        ov123=v1232(ls123,Al,1,2,3,7,jm61,jm62,A1,B1)*t1^3
                        ov135=v1232(ls123,Al,1,3,5,7,jm61,jm62,A1,B1)*t1^3
                        ov156=v1232(ls123,Al,1,5,6,7,jm61,jm62,A1,B1)*t1^3
                        ov126=v1232(ls123,Al,1,2,6,7,jm61,jm62,A1,B1)*t1^3
                        ov345=v1232(ls123,Al,3,4,5,7,jm61,jm62,A1,B1)*t1^3
                        ov456=v1232(ls123,Al,4,5,6,7,jm61,jm62,A1,B1)*t1^3
                        ov246=v1232(ls123,Al,2,4,6,7,jm61,jm62,A1,B1)*t1^3
                        ov234=v1232(ls123,Al,2,3,4,7,jm61,jm62,A1,B1)*t1^3
                        vol1q=im/4*(ov123+ov135+ov156-ov126+ov345-ov456-ov246-ov234)
                        vol1qM=Matrix(vol1q)
                        (e1,ev)=eigen(vol1qM)
                        av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                        vol1qM=nothing
                        operator_v=av1
                        operator_q1=(vol1q)
                        operator_q2=(vol1q)^(2)
                        operator_q4=(vol1q)^(2*2)
                        lls1=Array{Float64,2}(undef,a,6)
                        lls2=Array{Float64,2}(undef,a,6)
                        lls3=Array{Float64,2}(undef,a,6)
                        lls4=Array{Float64,2}(undef,a,6)
                        lls5=Array{Float64,2}(undef,a,6)
                        ij1=Array{Int32}(undef, a, 2)
                        lls1[1:a,1].=ls123[1]
                        lls1[1:a,2].=ls123[2]
                        lls1[1:a,3].=view(Al,1:a,2)
                        lls2[1:a,1].=view(Al,1:a,2)
                        lls2[1:a,2].=ls123[3]
                        lls2[1:a,3].=view(Al,1:a,3)
                        lls3[1:a,1].=view(Al,1:a,4)
                        lls3[1:a,2].=ls123[4]
                        lls3[1:a,3].=view(Al,1:a,3)
                        lls4[1:a,1].=view(Al,1:a,5)
                        lls4[1:a,2].=ls123[5]
                        lls4[1:a,3].=view(Al,1:a,4)
                        lls5[1:a,1].=ls123[6]
                        lls5[1:a,2].=ls123[7]
                        lls5[1:a,3].=view(Al,1:a,5)
                        i2=Array{Float64,1}(undef,a)
                        ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                        ij1=Array{Int32}(undef, a, 2)
                        bi=Array{Float64,1}(undef,a)
                        ci=Array{Float64,1}(undef,a)
                        di=Array{Float64,1}(undef,a)
                        ei=Array{Float64,1}(undef,a)
                        ja1=Array{Float64,1}(undef,a)
                        ja2=Array{Float64,1}(undef,a)
                        ja3=Array{Float64,1}(undef,a)
                        ja=Array{Float64,1}(undef,a)
                        i=Array{Float64,1}(undef,a)
                        j=Array{Float64,1}(undef,a)
                        rs1=Array{Int32}(undef, a, 2)
                        rs11=Array{Int32}(undef, a)
                        rs12=Array{Int32}(undef, a)
                        rs2=zeros(a)
                        ssum1=zeros(a,a)*im
                        sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                        state=zeros(a,Int32(2*ls123[7]+1))*im
                        for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                            m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                            m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                            m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                            m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                            m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                            m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                            s1 = -m1-m2
                            s2 = s1-m3
                            s3 = s2-m4
                            s4 = s3-m5
                            m7 = s4-m6
                            if abs(m7)<=ls123[7]
                                @view(lls1[1:a,4]).=m1
                                @view(lls1[1:a,5]).=m2
                                @view(lls1[1:a,6]).=s1
                                @view(lls2[1:a,4]).=-s1
                                @view(lls2[1:a,5]).=m3
                                @view(lls2[1:a,6]).=s2
                                @view(lls3[1:a,4]).=s3
                                @view(lls3[1:a,5]).=m4  
                                @view(lls3[1:a,6]).=-s2
                                @view(lls4[1:a,4]).=s4
                                @view(lls4[1:a,5]).=m5  
                                @view(lls4[1:a,6]).=-s3
                                @view(lls5[1:a,4]).=m6
                                @view(lls5[1:a,5]).=m7
                                @view(lls5[1:a,6]).=-s4
                                # print(ls123,[m1,m2,m3,m4,m5,m6])
                                com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                                Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                                inter5n=i2[:].*sA*com
                                state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                            end
                        end
                        for lis1 in 1:Int(2*ls123[7]+1)
                            ssum1.=ssum1.+conj(state[:,lis1])*transpose(state[:,lis1])
                        end
                        for lisi in 1:a
                            for lisj in 1:a
                                at1[i1nn,1]=at1[i1nn,1]+operator_v[lisi,lisj]*ssum1[lisi,lisj]
                                if lisi==lisj
                                    at1[i1nn,2]=at1[i1nn,2]+ssum1[lisi,lisj]
                                end
                                at1[i1nn,3]=at1[i1nn,3]+operator_q1[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,4]=at1[i1nn,4]+operator_q2[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,5]=at1[i1nn,5]+operator_q4[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,6]=at1[i1nn,6]+ov123[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,7]=at1[i1nn,7]+ov135[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,8]=at1[i1nn,8]+ov156[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,9]=at1[i1nn,9]+ov126[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,10]=at1[i1nn,10]+ov345[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,11]=at1[i1nn,11]+ov456[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,12]=at1[i1nn,12]+ov246[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,13]=at1[i1nn,13]+ov234[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                            end
                        end
                        # for mj1 in 1:Int32(2*ls123[7]+1)
                        #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                        #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                        #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                        #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                        #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                        #     at1[i1nn,10]=at1[i1nn,10]+adjoint(state[:,mj1])*ov123*state[:,mj1]
                        #     at1[i1nn,11]=at1[i1nn,11]+adjoint(state[:,mj1])*ov135*state[:,mj1]
                        #     at1[i1nn,12]=at1[i1nn,12]+adjoint(state[:,mj1])*ov156*state[:,mj1]
                        #     at1[i1nn,13]=at1[i1nn,13]+adjoint(state[:,mj1])*ov126*state[:,mj1]
                        #     at1[i1nn,14]=at1[i1nn,14]+adjoint(state[:,mj1])*ov345*state[:,mj1]
                        #     at1[i1nn,15]=at1[i1nn,15]+adjoint(state[:,mj1])*ov456*state[:,mj1]
                        #     at1[i1nn,16]=at1[i1nn,16]+adjoint(state[:,mj1])*ov246*state[:,mj1]
                        #     at1[i1nn,17]=at1[i1nn,17]+adjoint(state[:,mj1])*ov234*state[:,mj1]
                        # end
                        state=nothing
                else
                    at1[i1nn,:].=at1[i1nn,:].+im*0.0
                end
                Al=nothing
                a=nothing
            end
    end
    return 0
end

function tran123gv7New2test(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
            ls123=[10,10,10,10,10,10,10]
            if ls123[7]<=ls123[1]+ls123[2]+ls123[3]+ls123[1]+ls123[2]+ls123[3]
                g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
                g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
                g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
                Al=Fai7(ls123)
                a=size(Al,1)
                if size(Al)[1]!=0
                        ov123=v1232(ls123,Al,1,2,3,7,jm61,jm62,A1,B1)*t1^3
                        ov135=v1232(ls123,Al,1,3,5,7,jm61,jm62,A1,B1)*t1^3
                        ov156=v1232(ls123,Al,1,5,6,7,jm61,jm62,A1,B1)*t1^3
                        ov126=v1232(ls123,Al,1,2,6,7,jm61,jm62,A1,B1)*t1^3
                        ov345=v1232(ls123,Al,3,4,5,7,jm61,jm62,A1,B1)*t1^3
                        ov456=v1232(ls123,Al,4,5,6,7,jm61,jm62,A1,B1)*t1^3
                        ov246=v1232(ls123,Al,2,4,6,7,jm61,jm62,A1,B1)*t1^3
                        ov234=v1232(ls123,Al,2,3,4,7,jm61,jm62,A1,B1)*t1^3
                        vol1q=im/4*(ov123+ov135+ov156-ov126+ov345-ov456-ov246-ov234)
                        vol1qM=vol1q
                        (e1,ev)=eigs(vol1qM)
                        av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                        vol1qM=nothing
                        operator_v=av1
                        operator_q1=(vol1q)
                        operator_q2=(vol1q)^(2)
                        operator_q4=(vol1q)^(2*2)
                        lls1=Array{Float64,2}(undef,a,6)
                        lls2=Array{Float64,2}(undef,a,6)
                        lls3=Array{Float64,2}(undef,a,6)
                        lls4=Array{Float64,2}(undef,a,6)
                        lls5=Array{Float64,2}(undef,a,6)
                        ij1=Array{Int32}(undef, a, 2)
                        lls1[1:a,1].=ls123[1]
                        lls1[1:a,2].=ls123[2]
                        lls1[1:a,3].=view(Al,1:a,2)
                        lls2[1:a,1].=view(Al,1:a,2)
                        lls2[1:a,2].=ls123[3]
                        lls2[1:a,3].=view(Al,1:a,3)
                        lls3[1:a,1].=view(Al,1:a,4)
                        lls3[1:a,2].=ls123[4]
                        lls3[1:a,3].=view(Al,1:a,3)
                        lls4[1:a,1].=view(Al,1:a,5)
                        lls4[1:a,2].=ls123[5]
                        lls4[1:a,3].=view(Al,1:a,4)
                        lls5[1:a,1].=ls123[6]
                        lls5[1:a,2].=ls123[7]
                        lls5[1:a,3].=view(Al,1:a,5)
                        i2=Array{Float64,1}(undef,a)
                        ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                        ij1=Array{Int32}(undef, a, 2)
                        bi=Array{Float64,1}(undef,a)
                        ci=Array{Float64,1}(undef,a)
                        di=Array{Float64,1}(undef,a)
                        ei=Array{Float64,1}(undef,a)
                        ja1=Array{Float64,1}(undef,a)
                        ja2=Array{Float64,1}(undef,a)
                        ja3=Array{Float64,1}(undef,a)
                        ja=Array{Float64,1}(undef,a)
                        i=Array{Float64,1}(undef,a)
                        j=Array{Float64,1}(undef,a)
                        rs1=Array{Int32}(undef, a, 2)
                        rs11=Array{Int32}(undef, a)
                        rs12=Array{Int32}(undef, a)
                        rs2=zeros(a)
                        ssum1=zeros(a,a)*im
                        sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                        state=zeros(a,Int32(2*ls123[7]+1))*im
                        for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                            m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                            m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                            m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                            m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                            m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                            m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                            s1 = -m1-m2
                            s2 = s1-m3
                            s3 = s2-m4
                            s4 = s3-m5
                            m7 = s4-m6
                            if abs(m7)<=ls123[7]
                                @view(lls1[1:a,4]).=m1
                                @view(lls1[1:a,5]).=m2
                                @view(lls1[1:a,6]).=s1
                                @view(lls2[1:a,4]).=-s1
                                @view(lls2[1:a,5]).=m3
                                @view(lls2[1:a,6]).=s2
                                @view(lls3[1:a,4]).=s3
                                @view(lls3[1:a,5]).=m4  
                                @view(lls3[1:a,6]).=-s2
                                @view(lls4[1:a,4]).=s4
                                @view(lls4[1:a,5]).=m5  
                                @view(lls4[1:a,6]).=-s3
                                @view(lls5[1:a,4]).=m6
                                @view(lls5[1:a,5]).=m7
                                @view(lls5[1:a,6]).=-s4
                                # print(ls123,[m1,m2,m3,m4,m5,m6])
                                com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                                Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                                inter5n=i2[:].*sA*com
                                state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                            end
                        end
                        for lis1 in 1:Int(2*ls123[7]+1)
                            ssum1.=ssum1.+conj(state[:,lis1])*transpose(state[:,lis1])
                        end
                        for lisi in 1:a
                            for lisj in 1:a
                                at1[i1nn,1]=at1[i1nn,1]+operator_v[lisi,lisj]*ssum1[lisi,lisj]
                                if lisi==lisj
                                    at1[i1nn,2]=at1[i1nn,2]+ssum1[lisi,lisj]
                                end
                                at1[i1nn,3]=at1[i1nn,3]+operator_q1[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,4]=at1[i1nn,4]+operator_q2[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,5]=at1[i1nn,5]+operator_q4[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,6]=at1[i1nn,6]+ov123[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,7]=at1[i1nn,7]+ov135[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,8]=at1[i1nn,8]+ov156[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,9]=at1[i1nn,9]+ov126[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,10]=at1[i1nn,10]+ov345[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,11]=at1[i1nn,11]+ov456[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,12]=at1[i1nn,12]+ov246[lisi,lisj]*ssum1[lisi,lisj]
                                at1[i1nn,13]=at1[i1nn,13]+ov234[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                            end
                        end
                        # for mj1 in 1:Int32(2*ls123[7]+1)
                        #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                        #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                        #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                        #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                        #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                        #     at1[i1nn,10]=at1[i1nn,10]+adjoint(state[:,mj1])*ov123*state[:,mj1]
                        #     at1[i1nn,11]=at1[i1nn,11]+adjoint(state[:,mj1])*ov135*state[:,mj1]
                        #     at1[i1nn,12]=at1[i1nn,12]+adjoint(state[:,mj1])*ov156*state[:,mj1]
                        #     at1[i1nn,13]=at1[i1nn,13]+adjoint(state[:,mj1])*ov126*state[:,mj1]
                        #     at1[i1nn,14]=at1[i1nn,14]+adjoint(state[:,mj1])*ov345*state[:,mj1]
                        #     at1[i1nn,15]=at1[i1nn,15]+adjoint(state[:,mj1])*ov456*state[:,mj1]
                        #     at1[i1nn,16]=at1[i1nn,16]+adjoint(state[:,mj1])*ov246*state[:,mj1]
                        #     at1[i1nn,17]=at1[i1nn,17]+adjoint(state[:,mj1])*ov234*state[:,mj1]
                        # end
                        state=nothing
                else
                    at1[i1nn,:].=at1[i1nn,:].+im*0.0
                end
                Al=nothing
                a=nothing
            end
    return 0
end

function tran123gv7New2test2(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls11=Vector{Matrix{Float64}}()
    ls21=[]
    seq7oModN(js1,js2,js3,ls11,ls21)
    @everywhere ls1=$ls11
    @everywhere ls2=$ls21
    @sync @distributed for i1nn in 1:size(ls1,1)
        ls123=ls1[Int(ls2[i1nn])][:]      
        if ls123[7]<=ls123[1]+ls123[2]+ls123[3]+ls123[1]+ls123[2]+ls123[3]
            g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
            g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
            g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
            Al=Fai7(ls123)
            a=size(Al,1)
            if size(Al)[1]!=0
                if size(Al)>5000
                    ov123=v1232(ls123,Al,1,2,3,7,jm61,jm62,A1,B1)*t1^3
                    ov135=v1232(ls123,Al,1,3,5,7,jm61,jm62,A1,B1)*t1^3
                    ov156=v1232(ls123,Al,1,5,6,7,jm61,jm62,A1,B1)*t1^3
                    ov126=v1232(ls123,Al,1,2,6,7,jm61,jm62,A1,B1)*t1^3
                    ov345=v1232(ls123,Al,3,4,5,7,jm61,jm62,A1,B1)*t1^3
                    ov456=v1232(ls123,Al,4,5,6,7,jm61,jm62,A1,B1)*t1^3
                    ov246=v1232(ls123,Al,2,4,6,7,jm61,jm62,A1,B1)*t1^3
                    ov234=v1232(ls123,Al,2,3,4,7,jm61,jm62,A1,B1)*t1^3
                    vol1q=im/4*(ov123+ov135+ov156-ov126+ov345-ov456-ov246-ov234)
                    vol1qM=vol1q
                    (e1,ev)=eigs(vol1qM)
                    av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                    vol1qM=nothing
                    operator_v=av1
                    operator_q1=(vol1q)
                    operator_q2=(vol1q)^(2)
                    operator_q4=(vol1q)^(2*2)
                    lls1=Array{Float64,2}(undef,a,6)
                    lls2=Array{Float64,2}(undef,a,6)
                    lls3=Array{Float64,2}(undef,a,6)
                    lls4=Array{Float64,2}(undef,a,6)
                    lls5=Array{Float64,2}(undef,a,6)
                    ij1=Array{Int32}(undef, a, 2)
                    lls1[1:a,1].=ls123[1]
                    lls1[1:a,2].=ls123[2]
                    lls1[1:a,3].=view(Al,1:a,2)
                    lls2[1:a,1].=view(Al,1:a,2)
                    lls2[1:a,2].=ls123[3]
                    lls2[1:a,3].=view(Al,1:a,3)
                    lls3[1:a,1].=view(Al,1:a,4)
                    lls3[1:a,2].=ls123[4]
                    lls3[1:a,3].=view(Al,1:a,3)
                    lls4[1:a,1].=view(Al,1:a,5)
                    lls4[1:a,2].=ls123[5]
                    lls4[1:a,3].=view(Al,1:a,4)
                    lls5[1:a,1].=ls123[6]
                    lls5[1:a,2].=ls123[7]
                    lls5[1:a,3].=view(Al,1:a,5)
                    i2=Array{Float64,1}(undef,a)
                    ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                    ij1=Array{Int32}(undef, a, 2)
                    bi=Array{Float64,1}(undef,a)
                    ci=Array{Float64,1}(undef,a)
                    di=Array{Float64,1}(undef,a)
                    ei=Array{Float64,1}(undef,a)
                    ja1=Array{Float64,1}(undef,a)
                    ja2=Array{Float64,1}(undef,a)
                    ja3=Array{Float64,1}(undef,a)
                    ja=Array{Float64,1}(undef,a)
                    i=Array{Float64,1}(undef,a)
                    j=Array{Float64,1}(undef,a)
                    rs1=Array{Int32}(undef, a, 2)
                    rs11=Array{Int32}(undef, a)
                    rs12=Array{Int32}(undef, a)
                    rs2=zeros(a)
                    ssum1=zeros(a,a)*im
                    sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                    state=zeros(a,Int32(2*ls123[7]+1))*im
                    for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                        m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                        m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                        m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                        m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                        m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                        m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                        s1 = -m1-m2
                        s2 = s1-m3
                        s3 = s2-m4
                        s4 = s3-m5
                        m7 = s4-m6
                        if abs(m7)<=ls123[7]
                            @view(lls1[1:a,4]).=m1
                            @view(lls1[1:a,5]).=m2
                            @view(lls1[1:a,6]).=s1
                            @view(lls2[1:a,4]).=-s1
                            @view(lls2[1:a,5]).=m3
                            @view(lls2[1:a,6]).=s2
                            @view(lls3[1:a,4]).=s3
                            @view(lls3[1:a,5]).=m4  
                            @view(lls3[1:a,6]).=-s2
                            @view(lls4[1:a,4]).=s4
                            @view(lls4[1:a,5]).=m5  
                            @view(lls4[1:a,6]).=-s3
                            @view(lls5[1:a,4]).=m6
                            @view(lls5[1:a,5]).=m7
                            @view(lls5[1:a,6]).=-s4
                            # print(ls123,[m1,m2,m3,m4,m5,m6])
                            com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                            Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                            inter5n=i2[:].*sA*com
                            state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                        end
                    end
                    for lis1 in 1:Int(2*ls123[7]+1)
                        ssum1.=ssum1.+conj(state[:,lis1])*transpose(state[:,lis1])
                    end
                    for lisi in 1:a
                        for lisj in 1:a
                            at1[i1nn,1]=at1[i1nn,1]+operator_v[lisi,lisj]*ssum1[lisi,lisj]
                            if lisi==lisj
                                at1[i1nn,2]=at1[i1nn,2]+ssum1[lisi,lisj]
                            end
                            at1[i1nn,3]=at1[i1nn,3]+operator_q1[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,4]=at1[i1nn,4]+operator_q2[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,5]=at1[i1nn,5]+operator_q4[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,6]=at1[i1nn,6]+ov123[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,7]=at1[i1nn,7]+ov135[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,8]=at1[i1nn,8]+ov156[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,9]=at1[i1nn,9]+ov126[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,10]=at1[i1nn,10]+ov345[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,11]=at1[i1nn,11]+ov456[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,12]=at1[i1nn,12]+ov246[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,13]=at1[i1nn,13]+ov234[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                        end
                    end
                    # for mj1 in 1:Int32(2*ls123[7]+1)
                    #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                    #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                    #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                    #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                    #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                    #     at1[i1nn,10]=at1[i1nn,10]+adjoint(state[:,mj1])*ov123*state[:,mj1]
                    #     at1[i1nn,11]=at1[i1nn,11]+adjoint(state[:,mj1])*ov135*state[:,mj1]
                    #     at1[i1nn,12]=at1[i1nn,12]+adjoint(state[:,mj1])*ov156*state[:,mj1]
                    #     at1[i1nn,13]=at1[i1nn,13]+adjoint(state[:,mj1])*ov126*state[:,mj1]
                    #     at1[i1nn,14]=at1[i1nn,14]+adjoint(state[:,mj1])*ov345*state[:,mj1]
                    #     at1[i1nn,15]=at1[i1nn,15]+adjoint(state[:,mj1])*ov456*state[:,mj1]
                    #     at1[i1nn,16]=at1[i1nn,16]+adjoint(state[:,mj1])*ov246*state[:,mj1]
                    #     at1[i1nn,17]=at1[i1nn,17]+adjoint(state[:,mj1])*ov234*state[:,mj1]
                    # end
                    state=nothing
                else
                    ov123=v1232(ls123,Al,1,2,3,7,jm61,jm62,A1,B1)*t1^3
                    ov135=v1232(ls123,Al,1,3,5,7,jm61,jm62,A1,B1)*t1^3
                    ov156=v1232(ls123,Al,1,5,6,7,jm61,jm62,A1,B1)*t1^3
                    ov126=v1232(ls123,Al,1,2,6,7,jm61,jm62,A1,B1)*t1^3
                    ov345=v1232(ls123,Al,3,4,5,7,jm61,jm62,A1,B1)*t1^3
                    ov456=v1232(ls123,Al,4,5,6,7,jm61,jm62,A1,B1)*t1^3
                    ov246=v1232(ls123,Al,2,4,6,7,jm61,jm62,A1,B1)*t1^3
                    ov234=v1232(ls123,Al,2,3,4,7,jm61,jm62,A1,B1)*t1^3
                    vol1q=im/4*(ov123+ov135+ov156-ov126+ov345-ov456-ov246-ov234)
                    vol1qM=Matrix(vol1q)
                    (e1,ev)=eigen(vol1qM)
                    av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                    vol1qM=nothing
                    operator_v=av1
                    operator_q1=(vol1q)
                    operator_q2=(vol1q)^(2)
                    operator_q4=(vol1q)^(2*2)
                    lls1=Array{Float64,2}(undef,a,6)
                    lls2=Array{Float64,2}(undef,a,6)
                    lls3=Array{Float64,2}(undef,a,6)
                    lls4=Array{Float64,2}(undef,a,6)
                    lls5=Array{Float64,2}(undef,a,6)
                    ij1=Array{Int32}(undef, a, 2)
                    lls1[1:a,1].=ls123[1]
                    lls1[1:a,2].=ls123[2]
                    lls1[1:a,3].=view(Al,1:a,2)
                    lls2[1:a,1].=view(Al,1:a,2)
                    lls2[1:a,2].=ls123[3]
                    lls2[1:a,3].=view(Al,1:a,3)
                    lls3[1:a,1].=view(Al,1:a,4)
                    lls3[1:a,2].=ls123[4]
                    lls3[1:a,3].=view(Al,1:a,3)
                    lls4[1:a,1].=view(Al,1:a,5)
                    lls4[1:a,2].=ls123[5]
                    lls4[1:a,3].=view(Al,1:a,4)
                    lls5[1:a,1].=ls123[6]
                    lls5[1:a,2].=ls123[7]
                    lls5[1:a,3].=view(Al,1:a,5)
                    i2=Array{Float64,1}(undef,a)
                    ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                    ij1=Array{Int32}(undef, a, 2)
                    bi=Array{Float64,1}(undef,a)
                    ci=Array{Float64,1}(undef,a)
                    di=Array{Float64,1}(undef,a)
                    ei=Array{Float64,1}(undef,a)
                    ja1=Array{Float64,1}(undef,a)
                    ja2=Array{Float64,1}(undef,a)
                    ja3=Array{Float64,1}(undef,a)
                    ja=Array{Float64,1}(undef,a)
                    i=Array{Float64,1}(undef,a)
                    j=Array{Float64,1}(undef,a)
                    rs1=Array{Int32}(undef, a, 2)
                    rs11=Array{Int32}(undef, a)
                    rs12=Array{Int32}(undef, a)
                    rs2=zeros(a)
                    ssum1=zeros(a,a)*im
                    sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                    state=zeros(a,Int32(2*ls123[7]+1))*im
                    for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                        m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                        m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                        m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                        m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                        m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                        m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                        s1 = -m1-m2
                        s2 = s1-m3
                        s3 = s2-m4
                        s4 = s3-m5
                        m7 = s4-m6
                        if abs(m7)<=ls123[7]
                            @view(lls1[1:a,4]).=m1
                            @view(lls1[1:a,5]).=m2
                            @view(lls1[1:a,6]).=s1
                            @view(lls2[1:a,4]).=-s1
                            @view(lls2[1:a,5]).=m3
                            @view(lls2[1:a,6]).=s2
                            @view(lls3[1:a,4]).=s3
                            @view(lls3[1:a,5]).=m4  
                            @view(lls3[1:a,6]).=-s2
                            @view(lls4[1:a,4]).=s4
                            @view(lls4[1:a,5]).=m5  
                            @view(lls4[1:a,6]).=-s3
                            @view(lls5[1:a,4]).=m6
                            @view(lls5[1:a,5]).=m7
                            @view(lls5[1:a,6]).=-s4
                            # print(ls123,[m1,m2,m3,m4,m5,m6])
                            com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                            Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                            inter5n=i2[:].*sA*com
                            state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                        end
                    end
                    for lis1 in 1:Int(2*ls123[7]+1)
                        ssum1.=ssum1.+conj(state[:,lis1])*transpose(state[:,lis1])
                    end
                    for lisi in 1:a
                        for lisj in 1:a
                            at1[i1nn,1]=at1[i1nn,1]+operator_v[lisi,lisj]*ssum1[lisi,lisj]
                            if lisi==lisj
                                at1[i1nn,2]=at1[i1nn,2]+ssum1[lisi,lisj]
                            end
                            at1[i1nn,3]=at1[i1nn,3]+operator_q1[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,4]=at1[i1nn,4]+operator_q2[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,5]=at1[i1nn,5]+operator_q4[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,6]=at1[i1nn,6]+ov123[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,7]=at1[i1nn,7]+ov135[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,8]=at1[i1nn,8]+ov156[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,9]=at1[i1nn,9]+ov126[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,10]=at1[i1nn,10]+ov345[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,11]=at1[i1nn,11]+ov456[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,12]=at1[i1nn,12]+ov246[lisi,lisj]*ssum1[lisi,lisj]
                            at1[i1nn,13]=at1[i1nn,13]+ov234[lisi,lisj]*ssum1[lisi,lisj]
    #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
    #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                        end
                    end
                    # for mj1 in 1:Int32(2*ls123[7]+1)
                    #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                    #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                    #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                    #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                    #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                    #     at1[i1nn,10]=at1[i1nn,10]+adjoint(state[:,mj1])*ov123*state[:,mj1]
                    #     at1[i1nn,11]=at1[i1nn,11]+adjoint(state[:,mj1])*ov135*state[:,mj1]
                    #     at1[i1nn,12]=at1[i1nn,12]+adjoint(state[:,mj1])*ov156*state[:,mj1]
                    #     at1[i1nn,13]=at1[i1nn,13]+adjoint(state[:,mj1])*ov126*state[:,mj1]
                    #     at1[i1nn,14]=at1[i1nn,14]+adjoint(state[:,mj1])*ov345*state[:,mj1]
                    #     at1[i1nn,15]=at1[i1nn,15]+adjoint(state[:,mj1])*ov456*state[:,mj1]
                    #     at1[i1nn,16]=at1[i1nn,16]+adjoint(state[:,mj1])*ov246*state[:,mj1]
                    #     at1[i1nn,17]=at1[i1nn,17]+adjoint(state[:,mj1])*ov234*state[:,mj1]
                    # end
                    state=nothing
                end
            else
                at1[i1nn,:].=at1[i1nn,:].+im*0.0
            end
            Al=nothing
            a=nothing
        end
    end
return 0
end

function tran123gv7Newtest2(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls11=Vector{Matrix{Float64}}()
    ls21=[]
    seq7oModN(js1,js2,js3,ls11,ls21)
    @everywhere ls1=$ls11
    @everywhere ls2=$ls21
    @sync @distributed for i1nn in 1:size(ls1,1)
            ls123=ls1[Int(ls2[i1nn])][:]
            if ls123[7]<=ls123[1]+ls123[2]+ls123[3]+ls123[1]+ls123[2]+ls123[3]
                g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
                g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
                g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
                Al=Fai7(ls123)
                a=size(Al,1)
                if size(Al)[1]!=0
                        # vol1q=im/4*(v1232(ls123,Al,1,2,3,7,jm61,jm62,B1,B2)+v1232(ls123,Al,1,3,5,7,jm61,jm62,B1,B2)+v1232(ls123,Al,1,5,6,7,jm61,jm62,B1,B2)-v1232(ls123,Al,1,2,6,7,jm61,jm62,B1,B2)+v1232(ls123,Al,3,4,5,7,jm61,jm62,B1,B2)-v1232(ls123,Al,4,5,6,7,jm61,jm62,B1,B2)-v1232(ls123,Al,2,4,6,7,jm61,jm62,B1,B2)-v1232(ls123,Al,2,3,4,7,jm61,jm62,B1,B2))
                        # (e1,ev)=eigen(vol1q)
                        # av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                        # operator_v=av1*(t1)^(3/2)
                        # operator_q1=(vol1q*t1^3)
                        # operator_q2=(vol1q*t1^3)^(2)
                        # operator_q4=(vol1q*t1^3)^(2*2)
                        # operator_q6=(vol1q*t1^3)^(2*3)
                        # operator_q8=(vol1q*t1^3)^(2*4)
                        # operator_q10=(vol1q*t1^3)^(2*5)
                        # operator_q12=(vol1q*t1^3)^(2*6)
                        lls1=Array{Float64,2}(undef,a,6)
                        lls2=Array{Float64,2}(undef,a,6)
                        lls3=Array{Float64,2}(undef,a,6)
                        lls4=Array{Float64,2}(undef,a,6)
                        lls5=Array{Float64,2}(undef,a,6)
                        ij1=Array{Int32}(undef, a, 2)
                        lls1[1:a,1].=ls123[1]
                        lls1[1:a,2].=ls123[2]
                        lls1[1:a,3].=view(Al,1:a,2)
                        lls2[1:a,1].=view(Al,1:a,2)
                        lls2[1:a,2].=ls123[3]
                        lls2[1:a,3].=view(Al,1:a,3)
                        lls3[1:a,1].=view(Al,1:a,4)
                        lls3[1:a,2].=ls123[4]
                        lls3[1:a,3].=view(Al,1:a,3)
                        lls4[1:a,1].=view(Al,1:a,5)
                        lls4[1:a,2].=ls123[5]
                        lls4[1:a,3].=view(Al,1:a,4)
                        lls5[1:a,1].=ls123[6]
                        lls5[1:a,2].=ls123[7]
                        lls5[1:a,3].=view(Al,1:a,5)
                        i2=Array{Float64,1}(undef,a)
                        ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                        ij1=Array{Int32}(undef, a, 2)
                        bi=Array{Float64,1}(undef,a)
                        ci=Array{Float64,1}(undef,a)
                        di=Array{Float64,1}(undef,a)
                        ei=Array{Float64,1}(undef,a)
                        ja1=Array{Float64,1}(undef,a)
                        ja2=Array{Float64,1}(undef,a)
                        ja3=Array{Float64,1}(undef,a)
                        ja=Array{Float64,1}(undef,a)
                        i=Array{Float64,1}(undef,a)
                        j=Array{Float64,1}(undef,a)
                        rs1=Array{Int32}(undef, a, 2)
                        rs11=Array{Int32}(undef, a)
                        rs12=Array{Int32}(undef, a)
                        rs2=zeros(a)
                        sA=sqrt.(2*Al[:,2].+1).*sqrt.(2*Al[:,3].+1).*sqrt.(2*Al[:,4].+1).*sqrt.(2*Al[:,5].+1).*sqrt.(2*ls123[7]+1)
                        state=zeros(a,Int32(2*ls123[7]+1))*im
                        for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)*(2*ls123[6]+1)
                            m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                            m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                            m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                            m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                            m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                            m6=ls123[6]-rem(fld(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1)),(2*ls123[6]+1))
                            s1 = -m1-m2
                            s2 = s1-m3
                            s3 = s2-m4
                            s4 = s3-m5
                            m7 = s4-m6
                            if abs(m7)<=ls123[7]
                                @view(lls1[1:a,4]).=m1
                                @view(lls1[1:a,5]).=m2
                                @view(lls1[1:a,6]).=s1
                                @view(lls2[1:a,4]).=-s1
                                @view(lls2[1:a,5]).=m3
                                @view(lls2[1:a,6]).=s2
                                @view(lls3[1:a,4]).=s3
                                @view(lls3[1:a,5]).=m4  
                                @view(lls3[1:a,6]).=-s2
                                @view(lls4[1:a,4]).=s4
                                @view(lls4[1:a,5]).=m5  
                                @view(lls4[1:a,6]).=-s3
                                @view(lls5[1:a,4]).=m6
                                @view(lls5[1:a,5]).=m7
                                @view(lls5[1:a,6]).=-s4
                                # print(ls123,[m1,m2,m3,m4,m5,m6])
                                com=TransformM6Final2(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                                Inter7New(lls1, lls2, lls3, lls4, lls5, aj1, jm, a, s1, s2, s3, s4, rs1, rs2, rs11, rs12, ij1, ja1, ja2, ja3, ja, i, j, i2, bi, ci, di, ei)
                                inter5n=i2[:].*sA
                                state[:,Int32(ls123[7]+m7+1)].=state[:,Int32(ls123[7]+m7+1)].+inter5n
                            end
                        end
                        # for mj1 in 1:Int32(2*ls123[7]+1)
                        #     at1[i1nn,1]=at1[i1nn,1]+adjoint(state[:,mj1])*operator_v*state[:,mj1]
                        #     at1[i1nn,2]=at1[i1nn,2]+adjoint(state[:,mj1])*state[:,mj1]
                        #     at1[i1nn,3]=at1[i1nn,3]+adjoint(state[:,mj1])*operator_q1*state[:,mj1]
                        #     at1[i1nn,4]=at1[i1nn,4]+adjoint(state[:,mj1])*operator_q2*state[:,mj1]
                        #     at1[i1nn,5]=at1[i1nn,5]+adjoint(state[:,mj1])*operator_q4*state[:,mj1]
                        #     at1[i1nn,6]=at1[i1nn,6]+adjoint(state[:,mj1])*operator_q6*state[:,mj1]
                        #     at1[i1nn,7]=at1[i1nn,7]+adjoint(state[:,mj1])*operator_q8*state[:,mj1]
                        #     at1[i1nn,8]=at1[i1nn,8]+adjoint(state[:,mj1])*operator_q10*state[:,mj1]
                        #     at1[i1nn,9]=at1[i1nn,9]+adjoint(state[:,mj1])*operator_q12*state[:,mj1]
                        # end
                        state=nothing
                else
                    at1[i1nn,:].=at1[i1nn,:].+im*0.0
                end
                Al=nothing
                a=nothing
            end
    end
    return 0
end

# function tran123gv7(js1,js2,js3,js4,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
#     ls1=zeros((js1+1)*(js2+1)*(js3+1),7)
#     ls2=zeros((js1+1)*(js2+1)*(js3+1))
#     seq7oMod(js1,js2,js3,ls1,ls2)
#     @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js1+1)*(js2+1)*(js3+1)
#         for ij5 in 1:(js1+js2+js3+js1+js2+js3+1)
#             ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij5-1)/2]
#             g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
#             g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
#             g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
#             Al=Fai7(ls123)
#             a=size(Al,1)
#             if size(Al)[1]!=0
#                     vol1q=im/4*(v1232(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v1232(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v1232(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v1232(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v1232(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v1232(ls123,Al,2,3,4,6,jm61,jm62,B1,B2))
#                     (e1,ev)=eigen(vol1q)
#                     av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
#                     operator_v=av1*(t1)^(3/2)
#                     operator_q1=(vol1q*t1^3)
#                     operator_q2=(vol1q*t1^3)^(2)
#                     operator_q4=(vol1q*t1^3)^(2*2)
#                     operator_q6=(vol1q*t1^3)^(2*3)
#                     operator_q8=(vol1q*t1^3)^(2*4)
#                     operator_q10=(vol1q*t1^3)^(2*5)
#                     operator_q12=(vol1q*t1^3)^(2*6)
#                     ij1=Array{Int32}(undef, a, 2)
#                     lls1=Array{Float64,2}(undef,a,6)
#                     lls2=Array{Float64,2}(undef,a,6)
#                     lls3=Array{Float64,2}(undef,a,6)
#                     lls1[1:a,1].=ls1232[1]
#                     lls1[1:a,2].=ls1232[2]
#                     lls1[1:a,3].=view(Al2,1:a,2)
#                     lls2[1:a,1].=view(Al2,1:a,2)
#                     lls2[1:a,2].=ls1232[3]
#                     lls2[1:a,3].=view(Al2,1:a,3)
#                     lls3[1:a,1].=view(Al2,1:a,3)
#                     lls3[1:a,2].=ls1232[3]
#                     lls3[1:a,3].=view(Al2,1:a,4)
#                     lls4[1:a,1].=view(Al2,1:a,5)
#                     lls4[1:a,2].=ls1232[3]
#                     lls4[1:a,3].=view(Al2,1:a,4)
#                     lls5[1:a,1].=ls1232[7]
#                     lls5[1:a,2].=ls1232[6]
#                     lls5[1:a,3].=view(Al2,1:a,5)
#                     i2=zeros(a)
#                     state=zeros(size(a1)[1],a,Int(2*ls1232[4]+1))*im
#                     kc1=zeros(a,aa22)*im
#                     ssum1=zeros(a,a)*im
#                     for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1))
#                         mtag=kronmc3(ja,(ls1232[1:3]*2).+1)
#                         s1 = -mtag[1]-mtag[2]
#                         m4 = s1-mtag[3]
#                         if abs(m4)<=ls1232[4]
#                             @view(lls1[1:aa22,4]).=mtag[1]
#                             @view(lls1[1:aa22,5]).=mtag[2]
#                             @view(lls1[1:aa22,6]).=s1
#                             @view(lls2[1:aa22,4]).=-s1
#                             @view(lls2[1:aa22,5]).=mtag[3]  
#                             @view(lls2[1:aa22,6]).=m4
#                             @view(lls3[1:aa22,4]).=-s1
#                             @view(lls3[1:aa22,5]).=mtag[3]  
#                             @view(lls3[1:aa22,6]).=m4
#                             @view(lls4[1:aa22,4]).=-s1
#                             @view(lls4[1:aa22,5]).=mtag[3]  
#                             @view(lls4[1:aa22,6]).=m4
#                             @view(lls5[1:aa22,4]).=-s1
#                             @view(lls5[1:aa22,5]).=mtag[3]  
#                             @view(lls5[1:aa22,6]).=m4
#                             com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
#                             Inter7_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
#                             inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*Al2[:,3].+1).*sqrt.(2*Al2[:,4].+1).*sqrt.(2*Al2[:,5].+1).*sqrt.(2*ls123[7]+1)
#                             for jml in  1:size(a1)[1]
#                                 kron!(kc1,a1[jml,:,ja],transpose(inter5n))
#                                 state[jml,:,:,Int(ls1232[4]+m4+1)].=state[jml,:,:,Int(ls1232[4]+m4+1)].+kc1
#                             end
#                         end
#             end
#                     for lis1 in 1:Int(2*ls123[4]+1)*aa22*Int(2*ls1232[4]+1)
#                         is=Int(rem(lis1-1,aa22)+1)
#                         jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[4]+1))+1)
#                         jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[4]+1)),Int(2*ls1232[4]+1))+1)
#                         ssum1.=ssum1.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
#                     end
#                     for lisi in 1:a
#                         for lisj in 1:a
#                             at1[i1,1]=at1[i1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
#                             if lisi==lisj
#                                 at1[i1,2]=at1[i1,2]+ssum1[lisi,lisj]
#                             end
#                             at1[i1,3]=at1[i1,3]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
#                             for ii in 4:nn
#                                 at1[i1,ii]=at1[i1,ii]+((vol1q*t1^3/4*im)^(2*(ii-3)))[lisi,lisj]*ssum1[lisi,lisj]
#                             end
#         #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
#         #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
#         #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
#         #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
#                         end
#                     end
#                     kc1=nothing
#                     state=nothing
#                 end
#             end
#     else
#         at1[i1,:].=at1[i1,:].+im*0.0
#     end
#     Al=nothing
#     a=nothing
# end
#     end
#     return 0
# end

function tranVtest_s666overlap(js1,js2,js3,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm61,jm62,B1,B2)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls123=zeros(6)
    ls123[1]=js1
    ls123[2]=js2
    ls123[3]=js3
    ls123[4]=js1
    ls123[5]=js2
    ls123[6]=js3
        g1c=exp(-t1*(lb1(ls123[1])))*conj(glc(ls123[1],gf1))
        g2c=exp(-t1*(lb1(ls123[2])))*conj(glc(ls123[2],gf2))
        g3c=exp(-t1*(lb1(ls123[3])))*conj(glc(ls123[3],gf3))
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai6(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                vol1q=im/4*(v123(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v123(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v123(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,3,4,6,jm61,jm62,B1,B2))
                (e1,ev)=eigen(vol1q)
                maxe1=maximum(real.(e1))
                maxe1loc=findall(x->x==maxe1, real.(e1))[1]
                resoverlap=zeros(size(e1,1)+1)
                av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                lls3=Array{Float64,2}(undef,a,6)
                lls4=Array{Float64,2}(undef,a,6)
                # print(size(lls1))
                #Initialzation of ls
                lls1[:,1].=ls123[1]
                lls1[:,2].=ls123[2]
                lls1[:,3].=view(Al,:,2)
                lls2[:,1].=view(Al,:,2)
                lls2[:,2].=ls123[3]
                lls2[:,3].=view(Al,:,3)
                lls3[:,1].=view(Al,:,4)
                lls3[:,2].=ls123[4]
                lls3[:,3].=view(Al,:,3)
                lls4[:,1].=ls123[5]
                lls4[:,2].=ls123[6]
                lls4[:,3].=view(Al,1:a,4)
                i2=Array{Float64,1}(undef,a)
                ic1=sqrt(2*ls123[1]+1)*sqrt(2*ls123[2]+1)*sqrt(2*ls123[3]+1)
                bi=Array{Float64,1}(undef,a)
                ci=Array{Float64,1}(undef,a)
                di=Array{Float64,1}(undef,a)
                ja1=Array{Float64,1}(undef,a)
                ja2=Array{Float64,1}(undef,a)
                ja3=Array{Float64,1}(undef,a)
                ja=Array{Float64,1}(undef,a)
                i=Array{Float64,1}(undef,a)
                j=Array{Float64,1}(undef,a)
                rs1=Array{Int32}(undef, a, 2)
                rs11=Array{Int32}(undef, a)
                rs12=Array{Int32}(undef, a)
                rs2=zeros(a)
                sA=sqrt.((Al[:,2].*2).+1).*sqrt.((Al[:,3].*2).+1).*sqrt.((Al[:,4].*2).+1)
                state=zeros(a)*im
                for ma in 1:(2*ls123[1]+1)*(2*ls123[2]+1)*(2*ls123[3]+1)*(2*ls123[4]+1)*(2*ls123[5]+1)
                    #Block One
                    m1=ls123[1]-rem(ma-1,(2*ls123[1]+1))
                    m2=ls123[2]-rem(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1))
                    m3=ls123[3]-rem(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1))
                    m4=ls123[4]-rem(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1))
                    m5=ls123[5]-rem(fld(fld(fld(fld(ma-1,(2*ls123[1]+1)),(2*ls123[2]+1)),(2*ls123[3]+1)),(2*ls123[4]+1)),(2*ls123[5]+1))
                    s1 = -m1-m2
                    s2 = s1-m3
                    s3 = s2-m4
                    m6 = s3-m5
                    if abs(m6)<=ls123[6]
                    #Block Two    
                        lls1[1:a,4].=m1
                        lls1[1:a,5].=m2
                        lls1[1:a,6].=s1
                        lls2[1:a,4].=-s1
                        lls2[1:a,5].=m3  
                        lls2[1:a,6].=s2
                        lls3[1:a,4].=s3
                        lls3[1:a,5].=m4  
                        lls3[1:a,6].=-s2
                        lls4[1:a,4].=m5
                        lls4[1:a,5].=m6
                        lls4[1:a,6].=-s3
        #                 if j6+m6+1==1.5
        #                     print([j1 j2 j3 j4 j5 j6; m1 m2 m3 m4 m5 m6])
        #                 end
                    #Block Three
                        Inter6_mod2(lls1,lls2,lls3,lls4,aj1,jm,a,s1,s2,s3,rs1,rs2,rs11,rs12,ij1,ja1,ja2,ja3,ja,i,j, i2, bi, ci, di) 
                        com=TransformM6Final(ls123,g1c,g2c,g3c,[m1,m2,m3,m4,m5,m6])
                        state[:].=state[:].+(sA.*i2*ic1*com)
                    end
                end
                for is in 1:a
                    resoverlap[1:size(e1,1)].=resoverlap[1:size(e1,1)].+abs.(inv(ev)[:,:]*state[:,is]).^2
                    resoverlap[size(e1,1)+1]=resoverlap[size(e1,1)+1]+adjoint(state[:,is])*state[:,is]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
        print(resoverlap)
        print(e1)
        print(ev)
    return (resoverlap,e1,ev)
end

function tran123testp(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,hf1,hf2,hf3,hf4,jm,aj1,center,dist,at1,nn)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        a=size(Al,1)
        ij1=Array{Int32}(undef, a, 2)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                a2=Trans4ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                lls1=Array{Float64,2}(undef,a,6)
                lls2=Array{Float64,2}(undef,a,6)
                #Initialzation of ls
                lls1[1:a,1].=ls123[1]
                lls1[1:a,2].=ls123[2]
                lls1[1:a,3].=view(Al,1:a,2)
                lls2[1:a,1].=view(Al,1:a,2)
                lls2[1:a,2].=ls123[3]
                lls2[1:a,3].=ls123[4]
                i2=Array{Float64,1}(undef,a)
                state=zeros(a,a)*im
                for ja in 1:size(a1)[2]
                    mtag=kronmc1(ja,(ls123*2).+1)
#                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                    if mtag[4]==-mtag[1]-mtag[2]-mtag[3]
                        s1 = -mtag[1]-mtag[2]
                        @view(lls1[1:a,4]).=mtag[1]
                        @view(lls1[1:a,5]).=mtag[2]
                        @view(lls1[1:a,6]).=s1
                        @view(lls2[1:a,4]).=-s1
                        @view(lls2[1:a,5]).=mtag[3]  
                        @view(lls2[1:a,6]).=mtag[4]
                        Inter4_Compact(lls1, lls2, aj1, jm, a, s1, ij1, i2)
                        #println(a[:,ja])
                        for iia2 in 1:a
                            state[:,iia2].=state[:,iia2].+(a1[:,ja]*i2[iia2])*sqrt(2*Al[iia2,2]+1)
                        end
                    end
                end
                for is in 1:a
                    at1[i1,1]=at1[i1,1]+adjoint(state[:,is])*av1*(t1)^(3/2)*state[:,is]
                    at1[i1,2]=at1[i1,2]+adjoint(state[:,is])*state[:,is]
                    for ii in 3:nn
                        at1[i1,ii]=at1[i1,ii]+adjoint(state[:,is])*((vol1q*t1^3/8*im)^(ii-2))*state[:,is]
                    end
                    at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[:,is])*jl1*state[:,is]*t1^2
                    at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[:,is])*jl2*state[:,is]*t1^2
                    at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[:,is])*jl3*state[:,is]*t1^2
                    at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[:,is])*jl4*state[:,is]*t1^2
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,:].=im*0.0
            end
        else
            at1[i1,:].=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123gvm(js1,js2,js3,js4,js5,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
#     print(1)
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        for ij5 in 1:(js1+js2+js3+js4+1)
            for ij52 in 1:(js1+js2+js3+js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
                ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij5-1)/2]
                ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij52-1)/2]
                #ls123=ls1[i1,:]
                #@show ls123
                #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
                #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
                Al=Fai5(ls123)
                Al2=Fai5(ls1232)
                a=size(Al,1)
                aa22=size(Al2,1)
                ij1=Array{Int32}(undef, aa22, 2)
                if size(Al)[1]!=0
                    if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                        a1=Trans5ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
    #                     print(a1)
                        vol1q=v123(ls123,Al,1,2,3,5,jm1,jm2,A1,B1)-v123(ls123,Al,2,3,4,5,jm1,jm2,A1,B1)-v123(ls123,Al,1,2,4,5,jm1,jm2,A1,B1)+v123(ls123,Al,1,3,4,5,jm1,jm2,A1,B1)
    #                     print(vol1q)
                        if ls123[1]!=0
                            jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                        else
                            jl1=0
                        end
                        if ls123[2]!=0
                            jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                        else
                            jl2=0
                        end
                        if ls123[3]!=0
                            jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                        else
                            jl3=0
                        end
                        if ls123[4]!=0
                            jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                        else
                            jl4=0
                        end
                        (e1,ev)=eigen(vol1q)
                        av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                        lls1=Array{Float64,2}(undef,aa22,6)
                        lls2=Array{Float64,2}(undef,aa22,6)
                        lls3=Array{Float64,2}(undef,aa22,6)
                        #Initialzation of ls
                        lls1[1:aa22,1].=ls1232[1]
                        lls1[1:aa22,2].=ls1232[2]
                        lls1[1:aa22,3].=view(Al2,1:aa22,2)
                        lls2[1:aa22,1].=view(Al2,1:aa22,2)
                        lls2[1:aa22,2].=ls1232[3]
                        lls2[1:aa22,3].=view(Al2,1:aa22,3)
                        lls3[1:aa22,1].=view(Al2,1:aa22,3)
                        lls3[1:aa22,2].=ls1232[4]
                        lls3[1:aa22,3].=ls1232[5]
                        i2=zeros(aa22)
                        state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[5]+1))*im
                        for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1)*(2*ls1232[4]+1))
                            mtag=kronmc1(ja,(ls1232[1:4]*2).+1)
#                             println("m=",mtag)
        #                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                                s1 = -mtag[1]-mtag[2]
                                s2 = s1-mtag[3]
                                m5 = s2-mtag[4]
                                if abs(m5)<=ls1232[5]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=s2
                                    @view(lls3[1:aa22,4]).=-s2
                                    @view(lls3[1:aa22,5]).=mtag[4]  
                                    @view(lls3[1:aa22,6]).=m5
                                    Inter5_Compact(lls1, lls2, lls3, aj1, jm, aa22, s1, s2, ij1, i2)
                                    #println(a[:,ja])
                                        for iia2 in 1:aa22
                                            state[:,iia2].=state[:,iia2].+(a1[jml,:,ja]*i2[iia2])*sqrt(2*Al2[iia2,2]+1)*sqrt(2*Al2[iia2,3]+1)
                                        end
                                end
                        end
                            for is in 1:aa22
                                for jml1 in  1:size(a1)[1]
                                    for jml2 in  1:Int(2*ls1232[5]+1)
                                        at1[i1,1]=at1[i1,1]+adjoint(state[jml1,:,is,jml2])*av1*(t1)^(3/2)*state[jml1,:,is,jml2]
                                        at1[i1,2]=at1[i1,2]+adjoint(state[jml1,:,is,jml2])*state[jml1,:,is,jml2]
                                        for ii in 3:nn
                                            at1[i1,ii]=at1[i1,ii]+adjoint(state[jml1,:,is,jml2])*((vol1q*t1^3/8*im)^(2*(ii-2)))*state[jml1,:,is,jml2]
                                        end
                                        at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                    end
                                end
                            end
    #                     ja=1
    #                     at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
    #                     at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
    #                     println(ls123," ",at1[i1,2])
                    else
                        at1[i1,:].=at1[i1,:].+im*0.0
                    end
                else
                    at1[i1,:].=at1[i1,:].+im*0.0
                end
                Al=nothing
                a=nothing
#                 println("j2=",ls1232)
            end
#             println("j1=",ij5)
        end
    end
    return 0
end

function tran123gv3(js1,js2,js3,js4,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls1=zeros((js1+1)*(js2+1)*(js3+1),3)
    ls2=zeros((js1+1)*(js2+1)*(js3+1))
    seq3o(js1,js2,js3,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)
        for ij5 in 1:(js1+js2+js3+1)
            ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij5-1)/2]
            Al=Fai4(ls123)
            a=size(Al,1)
            if size(Al)[1]!=0
                if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                    a1=Trans3ott(ls123,t1,gf1,gf2,gf3,jm,aj1,Al,size(Al,1))
                    vol1q=v123(ls123,Al,1,2,3,4,jm1,jm2,A1,B1)
                    (e1,ev)=eigen(vol1q)
                    av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                    for ij52 in 1:(js1+js2+js3+js4+1)
                        ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij52-1)/2]
                        Al2=Fai4(ls1232)
                        aa22=size(Al2,1)
                        if aa22!=0
                            ij1=Array{Int32}(undef, aa22, 2)
                            lls1=Array{Float64,2}(undef,aa22,6)
                            lls2=Array{Float64,2}(undef,aa22,6)
                            lls3=Array{Float64,2}(undef,aa22,6)
                            lls1[1:aa22,1].=ls1232[1]
                            lls1[1:aa22,2].=ls1232[2]
                            lls1[1:aa22,3].=view(Al2,1:aa22,2)
                            lls2[1:aa22,1].=view(Al2,1:aa22,2)
                            lls2[1:aa22,2].=ls1232[3]
                            lls2[1:aa22,3].=ls1232[4]
                            i2=zeros(aa22)
                            state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            kc1=zeros(a,aa22)*im
                            ssum1=zeros(a,a)*im
                            for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1))
                                mtag=kronmc3(ja,(ls1232[1:3]*2).+1)
                                s1 = -mtag[1]-mtag[2]
                                m4 = s1-mtag[3]
                                if abs(m4)<=ls1232[4]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=m4
                                    Inter4_Compact(lls1, lls2, aj1, jm, aa22, s1, ij1, i2)
                                    inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*ls1232[4]+1)
                                    for jml in  1:size(a1)[1]
                                        kron!(kc1,a1[jml,:,ja],transpose(inter5n))
                                        state[jml,:,:,Int(ls1232[4]+m4+1)].=state[jml,:,:,Int(ls1232[4]+m4+1)].+kc1
                                    end
                                end
                            end
                            for lis1 in 1:Int(2*ls123[4]+1)*aa22*Int(2*ls1232[4]+1)
                                is=Int(rem(lis1-1,aa22)+1)
                                jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[4]+1))+1)
                                jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[4]+1)),Int(2*ls1232[4]+1))+1)
                                ssum1.=ssum1.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                            end
                            for lisi in 1:a
                                for lisj in 1:a
                                    at1[i1,1]=at1[i1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
                                    if lisi==lisj
                                        at1[i1,2]=at1[i1,2]+ssum1[lisi,lisj]
                                    end
                                    at1[i1,3]=at1[i1,3]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
                                    for ii in 4:nn
                                        at1[i1,ii]=at1[i1,ii]+((vol1q*t1^3/4*im)^(2*(ii-3)))[lisi,lisj]*ssum1[lisi,lisj]
                                    end
        #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                end
                            end
                            kc1=nothing
                            state=nothing
                        end
                    end
                else
                    at1[i1,:].=at1[i1,:].+im*0.0
                end
            else
                at1[i1,:].=at1[i1,:].+im*0.0
            end
            Al=nothing
            a=nothing
        end
    end
    return 0
end

function tran123gv3j(js1,js2,js3,js4,t1,gf1,gf2,gf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls1=zeros((js1+1)*(js2+1)*(js3+1),3)
    ls2=zeros((js1+1)*(js2+1)*(js3+1))
    seq3o(js1,js2,js3,ls1,ls2)
        for ij5 in 1:(js1+js2+js3+1)
            ls123=[js1/2 js2/2 js3/2 (ij5-1)/2]
            Al=Fai4(ls123)
            a=size(Al,1)
            if size(Al)[1]!=0
                if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                    a1=Trans3ott(ls123,t1,gf1,gf2,gf3,jm,aj1,Al,size(Al,1))
                    vol1q=qa(ls123)
                    (e1,ev)=eigen(vol1q)
                    av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                    for ij52 in 1:(js1+js2+js3+js4+1)
                        ls1232=[js1/2 js2/2 js3/2 (ij52-1)/2]
                        Al2=Fai4(ls1232)
                        aa22=size(Al2,1)
                        if aa22!=0
                            ij1=Array{Int32}(undef, aa22, 2)
                            lls1=Array{Float64,2}(undef,aa22,6)
                            lls2=Array{Float64,2}(undef,aa22,6)
                            lls3=Array{Float64,2}(undef,aa22,6)
                            lls1[1:aa22,1].=ls1232[1]
                            lls1[1:aa22,2].=ls1232[2]
                            lls1[1:aa22,3].=view(Al2,1:aa22,2)
                            lls2[1:aa22,1].=view(Al2,1:aa22,2)
                            lls2[1:aa22,2].=ls1232[3]
                            lls2[1:aa22,3].=ls1232[4]
                            i2=zeros(aa22)
                            state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            kc1=zeros(a,aa22)*im
                            ssum1=zeros(a,a)*im
                            for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1))
                                mtag=kronmc3(ja,(ls1232[1:3]*2).+1)
                                s1 = -mtag[1]-mtag[2]
                                m4 = s1-mtag[3]
                                if abs(m4)<=ls1232[4]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=m4
                                    Inter4_Compact(lls1, lls2, aj1, jm, aa22, s1, ij1, i2)
                                    inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*ls1232[4]+1)
                                    for jml in  1:size(a1)[1]
                                        kron!(kc1,a1[jml,:,ja],transpose(inter5n))
                                        state[jml,:,:,Int(ls1232[4]+m4+1)].=state[jml,:,:,Int(ls1232[4]+m4+1)].+kc1
                                    end
                                end
                            end
                            for lis1 in 1:Int(2*ls123[4]+1)*aa22*Int(2*ls1232[4]+1)
                                is=Int(rem(lis1-1,aa22)+1)
                                jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[4]+1))+1)
                                jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[4]+1)),Int(2*ls1232[4]+1))+1)
                                ssum1.=ssum1.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                            end
                            for lisi in 1:a
                                for lisj in 1:a
                                    at1[js1,1]=at1[js1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
                                    if lisi==lisj
                                        at1[js1,2]=at1[js1,2]+ssum1[lisi,lisj]
                                    end
                                    at1[js1,3]=at1[js1,3]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
                                    for ii in 4:nn
                                        at1[js1,ii]=at1[js1,ii]+((vol1q*t1^3/4*im)^(2*(ii-3)))[lisi,lisj]*ssum1[lisi,lisj]
                                    end
        #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                end
                            end
                            kc1=nothing
                            state=nothing
                        end
                    end
                else
                    at1[js1,:].=at1[js1,:].+im*0.0
                end
            else
                at1[js1,:].=at1[js1,:].+im*0.0
            end
            Al=nothing
            a=nothing
        end
    return 0
end

function tran123gv3p(js1,js2,js3,js4,t1,gf1,gf2,gf3,hf1,hf2,hf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls1=zeros((js1+1)*(js2+1)*(js3+1),3)
    ls2=zeros((js1+1)*(js2+1)*(js3+1))
    seq3o(js1,js2,js3,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)
        for ij5 in 1:(js1+js2+js3+1)
            ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij5-1)/2]
            Al=Fai4(ls123)
            a=size(Al,1)
            if size(Al)[1]!=0
                if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                    a1=Trans3ott(ls123,t1,gf1,gf2,gf3,jm,aj1,Al,size(Al,1))
                    a2=Trans3ott(ls123,t1,hf1,hf2,hf3,jm,aj1,Al,size(Al,1))
                    vol1q=qa(ls123)*im/4
                    (e1,ev)=eigen(vol1q)
                    av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                    for ij52 in 1:(js1+js2+js3+js4+1)
                        ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij52-1)/2]
                        Al2=Fai4(ls1232)
                        aa22=size(Al2,1)
                        if aa22!=0
                            ij1=Array{Int32}(undef, aa22, 2)
                            lls1=Array{Float64,2}(undef,aa22,6)
                            lls2=Array{Float64,2}(undef,aa22,6)
                            lls3=Array{Float64,2}(undef,aa22,6)
                            lls1[1:aa22,1].=ls1232[1]
                            lls1[1:aa22,2].=ls1232[2]
                            lls1[1:aa22,3].=view(Al2,1:aa22,2)
                            lls2[1:aa22,1].=view(Al2,1:aa22,2)
                            lls2[1:aa22,2].=ls1232[3]
                            lls2[1:aa22,3].=ls1232[4]
                            i2=zeros(aa22)
                            state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            statep=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            kc1=zeros(a,aa22)*im
                            kc2=zeros(a,aa22)*im
                            ssum1=zeros(a,a)*im
                            ssum2=zeros(a,a)*im
                            ssum3=zeros(a,a)*im
                            asum1=0
                            asum2=0
                            for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1))
                                mtag=kronmc3(ja,(ls1232[1:3]*2).+1)
                                s1 = -mtag[1]-mtag[2]
                                m4 = s1-mtag[3]
                                if abs(m4)<=ls1232[4]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=m4
                                    Inter4_Compact(lls1, lls2, aj1, jm, aa22, s1, ij1, i2)
                                    inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*ls1232[4]+1)
                                    for jml in  1:size(a1)[1]
                                        kron!(kc1,a1[jml,:,ja],transpose(inter5n))
                                        kron!(kc2,a2[jml,:,ja],transpose(inter5n))
                                        state[jml,:,:,Int(ls1232[4]+m4+1)].=state[jml,:,:,Int(ls1232[4]+m4+1)].+kc1
                                        statep[jml,:,:,Int(ls1232[4]+m4+1)].=statep[jml,:,:,Int(ls1232[4]+m4+1)].+kc2
                                    end
                                end
                            end
                            for lis1 in 1:Int(2*ls123[4]+1)*aa22*Int(2*ls1232[4]+1)
                                is=Int(rem(lis1-1,aa22)+1)
                                jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[4]+1))+1)
                                jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[4]+1)),Int(2*ls1232[4]+1))+1)
                                ssum1.=ssum1.+conj(statep[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                                ssum2.=ssum2.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                                ssum3.=ssum3.+conj(statep[jml1,:,is,jml2])*transpose(statep[jml1,:,is,jml2])
                            end
                            for lisi in 1:a
                                for lisj in 1:a
                                    at1[i1,1]=at1[i1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
                                    if lisi==lisj
                                        at1[i1,2]=at1[i1,2]+ssum2[lisi,lisj]
                                        at1[i1,3]=at1[i1,3]+ssum3[lisi,lisj]
                                    end
                                    at1[i1,4]=at1[i1,4]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
                                    for ii in 5:nn
                                        at1[i1,ii]=at1[i1,ii]+((vol1q*t1^3/4*im)^(2*(ii-4)))[lisi,lisj]*ssum1[lisi,lisj]
                                    end
        #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                end
                            end
                            kc1=nothing
                            kc2=nothing
                            state=nothing
                            statep=nothing
                        end
                    end
                else
                    at1[i1,:].=at1[i1,:].+im*0.0
                end
            else
                at1[i1,:].=at1[i1,:].+im*0.0
            end
            Al=nothing
            a=nothing
        end
    end
    return 0
end

function tran123gv3ptest(js1,js2,js3,js4,t1,gf1,gf2,gf3,hf1,hf2,hf3,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls1=zeros((js1+1)*(js2+1)*(js3+1),3)
    ls2=zeros((js1+1)*(js2+1)*(js3+1))
    seq3o(js1,js2,js3,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)
        for ij5 in 1:(js1+js2+js3+1)
            ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij5-1)/2]
            Al=Fai4(ls123)
            a=size(Al,1)
            if size(Al)[1]!=0
                if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                    a1=Trans3ott(ls123,t1,gf1,gf2,gf3,jm,aj1,Al,size(Al,1))
                    a2=Trans3ott(ls123,t1,hf1,hf2,hf3,jm,aj1,Al,size(Al,1))
                    vol1q=qa(ls123)*im/4
                    (e1,ev)=eigen(vol1q)
                    av1=ev*diagm((sqrt.(abs.(e1))))*inv(ev)
                    for ij52 in 1:(js1+js2+js3+js4+1)
                        ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] (ij52-1)/2]
                        Al2=Fai4(ls1232)
                        aa22=size(Al2,1)
                        if aa22!=0
                            ij1=Array{Int32}(undef, aa22, 2)
                            lls1=Array{Float64,2}(undef,aa22,6)
                            lls2=Array{Float64,2}(undef,aa22,6)
                            lls3=Array{Float64,2}(undef,aa22,6)
                            lls1[1:aa22,1].=ls1232[1]
                            lls1[1:aa22,2].=ls1232[2]
                            lls1[1:aa22,3].=view(Al2,1:aa22,2)
                            lls2[1:aa22,1].=view(Al2,1:aa22,2)
                            lls2[1:aa22,2].=ls1232[3]
                            lls2[1:aa22,3].=ls1232[4]
                            i2=zeros(aa22)
                            state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            statep=zeros(size(a1)[1],a,aa22,Int(2*ls1232[4]+1))*im
                            kc1=zeros(a,aa22)*im
                            kc2=zeros(a,aa22)*im
                            ssum1=zeros(a,a)*im
                            ssum2=zeros(a,a)*im
                            ssum3=zeros(a,a)*im
                            asum1=0
                            asum2=0
                            for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1))
                                mtag=kronmc3(ja,(ls1232[1:3]*2).+1)
                                s1 = -mtag[1]-mtag[2]
                                m4 = s1-mtag[3]
                                if abs(m4)<=ls1232[4]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=m4
                                    Inter4_Compact(lls1, lls2, aj1, jm, aa22, s1, ij1, i2)
                                    inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*ls1232[4]+1)
                                    for jml in  1:size(a1)[1]
                                        kron!(kc1,a1[jml,:,ja],transpose(inter5n))
                                        kron!(kc2,a2[jml,:,ja],transpose(inter5n))
                                        state[jml,:,:,Int(ls1232[4]+m4+1)].=state[jml,:,:,Int(ls1232[4]+m4+1)].+kc1
                                        statep[jml,:,:,Int(ls1232[4]+m4+1)].=statep[jml,:,:,Int(ls1232[4]+m4+1)].+kc2
                                    end
                                end
                            end
                            for lis1 in 1:Int(2*ls123[4]+1)*aa22*Int(2*ls1232[4]+1)
                                is=Int(rem(lis1-1,aa22)+1)
                                jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[4]+1))+1)
                                jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[4]+1)),Int(2*ls1232[4]+1))+1)
                                ssum1.=ssum1.+conj(statep[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                                ssum2.=ssum2.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                                ssum3.=ssum3.+conj(statep[jml1,:,is,jml2])*transpose(statep[jml1,:,is,jml2])
                            end
                            for lisi in 1:a
                                for lisj in 1:a
                                    at1[i1,1]=at1[i1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
                                    if lisi==lisj
                                        at1[i1,2]=at1[i1,2]+ssum2[lisi,lisj]
                                        at1[i1,3]=at1[i1,3]+ssum3[lisi,lisj]
                                    end
                                    at1[i1,4]=at1[i1,4]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
                                    for ii in 5:nn
                                        at1[i1,ii]=at1[i1,ii]+((vol1q*t1^3/4*im)^(2*(ii-4)))[lisi,lisj]*ssum1[lisi,lisj]
                                    end
        #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                end
                            end
                            kc1=nothing
                            kc2=nothing
                            state=nothing
                            statep=nothing
                        end
                    end
                else
                    at1[i1,:].=at1[i1,:].+im*0.0
                end
            else
                at1[i1,:].=at1[i1,:].+im*0.0
            end
            Al=nothing
            a=nothing
        end
    end
    return 0
end

function tran123gv(js1,js2,js3,js4,js5,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for ii in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)*(js1+js2+js3+js4+1)
        ij5=Int(rem(ii-1,(js1+js2+js3+js4+1))+1)
        i1=Int(rem(fld(ii-1,(js1+js2+js3+js4+1)),(js1+1)*(js2+1)*(js3+1)*(js4+1))+1)
        ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij5-1)/2]
        Al=Fai5(ls123)
        a=size(Al,1)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans5ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                vol1q=v123(ls123,Al,1,2,3,5,jm1,jm2,A1,B1)
                vol2q=v123(ls123,Al,2,3,4,5,jm1,jm2,A1,B1)
                vol3q=v123(ls123,Al,1,2,4,5,jm1,jm2,A1,B1)
                vol4q=v123(ls123,Al,1,3,4,5,jm1,jm2,A1,B1)
                if ls123[1]!=0
                    jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                else
                    jl1=0
                end
                if ls123[2]!=0
                    jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                else
                    jl2=0
                end
                if ls123[3]!=0
                    jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                else
                    jl3=0
                end
                if ls123[4]!=0
                    jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                else
                    jl4=0
                end
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ij52 in 1:(js1+js2+js3+js4+1)
                    ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij52-1)/2]
                    Al2=Fai5(ls1232)
                    aa22=size(Al2,1)
                    if aa22!=0
                        ij1=Array{Int32}(undef, aa22, 2)
                        lls1=Array{Float64,2}(undef,aa22,6)
                        lls2=Array{Float64,2}(undef,aa22,6)
                        lls3=Array{Float64,2}(undef,aa22,6)
                        lls1[1:aa22,1].=ls1232[1]
                        lls1[1:aa22,2].=ls1232[2]
                        lls1[1:aa22,3].=view(Al2,1:aa22,2)
                        lls2[1:aa22,1].=view(Al2,1:aa22,2)
                        lls2[1:aa22,2].=ls1232[3]
                        lls2[1:aa22,3].=view(Al2,1:aa22,3)
                        lls3[1:aa22,1].=view(Al2,1:aa22,3)
                        lls3[1:aa22,2].=ls1232[4]
                        lls3[1:aa22,3].=ls1232[5]
                        i2=zeros(aa22)
                        state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[5]+1))*im
                        kc1=zeros(a,aa22)*im
                        ssum1=zeros(a,a)*im
                        for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1)*(2*ls1232[4]+1))
                            mtag=kronmc1(ja,(ls1232[1:4]*2).+1)
                            s1 = -mtag[1]-mtag[2]
                            s2 = s1-mtag[3]
                            m5 = s2-mtag[4]
                            if abs(m5)<=ls1232[5]
                                @view(lls1[1:aa22,4]).=mtag[1]
                                @view(lls1[1:aa22,5]).=mtag[2]
                                @view(lls1[1:aa22,6]).=s1
                                @view(lls2[1:aa22,4]).=-s1
                                @view(lls2[1:aa22,5]).=mtag[3]  
                                @view(lls2[1:aa22,6]).=s2
                                @view(lls3[1:aa22,4]).=-s2
                                @view(lls3[1:aa22,5]).=mtag[4]  
                                @view(lls3[1:aa22,6]).=m5
                                Inter5_Compact(lls1, lls2, lls3, aj1, jm, aa22, s1, s2, ij1, i2)
                                inter5n=i2[:].*sqrt.(2*Al2[:,2].+1).*sqrt.(2*Al2[:,3].+1).*sqrt.(2*ls1232[5]+1)
                                for jml in  1:size(a1)[1]
                                    kron!(kc1,a1[jml,:,ja],transpose(inter5n))
                                    state[jml,:,:,Int(ls1232[5]+m5+1)].=state[jml,:,:,Int(ls1232[5]+m5+1)].+kc1
                                end
                            end
                        end
                            for lis1 in 1:Int(2*ls123[5]+1)*aa22*Int(2*ls1232[5]+1)
                                is=Int(rem(lis1-1,aa22)+1)
                                jml1=Int(rem(fld(lis1-1,aa22),Int(2*ls123[5]+1))+1)
                                jml2=Int(rem(fld(fld(lis1-1,aa22),Int(2*ls123[5]+1)),Int(2*ls1232[5]+1))+1)
                                ssum1.=ssum1.+conj(state[jml1,:,is,jml2])*transpose(state[jml1,:,is,jml2])
                            end
                            for lisi in 1:a
                                for lisj in 1:a
                                    at1[i1,1]=at1[i1,1]+av1[lisi,lisj]*(t1)^(3/2)*ssum1[lisi,lisj]
                                    if lisi==lisj
                                        at1[i1,2]=at1[i1,2]+ssum1[lisi,lisj]
                                    end
                                    at1[i1,3]=at1[i1,3]+(vol1q*t1^3/4*im)[lisi,lisj]*ssum1[lisi,lisj]
                                    for ii in 4:nn
                                        at1[i1,ii]=at1[i1,ii]+((vol1q*t1^3/4*im)^(2*(ii-2)))[lisi,lisj]*ssum1[lisi,lisj]
                                    end
        #                             at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
        #                             at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                end
                            end
                        kc1=nothing
                        state=nothing
                    end
                end
            else
                at1[i1,:].=at1[i1,:].+im*0.0
            end
        else
            at1[i1,:].=at1[i1,:].+im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran123gv111(js1,js2,js3,js4,js5,t1,gf1,gf2,gf3,gf4,jm,aj1,center,dist,at1,nn,jm1,jm2,A1,B1)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
#     print(1)
    seq4o(js1,js2,js3,js4,ls1,ls2)
    for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        for ij5 in 1:(js1+js2+js3+js4+1)
            for ij52 in 1:(js1+js2+js3+js4+1)
#     for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
                ls123=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij5-1)/2]
                ls1232=[ls1[Int(ls2[i1]),1] ls1[Int(ls2[i1]),2] ls1[Int(ls2[i1]),3] ls1[Int(ls2[i1]),4] (ij52-1)/2]
                #ls123=ls1[i1,:]
                #@show ls123
                #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
                #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
                Al=Fai5(ls123)
                Al2=Fai5(ls1232)
                a=size(Al,1)
                aa22=size(Al2,1)
                ij1=Array{Int32}(undef, aa22, 2)
                if size(Al)[1]!=0
                    if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                        a1=Trans5ott(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
    #                     print(a1)
                        vol1q=v123(ls123,Al,1,2,3,5,jm1,jm2,A1,B1)-v123(ls123,Al,2,3,4,5,jm1,jm2,A1,B1)-v123(ls123,Al,1,2,4,5,jm1,jm2,A1,B1)+v123(ls123,Al,1,3,4,5,jm1,jm2,A1,B1)
    #                     print(vol1q)
                        if ls123[1]!=0
                            jl1=Diagonal([ls123[1]*(ls123[1]+1) for i in 1:size(Al)[1]])
                        else
                            jl1=0
                        end
                        if ls123[2]!=0
                            jl2=Diagonal([ls123[2]*(ls123[2]+1) for i in 1:size(Al)[1]])
                        else
                            jl2=0
                        end
                        if ls123[3]!=0
                            jl3=Diagonal([ls123[3]*(ls123[3]+1) for i in 1:size(Al)[1]])
                        else
                            jl3=0
                        end
                        if ls123[4]!=0
                            jl4=Diagonal([ls123[4]*(ls123[4]+1) for i in 1:size(Al)[1]])
                        else
                            jl4=0
                        end
                        (e1,ev)=eigen(vol1q)
                        av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                        lls1=Array{Float64,2}(undef,aa22,6)
                        lls2=Array{Float64,2}(undef,aa22,6)
                        lls3=Array{Float64,2}(undef,aa22,6)
                        #Initialzation of ls
                        lls1[1:aa22,1].=ls1232[1]
                        lls1[1:aa22,2].=ls1232[2]
                        lls1[1:aa22,3].=view(Al2,1:aa22,2)
                        lls2[1:aa22,1].=view(Al2,1:aa22,2)
                        lls2[1:aa22,2].=ls1232[3]
                        lls2[1:aa22,3].=view(Al2,1:aa22,3)
                        lls3[1:aa22,1].=view(Al2,1:aa22,3)
                        lls3[1:aa22,2].=ls1232[4]
                        lls3[1:aa22,3].=ls1232[5]
                        i2=zeros(aa22)
                        state=zeros(size(a1)[1],a,aa22,Int(2*ls1232[5]+1))*im
                        for ja in 1:Int((2*ls1232[1]+1)*(2*ls1232[2]+1)*(2*ls1232[3]+1)*(2*ls1232[4]+1))
                            mtag=kronmc1(ja,(ls1232[1:4]*2).+1)
#                             println("m=",mtag)
        #                     println(" jnum= ",i1," jlist= ",ls123," msize= ",size(a1)[2]," m= ",ja," mtag= ",mtag)
                                s1 = -mtag[1]-mtag[2]
                                s2 = s1-mtag[3]
                                m5 = s2-mtag[4]
                                if abs(m5)<=ls1232[5]
                                    @view(lls1[1:aa22,4]).=mtag[1]
                                    @view(lls1[1:aa22,5]).=mtag[2]
                                    @view(lls1[1:aa22,6]).=s1
                                    @view(lls2[1:aa22,4]).=-s1
                                    @view(lls2[1:aa22,5]).=mtag[3]  
                                    @view(lls2[1:aa22,6]).=s2
                                    @view(lls3[1:aa22,4]).=-s2
                                    @view(lls3[1:aa22,5]).=mtag[4]  
                                    @view(lls3[1:aa22,6]).=m5
                                    Inter5_Compact(lls1, lls2, lls3, aj1, jm, aa22, s1, s2, ij1, i2)
                                    #println(a[:,ja])
                                        for iia2 in 1:aa22
                                            state[iia2,Int(ls1232[5]+m5+1)]=state[iia2,Int(ls1232[5]+m5+1)]+i2[iia2]*sqrt(2*Al2[iia2,2]+1)*sqrt(2*Al2[iia2,3]+1)
                                        end
                                end
                        end
                            for is in 1:aa22
                                for jml1 in  1:size(a1)[1]
                                    for jml2 in  1:Int(2*ls1232[5]+1)
                                        at1[i1,1]=at1[i1,1]+adjoint(state[jml1,:,is,jml2])*av1*(t1)^(3/2)*state[jml1,:,is,jml2]
                                        at1[i1,2]=at1[i1,2]+adjoint(state[jml1,:,is,jml2])*state[jml1,:,is,jml2]
                                        for ii in 3:nn
                                            at1[i1,ii]=at1[i1,ii]+adjoint(state[jml1,:,is,jml2])*((vol1q*t1^3/8*im)^(ii-2))*state[jml1,:,is,jml2]
                                        end
                                        at1[i1,nn+1]=at1[i1,nn+1]+adjoint(state[jml1,:,is,jml2])*jl1*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+2]=at1[i1,nn+2]+adjoint(state[jml1,:,is,jml2])*jl2*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+3]=at1[i1,nn+3]+adjoint(state[jml1,:,is,jml2])*jl3*state[jml1,:,is,jml2]*t1^2
                                        at1[i1,nn+4]=at1[i1,nn+4]+adjoint(state[jml1,:,is,jml2])*jl4*state[jml1,:,is,jml2]*t1^2
                                    end
                                end
                            end
    #                     ja=1
    #                     at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
    #                     at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
    #                     println(ls123," ",at1[i1,2])
                    else
                        at1[i1,:].=at1[i1,:].+im*0.0
                    end
                else
                    at1[i1,:].=at1[i1,:].+im*0.0
                end
                Al=nothing
                a=nothing
#                 println("j2=",ls1232)
            end
#             println("j1=",ij5)
        end
    end
    return 0
end



function tran6(js1,js2,js3,js4,js5,js6,t1,gf1,gf2,gf3,gf4,gf5,gf6,jm,aj1,center,dist,at1,nn,jm61,jm62,B1,B2)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1),6)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1))
    seq6o(js1,js2,js3,js4,js5,js6,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai6(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans6o(ls123,t1,gf1,gf2,gf3,gf4,gf5,gf6,jm,aj1,Al,size(Al,1))
                vol1q=v123(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v123(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v123(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,3,4,6,jm61,jm62,B1,B2)
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*(t1)^(3/2)*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,1]=im*0.0
                at1[i1,2]=im*0.0
            end
        else
            at1[i1,1]=im*0.0
            at1[i1,2]=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function tran6Mod(js1,js2,js3,js4,js5,js6,t1,gf1,gf2,gf3,gf4,gf5,gf6,jm,aj1,center,dist,at1,jm61,jm62,B1,B2)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1),6)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1))
    seq7oMod(js1,js2,js3,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai6(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans6o(ls123,t1,gf1,gf2,gf3,gf4,gf5,gf6,jm,aj1,Al,size(Al,1))
                vol1q=v123(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v123(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v123(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,3,4,6,jm61,jm62,B1,B2)
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*(t1)^(3/2)*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    at1[i1,3]=at1[i1,3]+transpose(conj(a1[:,ja]))*(vol1q*t1^3/8*im)*a1[:,ja]
                    for ii in 4:nn
                        at1[i1,ii]=at1[i1,ii]+transpose(conj(a1[:,ja]))*(vol1q*t1^3/8*im)^(2*(ii-3))*a1[:,ja]
                    end
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,1]=im*0.0
                at1[i1,2]=im*0.0
            end
        else
            at1[i1,1]=im*0.0
            at1[i1,2]=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function p123(js1,js2,js3,js4,t1,gf1,gf2,gf3,gf4,hf1,hf2,hf3,hf4,jm,aj1,center,dist,at1)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1),4)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1))
    seq4o(js1,js2,js3,js4,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai4(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans4o(ls123,t1,gf1,gf2,gf3,gf4,jm,aj1,Al,size(Al,1))
                a2=Trans4o(ls123,t1,hf1,hf2,hf3,hf4,jm,aj1,Al,size(Al,1))
                vol1q=2*qa(ls123)
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(conj(a2[:,ja]))*av1*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a2[:,ja]))*a1[:,ja]
                    at1[i1,3]=at1[i1,3]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    at1[i1,4]=at1[i1,4]+transpose(conj(a2[:,ja]))*a2[:,ja]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,1]=im*0.0
                at1[i1,2]=im*0.0
            end
        else
            at1[i1,1]=im*0.0
            at1[i1,2]=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end

function p6(js1,js2,js3,js4,js5,js6,t1,gf1,gf2,gf3,gf4,gf5,gf6,hf1,hf2,hf3,hf4,hf5,hf6,jm,aj1,center,dist,at1,jm61,jm62,B1,B2)
    #at1=SharedArray{Complex{Float64}}((js1+1)*(js2+1)*(js3+1)*(js4+1),2)
    ls1=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1),6)
    ls2=zeros((js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1))
    seq6o(js1,js2,js3,js4,js5,js6,ls1,ls2)
    @sync @distributed for i1 in 1:(js1+1)*(js2+1)*(js3+1)*(js4+1)*(js5+1)*(js6+1)
        ls123=ls1[Int(ls2[i1]),:]
        #ls123=ls1[i1,:]
        #@show ls123
        #println(ls123," ",center," ",sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)," ",i1)
        #ls123=[mod(i1,js1+1) mod(floor(i1/(js1+1)),js2+1) mod(floor(i1/(js1+1)/(js2+1)),js3+1) mod(floor(i1/(js1+1)/(js2+1)/(js3+1)),js4+1)]/2
        Al=Fai6(ls123)
        if size(Al)[1]!=0
            if sqrt((ls123[1]-center)^2+(ls123[2]-center)^2+(ls123[3]-center)^2)<=dist
                a1=Trans6o(ls123,t1,gf1,gf2,gf3,gf4,gf5,gf6,jm,aj1,Al,size(Al,1))
                a2=Trans6o(ls123,t1,hf1,hf2,hf3,hf4,hf5,hf6,jm,aj1,Al,size(Al,1))
                vol1q=v123(ls123,Al,1,2,3,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,3,5,6,jm61,jm62,B1,B2)+v123(ls123,Al,1,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,1,2,6,6,jm61,jm62,B1,B2)+v123(ls123,Al,3,4,5,6,jm61,jm62,B1,B2)-v123(ls123,Al,4,5,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,4,6,6,jm61,jm62,B1,B2)-v123(ls123,Al,2,3,4,6,jm61,jm62,B1,B2)
                (e1,ev)=eigen(vol1q)
                av1=ev*diagm((sqrt.(abs.(e1./4))))*inv(ev)
                for ja in 1:size(a1)[2]
                    #println(a[:,ja])
                    at1[i1,1]=at1[i1,1]+transpose(conj(a2[:,ja]))*av1*a1[:,ja]
                    at1[i1,2]=at1[i1,2]+transpose(conj(a2[:,ja]))*a1[:,ja]
                    at1[i1,3]=at1[i1,3]+transpose(conj(a1[:,ja]))*a1[:,ja]
                    at1[i1,4]=at1[i1,4]+transpose(conj(a2[:,ja]))*a2[:,ja]
                end
#                 ja=1
#                 at1[i1,1]=at1[i1,1]+transpose(conj(a1[:,ja]))*av1*a1[:,ja]
#                 at1[i1,2]=at1[i1,2]+transpose(conj(a1[:,ja]))*a1[:,ja]
            else
                at1[i1,1]=im*0.0
                at1[i1,2]=im*0.0
                at1[i1,3]=im*0.0
                at1[i1,4]=im*0.0
            end
        else
            at1[i1,1]=im*0.0
            at1[i1,2]=im*0.0
            at1[i1,3]=im*0.0
            at1[i1,4]=im*0.0
        end
        Al=nothing
        a=nothing
    end
    return 0
end
