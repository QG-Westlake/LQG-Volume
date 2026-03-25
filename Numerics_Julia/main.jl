using WignerSymbols
using LinearAlgebra
setprecision(50)
push!(LOAD_PATH, ".")
include("SL2C.jl")
include("threej.jl")
include("cuba.jl")
include("Trans.jl")
include("volume.jl")

function ME1(g1,g2,g3,g4,t,jc,jr)
    rs=0.0+im*0.0
    for j1 in 2*(jc-jr):2*(jc+jr)
        for j2 in 2*(jc-jr):2*(jc+jr)
            for j3 in 2*(jc-jr):2*(jc+jr)
                for j4 in 2*(jc-jr):2*(jc+jr)
                    rs=rs+transpose(conj(TransformM4([j1 j2 j3 j4]/2,t,g1,g2,g3,g4)))*sqrt(qa([j1 j2 j3 j4]/2))*TransformM4([j1 j2 j3 j4]/2,t,g1,g2,g3,g4)
                end
            end
        end
    end
    return rs
end

function ME2(g1,g2,g3,g4,t,jc,jr)
    rs=0.0+im*0.0
    for j1 in 2*(jc-jr):2*(jc+jr)
        println(j1/2)
        for j2 in 2*(jc-jr):2*(jc+jr)
            print(j2/2)
            for j3 in 2*(jc-jr):2*(jc+jr)
                for j4 in 2*(jc-jr):2*(jc+jr)
                    rs=rs+transpose(conj(TransformM4([j1 j2 j3 j4]/2,t,g1,g2,g3,g4)))*TransformM4([j1 j2 j3 j4]/2,t,g1,g2,g3,g4)
                end
            end
        end
        println()
    end
    return rs
end

function ME3(g1,g2,g3,t,j)
    return transpose(conj(TransformM6f(j,t,g1,g2,g3)))*vq1(j,6)*TransformM6f(j,t,g1,g2,g3)
end

#------------------------------------------------------------------------------------------------
jp=[5 5 5 5 5 5]
jv=[5 5 5 5 5 5]

#------------------------------------------------------------------------------------------------
# 4 vertex SL2C elements
#=
gl1=g1(1/2,-PList[1,:],SU2List[1,:])
gl2=g1(1/2,-PList[2,:],SU2List[2,:])
gl3=g1(1/2,-PList[3,:],SU2List[3,:])
gl4=g1(1/2,-PList[4,:],SU2List[4,:])
=#




#=
#------------------------------------------------------------------------------------------------
#3 flower graph SL2C elements
z1v=[0.0+im*0.0 0.0+im*0.0 0.0+im*jp[1]] #SL2C element for edge 1 in |psi>
z2v=[pi/2+im*0.0 0.0+im*jp[1] 0.0+im*0.0] #SL2C element for edge 2 in |psi>
z3v=[0.0+im*jp[1] pi/2+im*0.0 0.0+im*0.0] #SL2C element for edge 3 in |psi>

w1v=[0.0+im*0.0 0.0+im*0.0 0.0+im*jp[1]] #SL2C element for edge 1 in <psi|
w2v=[pi/2+im*0.0 0.0+im*jp[1] 0.0+im*0.0] #SL2C element for edge 2 in <psi|
w3v=[0.0+im*jp[1] pi/2+im*0.0 0.0+im*0.0] #SL2C element for edge 3 in <psi|

z1=sqrt(conj(z1v)*transpose(z1v))
z2=sqrt(conj(z2v)*transpose(z2v))
z3=sqrt(conj(z3v)*transpose(z3v))
w1=sqrt(conj(w1v)*transpose(w1v))
w2=sqrt(conj(w2v)*transpose(w2v))
w3=sqrt(conj(w3v)*transpose(w3v))
theta=0
chi=0

w1pv=[0; 0; w1]
w2pv=[0; w2*sin(chi); w2*cos(chi)]
z1pv=[0; z1*sin(theta); z1*cos(theta)]
z2pv=[0; 0; z2]




nn=0
t1=0.5
gf1=slce(z1v)
gf2=slce(z2v)
gf3=slce(z3v)
@time ME3(gf1,gf2,gf3,t1,jv)

@time vq1(jv,6)
=#
a=Fai6(jv)[4,:]

@time wigner3j(50, 50, 50, 5, 5, -10)

@time tjst(50, 50, 50, 5, 5, -10)

@time Inter6(jv,a,[0 0 0 0 0 2])

@time Inter6_Optimization(jv,a,[0 0 0 0 0 2])

#println(Fai6(jv))
#println(q123t(jv,Fai6(jv)[1,:],Fai6(jv)[1,:],6))


#println(glc(5,gl3))
#print("Start counting:")
#@time println(TransformM4([4 4 4 4],0.1,gl1,gl2,gl3,gl4))
#@time println(TransformM4p(jv,0.1,gl1,gl2,gl3,gl4))
#@time println(TransformM4(jv,0.1,gl1,gl2,gl3,gl4)-TransformM4p(jv,0.1,gl1,gl2,gl3,gl4))
#println(qa(jv))
#for i in 1:20
#    print(i/2,ME1(gl1,gl2,gl3,gl4,0.1,i/2,0),"\n")
#    print(i/2,ME2(gl1,gl2,gl3,gl4,0.1,i/2,0),"\n")
#end
#ME2(gl1,gl2,gl3,gl4,0.5,4,2)
#nttjsl(10)



#testa=pv([5,5,5,5],pi-asin(2*sqrt(2)/3))
#testb=rlist([pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2])
#@time a = TransformM4([5 5 5 5],0.1,g1(5,testa[3,:],[pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2]),g1(5,testa[3,:],[pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2]),g1(5,testa[3,:],[pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2]),g1(5,testa[3,:],[pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2]))
#println("\n",a)



#=test1a=pv([5,5,5,5],pi-asin(2*sqrt(2)/3))
test2b=rlist([pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2])
test3c=g1(5,testa[3,:],[pi/2,-(pi-asin(2*sqrt(2)/3)),-pi/2])
println(test1a,"\n")
println(test2b,"\n")
println(test3c,"\n")
=#
