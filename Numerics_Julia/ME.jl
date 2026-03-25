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