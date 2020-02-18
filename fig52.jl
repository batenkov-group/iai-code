include("./superres.jl")

using Plots
pyplot();

let ntests=15000,seed=0,Ωₘ=2000
    let ℓ=3,R=4,clusterFun=Ζ, amplitudeFun = imag1s
        let pertFun=𝔉ᵣ(Ωₘ,seed)
            let ϵR = 10.0.^LinRange(-14,-5,150), ΩR = NrangeLog(2,3,150), ΔR = 10.0.^LinRange(-4,-2,150)
                let res = mapreduce(transpose,vcat,
                        [runSingle(ℓ,R,clusterFun, pertFun,amplitudeFun, ϵR,ΩR,ΔR) for _ in 1:ntests])
                    let clr = map(cc,res[:,2]), sh = map(ss,res[:,2]),idxsSucc = findall(res[:,2].==1),idxsFail=findall(res[:,2].==0)
                        #p = scatter(log10.(res[:,1]),log10.(inv.(res[:,3])),
                        #    mc=clr,msc=:white,ms=7,shape=sh,markerstrokewidth=0.4,markeralpha=0.5,
                        #    xlabel="log SRF",ylabel="log 1/ϵ",
                        #    color_palette=:magma,lab="")
                        p=scatter(log10.(res[idxsSucc,1]),log10.(inv.(res[idxsSucc,3])),
                            mc=:blue,msc=:white,ms=7,shape=:dtriangle,markerstrokewidth=0.1,markeralpha=0.1,
                            xlabel="log SRF",ylabel="log 1/ϵ",lab="Success",color_palette=:magma)
                        scatter!(log10.(res[idxsFail,1]),log10.(inv.(res[idxsFail,3])),
                            mc=:red,msc=:white,ms=7,shape=:circle,markerstrokewidth=0.1,markeralpha=1.0,
                            xlabel="log SRF",ylabel="log 1/ϵ",lab="Fail",color_palette=:magma)
                        #res[:,2]
                        plot!(t->(2*ℓ-1)*t+4,w=2,l=:dash,lc=:black,
                            lab="\$ SRF^{2p-1}\$")
                        ylims!(5,13)
                        xlims!(0,2)
                        title!("p=$(ℓ),d=$(R), \$ {\\bf S_1} \$")
                        #savefig(p,"../figures/phase_tran_5.pdf")
                    end
                end
            end
        end
    end
end
