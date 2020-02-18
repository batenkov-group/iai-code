include("./superres.jl")

using Plots
pyplot();

let ntests=15000,seed=5,Ωₘ=100
    let ℓ=2,R=8,clusterFun=Ζ,amplitudeFun = imag1s
        let pertFun=𝔉ᵣ(Ωₘ,seed), ni = 6
            let ϵR = 10.0.^LinRange(-4,3,ntests), ΩR = LinRange(20,40,ntests), ΔR = LinRange(1e-4,2e-3,ntests)#10.0.^LinRange(-5,-1,ntests)
#            let ϵR = 10.0.^LinRange(-2,1,ntests), ΩR = NrangeLog(1.2,2,ntests), ΔR = 10.0.^LinRange(-5,-1,ntests)
                let res=mapreduce(transpose,vcat,
                        [runSingle_ind(ℓ,R,clusterFun, pertFun, amplitudeFun, ni, ϵR,ΩR,ΔR) for _ in 1:ntests])
                    let clr = map(cc,res[:,2]), sh = map(ss,res[:,2]),idxsSucc = findall(res[:,2].==1),idxsFail=findall(res[:,2].==0)
                        #p=scatter(log10.(res[:,1]),log10.(inv.(res[:,3])),
                        #    mc=clr,msc=:white,ms=7,shape=sh,markerstrokewidth=0.4,markeralpha=0.5,
                        #    xlabel="log SRF",ylabel="log 1/ϵ",lab="",color_palette=:magma)
                        p=scatter(log10.(res[idxsSucc,1]),log10.(inv.(res[idxsSucc,3])),
                            mc=:blue,msc=:white,ms=7,shape=:dtriangle,markerstrokewidth=0.1,markeralpha=0.1,
                            xlabel="log SRF",ylabel="log 1/ϵ",lab="Success",color_palette=:magma)
                        scatter!(log10.(res[idxsFail,1]),log10.(inv.(res[idxsFail,3])),
                            mc=:red,msc=:white,ms=7,shape=:circle,markerstrokewidth=0.1,markeralpha=1.0,
                            xlabel="log SRF",ylabel="log 1/ϵ",lab="Fail",color_palette=:magma)
                        #plot!([t->(2*ℓ-2)*t+5,t->(2*ℓ-1)*t+5],w=3,l=:dash,lc=[:black :blue],lab=["SRF^(2-2ℓ)" "SRF^(1-2ℓ)"])
                        title!("p=$(ℓ),d=$(R), j=$(ni), \$ {\\bf S_1} \$")
                        #ylims!(-.3,0.3)
                        xlims!(1.5,2.5)
                        savefig(p,"./figures/fig6.pdf")
                        #display(sum(res[:,2].==0))
                    end
                end
            end
        end
    end
end
