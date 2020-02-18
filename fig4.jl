include("./superres.jl")

using Plots
pyplot();

let ntests=100,ℓ=2,R=3,clusterFun=Ζ,amplitudeFun = opposite1s
    let pertFun=𝔉ᵪ(ℓ)
        let ϵR = 10.0.^LinRange(-14,-3,150), ΩR = NrangeLog(1,3,150), ΔR = 10.0.^LinRange(-4,-2,150)
            let res1 = [runSingle_Acc(ℓ,R,clusterFun, pertFun,amplitudeFun, ϵR,ΩR,ΔR) for _ in 1:ntests]
                let res=mapreduce(transpose,vcat,res1)
                    let ΔX=reduce(vcat,res[:,2]), ΔC = reduce(vcat,res[:,3]),Succ=reduce(vcat,res[:,4])
                        p1=scatter(log10.(res[:,1]),log10.(ΔX./Succ),
                                    shape=[:circle :x :square],color=[:red :black :blue],
                                    markeralpha=[0.6 1.0 0.6],markerstrokewidth=[0.1 1 0.1],
                                    msc=[:white :black :white],ms=[8 5 5],
                                    xlabel="log(SRF)",ylabel="log(\$ \\cal K \$)",
                                    lab=cat(["MP \$ {\\cal K}_{x,$(i)} \$" for i=1:R]...,dims=2));
                        plot!([t->(2*ℓ-2)*t+1],w=1,l=:dash,lc=[:black],lab="\$ SRF^{2p-2} \$ ");
                        p2=scatter(log10.(res[:,1]),log10.(ΔC./Succ),
                                shape=[:circle :x :square],color=[:red :black :blue],
                                markeralpha=[0.6 1.0 0.6],markerstrokewidth=[0.1 1 0.1],
                                msc=[:white :black :white],ms=[8 5 5],
                                xlabel="log(SRF)",ylabel="log(\$ \\cal K \$)",
                                lab=cat(["MP \$ {\\cal K}_{a,$(i)} \$" for i=1:R]...,dims=2));
                        plot!([t->(2*ℓ-1)*t+1],w=1,l=:dash,lc=[:black],lab="\$ SRF^{2p-1} \$ ");
                        titleStr = "p=$(ℓ),d=$(R), \$ {\\bf S_2} \$"
                        p=plot(p1,p2,layout=(1,2),size=(800,300),title=titleStr);
                        savefig(p,"./figures/fig4.pdf")
                    end
                end
            end
        end
    end
end
