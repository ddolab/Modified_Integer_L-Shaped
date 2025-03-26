using CairoMakie
using DataFrames
using JLD
using Typst_jll
ILS_results = load("75_15_64_24_24_Oct_2024_03-56-00results.jld")
MILS_results = load("75_15_64_24_23_Oct_2024_17-59-20results.jld")
inch = 96
pt = 4/3
f= Figure(fontsize=15pt)
#f = Figure(size=(6.5*inch, 2.3214285714*inch),fontsize=9*pt)
#f = Figure(size=(1600, 1200))
ax2=  Axis(f[1, 1], xlabel="Time (s)", ylabel="Objective Value", title="Convergence Plot")
ax1 = Axis(f[1, 2], limits=(nothing, nothing, 1.39e7, 1.43e7), xlabel="Time (s)", ylabel="Objective Value", title="Lower Bounds and Mixed-Integer Subproblems")
colsize!(f.layout, 1, Aspect(1, 4/3))
colsize!(f.layout, 2, Aspect(1, 4/3))
x1 = ILS_results["boundstimes"][1:end-3]
push!(x1, 10800)
y1 = -1 * ILS_results["boundsvec"][1:end-3]
push!(y1, -1 * ILS_results["boundsvec"][end-3])
x2 = ILS_results["incumbenttimes"]
push!(x2, 10800)
y2 = -1 * ILS_results["incumbentobjvals"]
push!(y2, -1 * ILS_results["incumbentobjvals"][end])

#scatter!(f[1, 1], x1[1:end-1], y1[1:end-1]; color= 2,colorrange=(1,4), markersize=15, label="ILS")
ilsub=stairs!(f[1, 1], x2, y2; step=:post, color= 2,colorrange=(1,4), label="ILS upper bound",linestyle=:dash) 
full_optpoints1 = findall(x -> x == "full_opt", ILS_results["cuttypes"][1:end-3])
Benderspoints1 = findall(x -> x == "Benders", MILS_results["cuttypes"])
scatter!(f[1, 2], x1[full_optpoints1], y1[full_optpoints1]; color= 1,colorrange=(1,4), marker=:xcross, markersize=15, label="ILS")
stairs!(f[1, 2], x2, y2; step=:post, color= 2,colorrange=(1,4), label="ILS")

x3 = MILS_results["boundstimes"]

y3 = -1 * MILS_results["boundsvec"]

x4 = MILS_results["incumbenttimes"]
push!(x4, x3[end])
y4 = -1 * MILS_results["incumbentobjvals"]
push!(y4, -1 * MILS_results["incumbentobjvals"][end])
ilslb=stairs!(f[1, 1], x1, y1; step=:post, color= 2,colorrange=(1,4), label="ILS lower bound")
milsub=stairs!(f[1, 1], x4, y4; step=:post, color= 4,colorrange=(1,4), label="Mod. ILS upper bound",linestyle=:dash)
milslb=stairs!(f[1, 1], x3, y3; step=:post, color= 4,colorrange=(1,4), label="Mod. ILS lower bound")
stairs!(f[1, 2], x3, y3; step=:post, color= 4,colorrange=(1,4), label="Mod. ILS")
stairs!(f[1, 2], x1, y1; step=:post, color= 2,colorrange=(1,4), label="ILS")
#scatter!(f[1,1], x3, y3; color = :red, markersize = 10, label = "Mod. ILS")
stairs!(f[1, 2], x4, y4; step=:post, color= 4,colorrange=(1,4), label="Mod. ILS")
Benderspoints = findall(x -> x == "Benders", MILS_results["cuttypes"])
suboptpoints = findall(x -> x == "subopt", MILS_results["cuttypes"])
full_optpoints = findall(x -> x == "full_opt", MILS_results["cuttypes"])

#scatter!(f[1,1], x3[Benderspoints], y3[Benderspoints]; color = :red, markersize = 10, label = "Mod. ILS")
#scatter!(f[1, 1], x3[suboptpoints], y3[suboptpoints]; color=3,colorrange=(1,4), marker=:circle, markersize=15, label="Mod. ILS")
#scatter!(f[1, 1], x3[full_optpoints], y3[full_optpoints]; color= 1,colorrange=(1,4), marker=:xcross, markersize=15, label="Mod. ILS")
circ=scatter!(f[1, 2], x3[suboptpoints], y3[suboptpoints]; color=3,colorrange=(1,4), marker=:circle, markersize=15, label="Mod. ILS")
cross=scatter!(f[1, 2], x3[full_optpoints], y3[full_optpoints]; color= 1,colorrange=(1,4), marker=:xcross, markersize=15, label="Mod. ILS")
colgap!(f.layout, 1, 100)
axislegend(ax2)
axislegend(ax1,[circ,cross],[ "Suboptimal MISP", "Optimal MISP"],"Source of Cuts",position= :lt)
#resize_to_layout!(f)

save("figure.png", f)
save("figure.svg", f)