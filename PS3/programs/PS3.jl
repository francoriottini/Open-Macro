include("deuda.jl") #Correr solo si no se va a necesitar el theta
include("deudatheta.jl") #Correr solo si se va a necesitar el theta


using Colors, ColorSchemes

nd = NoDefault()
vfi!(nd)

cons = Vector(undef, length(nd.ygrid))

cmap = ColorScheme(range(colorant"red", colorant"orange", length=21))

cons2 = [scatter(x=nd.gb[1:size(nd.gb)[1], i], y=nd.gc[1:size(nd.gc)[1], i], name= round(nd.ygrid[i], digits = 2), line_color = cmap[i]) for i in 1:length(nd.ygrid)]

layout = Layout(title="Consume as a function of debt in diferents incomes",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.02),
xaxis_title="Debt", 
yaxis_title = "consume",
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
legend_title_text = "y =",
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=12,
font_family="Lato",
hovermode="x",
)

plot1 = plot(cons2, layout)

savefig(plot1, "consumefunction.png")

#=====================Opcional de consumo ==============================#
cons3 = [scatter(x=nd.gb[1:size(nd.gb)[1], i], y= nd.gc[1:size(nd.gc)[1], i]/nd.ygrid[i], name= round(nd.ygrid[i], digits = 2), line_color = cmap[i]) for i in 1:length(nd.ygrid)]

layout = Layout(title="Consume over income as a function of debt in diferents incomes",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.02),
xaxis_title="Debt", 
yaxis_title = "consume over income",
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
legend_title_text = "y =",
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=12,
font_family="Lato",
hovermode="x",
)

plot2 = plot(cons3, layout)
savefig(plot2, "consumeoverincomefunction.png")

#==================Estaticas comparadas=============#
Q = [i for i in range(nd.pars[:β],0.99,10)]

pv = Vector{AbstractTrace}(undef, 10)

new_r = Vector(undef, 10)
for i in 1:length(Q)
    new_r[i] = (1/Q[i]) - 1
end

nd = NoDefault(r = new_r[1])
vfi!(nd);

gcvsincome = scatter(x = nd.ygrid, y = nd.gc[1,1:21], line_color = "#0098e9")


layout = Layout(title="Consume as a function of income",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.05),
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
legend_title_text = "q = 0.96",
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=16,
font_family="Lato",
hovermode="x",
)

plot(gcvsincome, layout)

#plot3 = plot(gcscatters, layout)
#savefig(plot3, "consumeandincome.png")

#==================Estaticas comparadas completo=======================#
Q = [i for i in range(nd.pars[:β],0.99,10)]

pv = Vector{AbstractTrace}(undef, 10)

new_r = Vector(undef, 10)
for i in 1:length(Q)
    new_r[i] = (1/Q[i]) - 1
end

for i in 1:length(pv)
    nd = NoDefault(r = new_r[i])
    vfi!(nd, verbose = false)
    leyenda = round(Q[i], digits = 3)
    pv[i] = scatter(x=nd.ygrid, y= nd.gc[1,1:21], line_color = cmap[i], name = "q = $leyenda")
end

layout = Layout(title="Consume as a function of income with differents interest rates (r = (1/q) -1)",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.02),
xaxis_title="Income", 
yaxis_title = "c(0,y)",
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=12,
font_family="Lato",
hovermode="x",
)

plot4 = plot(pv, layout)

savefig(plot4, "consumeandincome.png")

#=========================Movimientos en la varianza del ingreso ==================#
Σ=[i for i in range(0.01, 0.05, 10)]

pv2 = Vector{AbstractTrace}(undef, 10)

for i in 1:length(pv2)
    nd = NoDefault(σy = Σ[i]);
    vfi!(nd, verbose = false);
    leyenda = round(Σ[i], digits = 3);
    pv2[i] = scatter(x=nd.ygrid, y= nd.gc[1,:], line_color = cmap[i], name = "σy = $leyenda")
end

layoutsigma = Layout(title="Consume as a function of income with differents income standar deviations (σy)",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.02),
xaxis_title="Income", 
yaxis_title = "c(0,y)",
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=12,
font_family="Lato",
hovermode="x",
)

plot5 = plot(pv2, layoutsigma)

savefig(plot5, "consumeandincome_sigma.png")


#=================BONUS===================#
nd = NoDefault()
vfi!(nd, verbose = false)

cons = Vector{AbstractTrace}(undef, length(nd.ygrid))

for i in 1:length(nd.ygrid)
    cons[i] = scatter(x=nd.gc[1:size(nd.gc)[1], i], y=nd.ygrid, name= round(nd.ygrid[i], digits = 2), line_color = cmap[i])
end

my_colors = ["#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF",
"#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF",
"#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF"]

cmap = ColorScheme(range(colorant"red", colorant"orange", length=21))

layouttitazero = Layout(title="Consume as a function of income with θ = 0",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.05),
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
legend_title_text = "y =",
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=16,
font_family="Lato",
hovermode="x",
)

plot_θzero = plot(cons, layouttitazero)

savefig(plot_θzero, "consumeandincome_titazero.png")


#=============================Cambiamos\theta ===#

nd2 = NoDefault()
vfi!(nd2)

constita = Vector{AbstractTrace}(undef, length(nd2.ygrid))

for i in 1:length(nd2.ygrid)
    constita[i] = scatter(x=nd2.gc[1:size(nd2.gc)[1], i], y=nd2.ygrid, name= round(nd2.ygrid[i], digits = 2), line_color = cmap[i])
end

my_colors = ["#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF",
"#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF",
"#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF", "#357EBDFF", "#9632B8FF", "#B8B8B8FF"]

cmap = ColorScheme(range(colorant"red", colorant"orange", length=21))

layouttheta = Layout(title="Consume as a function of income with θ = 1",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.05),
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
legend_title_text = "y =",
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=16,
font_family="Lato",
hovermode="x",
)

plot(constita, layouttheta)

#=============================#
Θ=[i for i in range(1, 5, 10)]

pvtheta = Vector{AbstractTrace}(undef, 10)

for i in 1:length(pvtheta)
    nd = NoDefault(θ = Θ[i]);
    vfi!(nd, verbose = false);
    leyenda = round(Θ[i], digits = 3);
    pvtheta[i] = scatter(x=nd.ygrid, y= nd.gc[1,:], line_color = cmap[i], name = "θ = $leyenda")
end

layoutttheta = Layout(title="Consume as a function of income with differents robustness preferences (θ)",
width=1920 * 0.5, height=1080 * 0.5,
legend=attr(orientation="h", x=0.02),
xaxis_title="Income", 
yaxis_title = "c(0,y)",
xaxis=attr(zeroline=false, gridcolor="#434343"),
yaxis=attr(zeroline=false, gridcolor="#434343"),
paper_bgcolor="#272929", plot_bgcolor="#272929",
font_color="#F8F1E9", font_size=12,
font_family="Lato",
hovermode="x",
)

plottheta = plot(pvtheta, layoutttheta)

savefig(plottheta, "consumeandincomerobustness.png")