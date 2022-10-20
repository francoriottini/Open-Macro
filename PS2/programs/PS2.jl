include("itpcake.jl")

#=funciones de consumo =#

ce = CakeEating(Nk = 20);
vfi!(ce)
vfi_itp!(ce)

#B = [i for i in range(1,200,200)]

plot_c= scatter(x = ce.kgrid, y = ce.gc)

plot1 = plot([plot_c],  Layout(title="Modelo de la torta: consumo en función del capital", xaxis_title="torta", yaxis_title = "consumo"))

savefig(plot1, "Consume_Capital_vfi.png")

fracc_consumo = ce.gc ./ ce.kgrid

plot_fc = scatter(x = ce.kgrid, y = fracc_consumo)

plot2 = plot([plot_fc], Layout(title="Modelo de la torta: fracción consumida en función del capital", xaxis_title="torta", yaxis_title = "fracción consumida"))

savefig(plot2, "Fraccion_consumida.png")

fracc_ahorro = ce.gk ./ ce.kgrid

plot_fc2 = scatter(x = ce.gk, y = fracc_ahorro)

plot3 = plot([plot_fc2], Layout(title = "Modelo de la torta: ahorro en función del capital", xaxis_title = "torta", yaxis_title = "fracción de torta ahorrada"))

savefig(plot3, "Fraccion_ahorrada.png")

#= Simulador de torta =#

function CakeSimul!(ce::CakeEating, T::Int64 = 1000)

    C = Vector(undef, T)
    K = Vector(undef, T)

    for i in 1:T
        ce = CakeEating(Nk = 20)
        vfi_itp!(ce)
        C[i] = ce.gc
        K[i] = ce.gk
    end
    plot_1 = scatter(x = T, y = C)
    plot_2 = scatter(x = T, y = K)

    return plot([plot_1, plot_2])
end