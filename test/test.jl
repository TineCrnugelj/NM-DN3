using DN3
using Plots

# Diferencialna enacba
perioda_f(t, x) = [x[2], 4(1 - x[1]^2)*x[2] - x[1]]


# zacetna priblizka
x0 = 2.0
y0 = 0.0

periodaZacetniProblem = ZacetniProblemNDE(perioda_f, 0.0, [x0, y0])

# resitev s trapezno metodo
z_tr, t_tr = resiTrapez(periodaZacetniProblem, 40; stKorakov=20000)

plot(z_tr[:,1], z_tr[:,2])


# nicle
tocke = []
for i=1:length(t_tr)-1
    if z_tr[i, 1] * z_tr[i+1, 1] < 0
        push!(tocke, (t_tr[i] + t_tr[i+1])/2)
    end
end


#rezultat 10.2035236
perioda = round(tocke[end] - tocke[end-2]; digits = 10)


