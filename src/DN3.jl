module DN3

export  ZacetniProblemNDE, resiTrapez

"""
ZacetniProblemNDE(f, x0, y0)
    
Podatkovni tip za zacetni problem navadne diferencialne enacbe. 
`f` predstavlja funkcijo desnih strani DE y' = f(x, y).
`x0` je zacetna vrednost neodvisne spremenljivke.
`y0` je zacetna vrednost odvisne spremenljivke.
"""

struct ZacetniProblemNDE
    f 
    x0 
    y0 
end

"""
y, x = resiTrapez(zp::ZacetniProblemNDE, xk; n = n)

Resi zacetni problem za NDE s trapezno metodo na [zp.x0, xk] s fiksnim korakom.
Funkcija vzame zacetni problem in izracuna resitev v diskretnih tockah med x0 in xk z Eulerjevo metodo.
"""

function resiTrapez(zacPr::ZacetniProblemNDE, xk; stKorakov = 100)
    y = zeros(stKorakov+1, length(zacPr.y0))
    y[1,:] = zacPr.y0
    x = LinRange(zacPr.x0, xk, stKorakov+1)
    h = x[2] - x[1]

    for i=1:stKorakov
        # Eulerjeva metoda
        y[i+1, :] = y[i, :] + h * zacPr.f(x[i], y[i, :])
        y[i+1, :] = y[i, :] + h / 2 * (zacPr.f(x[i], y[i, :]) + zacPr.f(x[i+1], y[i+1, :]))  
    end
    
    return y, x
end


end # module DN3