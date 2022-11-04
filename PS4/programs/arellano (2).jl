include("deuda.jl")
abstract type Default <: Deuda end

struct Arellano <: Default
	pars::Dict{Symbol, Float64}

	bgrid::Vector{Float64}
	ygrid::Vector{Float64}
	Py::Matrix{Float64}

	v::Matrix{Float64}
	vR::Matrix{Float64}
	vD::Matrix{Float64}
	prob::Matrix{Float64}

	gc::Array{Float64, 3}
	gb::Array{Float64, 2}

	q::Matrix{Float64}
end

function Arellano(;
	β = 0.953,
	γ = 2,
	r = 0.017,
	ψ = 0.282,
	ψ0 = 0.282,
	ψ1 = 1,
	χ = 0.01,

	Δ = 0.1,
	defcost_OG = 1,

	ρy = 0.945,
	σy = 0.025,

	Nb = 200,
	Ny = 21,

	bmax = 0.6
	)

	pars = Dict(:β=>β, :γ=>γ, :r=>r, :ψ=>ψ, :ψ0 => ψ0, :ψ1 => ψ1, :χ=>χ, :Δ=>Δ, :ρy=>ρy, :σy=>σy, :defcost_OG => defcost_OG)

	ychain = tauchen(Ny, ρy, σy, 0, 2)

	Py = ychain.p
	ygrid = exp.(ychain.state_values)

	bgrid = range(0, bmax, length=Nb)

	v = zeros(Nb, Ny)
	vR = zeros(Nb, Ny)
	vD = zeros(Nb, Ny)
	prob = zeros(Nb, Ny)

	gc = zeros(Nb, Ny, 2)
	gb = zeros(Nb, Ny)

	q = ones(Nb, Ny)

	return Arellano(pars, bgrid, ygrid, Py, v, vR, vD, prob, gc, gb, q)
end

function logsumexp(a::AbstractVector{<:Real})
	m = maximum(a)
	return m + log.(sum(exp.(a .- m)))
end

function defcost(yv, dd::Arellano)
	if haskey(dd.pars, :defcost_OG) && dd.pars[:defcost_OG] == 1
		return defcost_OG(yv)
	else
		return defcost_lineal(yv, dd)
	end
end

defcost_OG(yv) = ifelse(yv <= 0.969, yv, 0.969)

defcost_lineal(yv, dd::Default) = defcost_lineal(yv, dd.pars[:Δ])
defcost_lineal(yv, Δ::Number) = yv * (1-Δ)

debtprice(dd::Default, bpv, yv, itp_q) = itp_q(bpv, yv)

function value_default(jb, jy, dd::Default)
	""" Calcula el valor de estar en default en el estado (b,y) """
	β, ψ, ψ0, ψ1 = (dd.pars[sym] for sym in (:β, :ψ, :ψ0, :ψ1))
	yv = dd.ygrid[jy]
	bv = dd.gb[jb, jy]

	# Consumo en default es el ingreso menos los costos de default
	c = defcost(yv, dd)

	# ψ es una función exógena (primera versión)
	#ψ = ψ0 + ψ1*(yv - mean(dd.ygrid))

	# ψ es una función del monto de deuda defaulteado (segunda versión)
	ψ = max(0, ψ0 - ψ1 * bv)


	# Valor de continuación tiene en cuenta la probabilidad ψ de reacceder a mercados
	Ev = 0.0
	for jyp in eachindex(dd.ygrid)
		prob = dd.Py[jy, jyp]
		Ev += prob * ( ψ * dd.v[jb, jyp] + (1-ψ) * dd.vD[jb, jyp] )
	end

	v = u(c, dd) + β * Ev
	
	return c, v
end

function vfi_iter!(new_v, itp_q, dd::Default)
	# Reconstruye la interpolación de la función de valor
	itp_v = make_itp(dd, dd.v)

	for jy in eachindex(dd.ygrid)
		for jb in eachindex(dd.bgrid)
		
			# En repago
			vp, c_star, b_star = opt_value(jb, jy, itp_q, itp_v, dd)

			# Guarda los valores para repago 
			dd.vR[jb, jy] = vp
			dd.gb[jb, jy] = b_star
			dd.gc[jb, jy, 1] = c_star

			# En default
			cD, vD = value_default(jb, jy, dd)
			dd.vD[jb, jy] = vD
			dd.gc[jb, jy, 2] = cD		

		end

	end

	χ = dd.pars[:χ]
	for jb in eachindex(dd.bgrid), jy in eachindex(dd.ygrid)
		# Valor de repagar y defaultear llegando a (b,y)
		vr = dd.vR[jb, jy]
		vd = dd.vD[jb, jy]

		# Probabilidad de default
		## Modo 1: valor extremo tipo 1 directo
		# pr = exp(vd / χ) / ( exp(vd / χ) + exp(vr / χ) )
		# V = χ * log( exp(vd/χ) + exp(vr/χ) )

		## Modo 2: valor extremo tipo X evitando comparar exponenciales de cosas grandes
		lse = logsumexp([vd/χ, vr/χ])
		lpr = vd / χ - lse
		pr = exp(lpr)
		V = χ * lse

		# Guarda el valor y la probabilidad de default al llegar a (b,y)
		new_v[jb, jy] = V
		dd.prob[jb, jy] = pr
	end
end

function q_iter!(new_q, dd::Arellano)
	""" Ecuación de Euler de los acreedores determinan el precio de la deuda dada la deuda, el ingreso, y el precio esperado de la deuda """
	r = dd.pars[:r]

	for jbp in eachindex(dd.bgrid), jy in eachindex(dd.ygrid)
		Eq = 0.0
		for jyp in eachindex(dd.ygrid)
			prob = dd.Py[jy, jyp]
			prob_def = dd.prob[jbp, jyp]

			# Si el país tiene acceso a mercados, emite y puede hacer default mañana
			rep_R = (1-prob_def) * 1 + prob_def * 0

			Eq += prob * rep_R
		end
		new_q[jbp, jy]  = Eq / (1+r)
	end
end

function update_q!(dd::Default; tol = 1e-8, maxiter = 2000, verbose = false)
	""" Itera sobre la ecuación de Bellman de los acreedores para encontrar el precio de la deuda """
	new_q = copy(dd.q)

	dist = 1+tol
	iter = 0
	while dist > tol && iter < maxiter
		iter += 1

		q_iter!(new_q, dd)

		dist = norm(new_q - dd.q) / max(1,norm(dd.q))

		dd.q .= new_q
	end
	verbose && print(iter)
end

function eqm!(dd::Default; tol = 1e-8, maxiter = 250, verbose=true)
	""" Itera sobre la mejor respuesta del país y los acreedores hasta encontrar políticas de default y consumo óptimas dados precios de la deuda que a su vez reflejen la probabilidad de default """
	dist = 1+tol
	iter = 0

	q_old = similar(dd.q)
	v_old = similar(dd.v)

	while dist > tol && iter < maxiter
		iter += 1
		print("Iteration $iter: ")

		q_old .= dd.q
		v_old .= dd.v

		# Problema del país dados los precios de la deuda
		vfi!(dd, verbose = false)
		dist_v = norm(dd.v - v_old) / max(1,norm(v_old))

		# Problema de los acreedores dada la probabilidad de default del país
		update_q!(dd, verbose = false)

		dist_q = norm(dd.q - q_old) / max(1,norm(q_old))
		dist = max(dist_q, dist_v)

		verbose && print("dist = $dist\n")
	end
end

function make_itp(dd::Default, y::Array{Float64,2})
	@assert size(y) == (length(dd.bgrid), length(dd.ygrid))

	knts = (dd.bgrid, dd.ygrid)
	interpolate(knts, y, Gridded(Linear()))
end

function mpe!(dd::Default; tol = 1e-8, maxiter = 500)
	new_v = similar(dd.v)
	new_q = similar(dd.q)

	dist = 1+tol
	iter = 0

	while dist > tol && iter < maxiter
		iter += 1

		print("Iteration $iter: ")

		# Actualiza el precio de la deuda
		q_iter!(new_q, dd)
		dist_q = norm(new_q - dd.q) / max(1, norm(dd.q))
		
		# Interpolación del precio de la deuda
		itp_q = make_itp(dd, new_q)

		# Actualiza la función de valor
		vfi_iter!(new_v, itp_q, dd)
		dist_v = norm(new_v - dd.v) / max(1, norm(dd.v))

		# Distancias
		dist = max(dist_q, dist_v)

		# Guardamos todo
		dd.v .= new_v
		dd.q .= new_q

		print("dist (v,q) = ($(round(dist_v, sigdigits=2)), $(round(dist_q, sigdigits=2)))\n")
	end
end

