using JuMP, LinearAlgebra, DelimitedFiles, Ipopt, CPUTime


function qubo(nome)
#readQ = readdlm("in3000_3.txt")
#readQ = readdlm("be120.3.3")
readQ = readdlm(nome)
#readQ = readdlm("gka/gka5f")
BigM = 50000
#n = 3000
#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0))




count = 1
n = size(readQ,1)
#=
Q = zeros(n,n)
for i=1:n
	for j=1:n
	Q[i,j] = readQ[count]
	count = count +1
	end
end
=#
 Q = -readQ
epsolon_1 = 10^(-3)
lambda = 1
rho = 5
c = 0.3
alpha = 0.98
old_x = rand(n)
old_y = rand(n)
error_phi = (n-sum((old_x[i]-old_y[i])^2 for i=1:n))
println(error_phi, " ", epsolon_1)
	@variable(model, 0<=x[i=1:n]<=1)
	#@variable(model, 0<=y[i=1:n]<=1)

for inter=1:100	#println("to entrando aqui??")
	#println(error_phi, " ", epsolon_1)
	if error_phi >= epsolon_1


	@objective(model, Min, sum(x[i]*x[j]*Q[i,j] for i=1:n, j=1:n) + alpha*sum(x[i]-x[i]^2 for i=1:n)+0.5*lambda*sum(c^2+c*(x[i]-x[i]^2) -c for i=1:n))
	#@constraint(model, sum(x[i] for i=1:n)<= 2500)

	for i=1:n
	set_start_value(x[i], old_x[i])
	end

	JuMP.optimize!(model)
	old_x = JuMP.value.(x)
	#old_y = JuMP.value.(y)
	lambda = lambda*rho

	phi = sum(abs(JuMP.value(x[i])^2-JuMP.value(x[i])) for i=1:n)
	alpha = alpha + lambda*phi
	error_phi = sqrt(phi)
		println( "iteracao ", inter, " Erro cometido ", error_phi)
	end
end

#=

for i=1:n
if JuMP.value(x[i]) > 0.1
	println("indice ", i, " valor x ", JuMP.value(x[i]), " valor y ", JuMP.value(y[i]))
end 
end
=#
obj_value = sum(round(JuMP.value(x[i]))*round(JuMP.value(x[j]))*Q[i,j] for i=1:n, j=1:n)
#println(Q)

f = open("resultados_x(x-1)_ipopt.txt", "a+")
println(f, nome, " ", n, " ",rho," ", c, " ",lambda, " ",epsolon_1," ", error_phi, " ", obj_value, " ", maximum(JuMP.value.(x)), " ", minimum(JuMP.value.(x)))
close(f)

return obj_value 
end

function main()
nomes_instancias = ["gka/gka10b", "gka/gka1d", "gka/gka2c", "gka/gka3b", "gka/gka4a", "gka/gka4f", "gka/gka5e", "gka/gka6d", "gka/gka8a",
"gka/gka10d", "gka/gka1e", "gka/gka2d", "gka/gka3c", "gka/gka4b", "gka/gka5a", "gka/gka5f", "gka/gka7a", "gka/gka8b",
"gka/gka1a", "gka/gka2e", "gka/gka3d", "gka/gka4c", "gka/gka5b", "gka/gka6a", "gka/gka7b", "gka/gka8d",
"gka/gka1b", "gka/gka2a", "gka/gka2f", "gka/gka3e", "gka/gka4d", "gka/gka5c", "gka/gka6b", "gka/gka7c", "gka/gka9b",
"gka/gka1c", "gka/gka2b", "gka/gka3a", "gka/gka3f", "gka/gka4e", "gka/gka5d", "gka/gka6c", "gka/gka7d", "gka/gka9d", "gka/gka1f"]
g = open("tempo_execucao_x(x-1)_ipopt.txt", "w+")
for i=1:length(nomes_instancias)
for k=1:200
tempo = @time @CPUtime qubo( nomes_instancias[i])

println(g, nomes_instancias[i], " ", tempo )
end
end
close(g)
end

main()
