using JuMP, LinearAlgebra, DelimitedFiles, Ipopt, CPUTime


function qubo(nome)

readQ = readdlm(nome)
BigM = 50000
#n = 3000
#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0))
#model = Model(Alpine.Optimizer)
count = 1

count = 1
n = size(readQ,1)

Q = readQ
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
	@variable(model, 0<=y[i=1:n]<=1)

for inter=1:100	#println("to entrando aqui??")
	#println(error_phi, " ", epsolon_1)
	if error_phi >= epsolon_1


	@objective(model, Min, -sum(x[i]*x[j]*Q[i,j] for i=1:n, j=1:n) + alpha*(n-sum((x[i]-y[i])^2 for i=1:n))+0.5*lambda*sum(c^2+c*(1-(x[i]-y[i])^2) -c for i=1:n))
	#@constraint(model, sum(x[i] for i=1:n)<= 2500)

	for i=1:n
	set_start_value(x[i], old_x[i])
	set_start_value(y[i], old_y[i])
	end

	JuMP.optimize!(model)
	old_x = JuMP.value.(x)
	old_y = JuMP.value.(y)
	lambda = lambda*rho

	phi = n-sum((old_x[i]-old_y[i])^2 for i=1:n)
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
obj_value = -sum(round(JuMP.value(x[i]))*round(JuMP.value(x[j]))*Q[i,j] for i=1:n, j=1:n)
#println(Q)

f = open("resultados_(x-y)_be.txt", "a+")
println(f, nome, " ", n, " ",rho," ", c, " ",lambda, " ",epsolon_1," ", error_phi, " ", obj_value, " ", maximum(JuMP.value.(x)), " ", minimum(JuMP.value.(x)))
close(f)

return obj_value 
end

function main()
nomes_instancias = ["be/be100.10","be/be100.2","be/be100.3","be/be100.4","be/be100.5","be/be100.6","be/be100.7","be/be100.1","be/be100.8","be/be100.9", "be/be120.3.10" ,"be/be120.8.4" ,"be/be150.3.8" ,"be/be200.3.2" ,"be/be200.8.6" ,"be/be120.3.1" ,"be/be120.8.5" ,"be/be150.3.9" ,"be/be200.3.3" ,
"be/be200.8.7" ,"be/be120.3.2" ,"be/be120.8.6" ,"be/be150.8.10" ,"be/be200.3.4" ,"be/be200.8.8" ,"be/be120.3.3" ,"be/be120.8.7" ,"be/be150.8.1" ,
"be/be200.3.5" ,"be/be200.8.9" ,"be/be120.3.4" ,"be/be120.8.8" ,"be/be150.8.2" ,"be/be200.3.6" ,"be/be250.10" ,"be/be120.3.5" ,"be/be120.8.9" ,"be/be150.8.3" ,
"be/be200.3.7" ,"be/be250.1" ,"be/be120.3.6" ,"be/be150.3.10" ,"be/be150.8.4" ,"be/be200.3.8" ,"be/be250.2" ,"be/be120.3.7" ,"be/be150.3.1" ,"be/be150.8.5" ,
"be/be200.3.9" ,"be/be250.3" ,"be/be120.3.8" ,"be/be150.3.2" ,"be/be150.8.6" ,"be/be200.8.10" ,"be/be250.4" ,"be/be120.3.9" ,"be/be150.3.3" ,"be/be150.8.7" ,
"be/be200.8.1" ,"be/be250.5" ,"be/be120.8.10" ,"be/be150.3.4" ,"be/be150.8.8" ,"be/be200.8.2" ,"be/be250.6" ,"be/be120.8.1" ,"be/be150.3.5" ,"be/be150.8.9" ,
"be/be200.8.3" ,"be/be250.7" ,"be/be120.8.2" ,"be/be150.3.6" ,"be/be200.3.10" ,"be/be200.8.4" ,"be/be250.8" ,"be/be120.8.3" ,"be/be150.3.7" ,"be/be200.3.1" ,
"be/be200.8.5" ,"be/be250.9"]

g = open("tempo_execucao_(x-y)_ipopt_be.txt", "w+")
for i=1:length(nomes_instancias)
for k=1:200
tempo = @time @CPUtime qubo( nomes_instancias[i])

println(g, nomes_instancias[i], " ", tempo )
end
end
close(g)
end

main()
