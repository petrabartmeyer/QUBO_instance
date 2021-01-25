using JuMP, LinearAlgebra, DelimitedFiles, Ipopt, CPUTime


function qubo(nome)

readQ = readdlm(nome)
BigM = 50000
#n = 3000
#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))


println(nome)
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
#model = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "OutputFlag" => 0))
#model = Model(Alpine.Optimizer)
count = 1


Q = zeros(readQ[1,1],readQ[1,1])
for i=2:size(readQ,1)
	Q[readQ[i,1],readQ[i,2]] = readQ[i,3]
	Q[readQ[i,2],readQ[i,1]] = readQ[i,3]
end

n= readQ[1,1]

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

f = open("resultados_(x-y)_bqp.txt", "a+")
println(f, nome, " ", n, " ",rho," ", c, " ",lambda, " ",epsolon_1," ", error_phi, " ", obj_value, " ", maximum(JuMP.value.(x)), " ", minimum(JuMP.value.(x)))
close(f)

return obj_value 
end

function main()
nomes_instancias =  ["beasley(bqp)/bqp1000-10.sparse.txt","beasley(bqp)/bqp100-8.sparse","beasley(bqp)/bqp250-3.sparse","beasley(bqp)/bqp500-8.sparse","beasley(bqp)/bqp1000-2.sparse.txt","beasley(bqp)/bqp100-9.sparse","beasley(bqp)/bqp250-4.sparse","beasley(bqp)/bqp500-9.sparse",
"beasley(bqp)/bqp1000-5.sparse.txt","beasley(bqp)/bqp2500-10.sparse.txt","beasley(bqp)/bqp250-5.sparse","beasley(bqp)/bqp50-10.sparse","beasley(bqp)/bqp1000-6.sparse.txt","beasley(bqp)/bqp2500-1.sparse.txt","beasley(bqp)/bqp250-6.sparse","beasley(bqp)/bqp50-1.sparse",
"beasley(bqp)/bqp1000-7.sparse.txt","beasley(bqp)/bqp2500-2.sparse.txt","beasley(bqp)/bqp250-7.sparse","beasley(bqp)/bqp50-2.sparse","beasley(bqp)/bqp1000-8.sparse.txt","beasley(bqp)/bqp2500-3.sparse.txt","beasley(bqp)/bqp250-8.sparse","beasley(bqp)/bqp50-3.sparse","beasley(bqp)/bqp1000-9.sparse.txt",
"beasley(bqp)/bqp2500-4.sparse.txt","beasley(bqp)/bqp250-9.sparse","beasley(bqp)/bqp50-4.sparse","beasley(bqp)/bqp100-10.sparse","beasley(bqp)/bqp2500-5.sparse.txt","beasley(bqp)/bqp500-10.sparse","beasley(bqp)/bqp50-5.sparse","beasley(bqp)/bqp100-1.sparse","beasley(bqp)/bqp2500-6.sparse.txt",
"beasley(bqp)/bqp500-1.sparse","beasley(bqp)/bqp50-6.sparse","beasley(bqp)/bqp100-2.sparse","beasley(bqp)/bqp2500-7.sparse.txt","beasley(bqp)/bqp500-2.sparse","beasley(bqp)/bqp50-7.sparse","beasley(bqp)/bqp100-3.sparse","beasley(bqp)/bqp2500-8.sparse.txt","beasley(bqp)/bqp500-3.sparse","beasley(bqp)/bqp50-8.sparse",
"beasley(bqp)/bqp100-4.sparse","beasley(bqp)/bqp2500-9.sparse.txt","beasley(bqp)/bqp500-4.sparse","beasley(bqp)/bqp50-9.sparse","beasley(bqp)/bqp100-5.sparse","beasley(bqp)/bqp250-10.sparse","beasley(bqp)/bqp500-5.sparse","beasley(bqp)/bqp100-6.sparse","beasley(bqp)/bqp250-1.sparse","beasley(bqp)/bqp500-6.sparse",
"beasley(bqp)/bqp100-7.sparse","beasley(bqp)/bqp250-2.sparse","beasley(bqp)/bqp500-7.sparse"]


g = open("tempo_execucao_(x-y)_ipopt_bqp.txt", "w+")
for i=1:length(nomes_instancias)
for k=1:200
tempo = @time @CPUtime qubo( nomes_instancias[i])

println(g, nomes_instancias[i], " ", tempo )
end
end
close(g)
end

main()
