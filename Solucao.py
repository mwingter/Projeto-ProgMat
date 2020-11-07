from ortools.linear_solver import pywraplp



#
 #	~ Trabalho: Parte 1 ~
 #
 #	Julia (Nº USP XXXXX)
 #	Michelle (Nº USP XXXXX)
 #	Matheus Carvalho Raimundo (Nº USP 10369014)
 #	Marcelo (Nº USP XXXXX)
 #
 #	Programação Matemática: SME-0110 2020.2
 #
 #	 _______ _______ _______
 #	|   |   |               \
 #	|   |   |      \    |___|
 #	|   |    \      |   |
 #	|_______ _______|___|
 #
#



# Primeiro, vamos instanciar a classe do nosso solver (OR-Tools) e declarar constantes usadas ao longo do código.
solver = pywraplp.Solver.CreateSolver('SCIP')
INFINITY = solver.infinity()



# Inicialização do programa.
print('Bem-vinde ao nosso programa.')
print()



# Definir variáveis usando nossa modelagem.
x = solver.IntVar(0.0, INFINITY, 'x')
y = solver.IntVar(0.0, INFINITY, 'y')
print('Número de variáveis =', solver.NumVariables())



# Definir restrições usando nossa modelagem.
solver.Add(x + 7 * y <= 17.5) # x + 7 * y <= 17.5.
solver.Add(x <= 3.5) # x <= 3.5.
print('Número de restrições =', solver.NumConstraints())



# Definir a função objetivo usando nossa modelagem.
solver.Maximize(x + 10 * y) # Max x + 10 * y.



# Resolver utilizando o OR-Tools.
status = solver.Solve()



# Verificar se há solução e imprimir resultados.
if status == pywraplp.Solver.OPTIMAL:
	print()
	print('== Solução ==')
	print('Valor da função objetivo:', solver.Objective().Value())
	print('x:', x.solution_value())
	print('y:', y.solution_value())
	print()
else:
	print()
	print('== Solução ==')
	print('O problema não possui uma solução ótima.')
	print()
