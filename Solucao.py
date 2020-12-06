#!python3

from ortools.linear_solver import pywraplp
import multiprocessing
import numpy as np
import time
import math
import sys
import re



#
 #	~ Trabalho: Parte 1 ~
 #
 #	Julia Carolina Frare Peixoto (Nº USP 10734727)
 #	Michelle Wingter da Silva (Nº USP 10783243)
 #	Matheus Carvalho Raimundo (Nº USP 10369014)
 #	Marcelo Duchêne (Nº USP 8596351)
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
REGEX_3NUMBERS = re.compile(r'^\s*[0-9]+\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.
REGEX_2NUMBERS = re.compile(r'^\s*([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.
MAX_EXECUTION_TIME = 10 # Tempo máximo de execução (em segundos).



# Inicialização do programa.
print('Bem-vinde ao nosso programa.')
print('Por favor, entre com as localizações de cada galáxia, uma por linha (formato: posição_x posição_y). Ao fim, entre com EOF (CTRL+D).')



# Fazer leitura da entrada (variável L).
L = [ ] # Armazenará as posições x e y das galáxias.
while(True):
	try:
		line = input("> ").strip()
		if line == "":
			continue
		grps = REGEX_3NUMBERS.findall(line)
		if len(grps) == 1:
			L.append(( float(grps[0][0]), float(grps[0][1]), ))
			continue
		grps = REGEX_2NUMBERS.findall(line)
		if len(grps) == 1:
			L.append(( float(grps[0][0]), float(grps[0][1]), ))
			continue
		print("Linha ignorada por ser incompatível. Por favor, use o seguinte formato: posição_x posição_y")
	except EOFError:
		break
print()
print('Número de galáxias =', len(L))

if len(L) < 1:
	print('== Solução ==')
	print('\tVocê não forneceu nenhuma galáxia como entrada.')
	print()
	sys.exit(0)
elif len(L) < 4:
	print('== Solução ==')
	print('\tVocê forneceu apenas %d galáxia/s como entrada. A solução é trivial.' % (len(L)))
	print()
	sys.exit(0)

# Calcular custos, ou seja, distâncias percorridas pelo telescópio (variável C).
C = [ ] # Armazenará os custos de locomoção do telescópio (da galáxia i para a galáxia j).
euclidian_distance = lambda p1, p2: math.sqrt(math.pow(p2[0] - p1[0], 2) + math.pow(p2[1] - p1[1], 2)) # Função de distância euclidiana.
for i in range(len(L)):
	C.append([ ])
	for j in range(len(L)):
		C[i].append( euclidian_distance(L[i], L[j]) )



# Definir as variáveis usando nossa modelagem.
Z = [ ] # Armazenará a variável Z (matriz).
for i in range(len(L)):
	Z.append([ ])
	for j in range(len(L)):
		Z[i].append( solver.IntVar(0.0, 1.0, 'Z_' + str(i) + '_' + str(j)) ) # Zij.
U = [ ] # Armazenará a variável extra u (vetor).
for i in range(1, len(L)):
	U.append( solver.IntVar(1.0, len(L), 'U_' + str(i)) )
print('Número de variáveis =', solver.NumVariables())



# Definir restrições usando nossa modelagem.
for i in range(len(L)):
	solver.Add(Z[i][i] == 0) # Zii = 0.
	solver.Add(solver.Sum([Z[i][j] for j in range(len(L))]) == 1) # Somatório de j = 1 até N de (Zij) = 1, para todo i.
	solver.Add(solver.Sum([Z[j][i] for j in range(len(L))]) == 1) # Somatório de j = 1 até N de (Zji) = 1, para todo i.
for i in range(1, len(L)):
	for j in range(1, len(L)):
		if i == j:
			continue
		solver.Add(U[i - 1] - U[j - 1] + len(L) * Z[i][j] <= len(L) - 1)
print('Número de restrições =', solver.NumConstraints())



# Definir a função objetivo usando nossa modelagem.
solver.Minimize(solver.Sum([np.dot(z_row, c_row) for (z_row, c_row) in zip(Z, C)])) # Min somatório de i,j em L de (Cij * Zij)



# Definir uma solução inicial para o solver.
variables = [ ]
values = [ ]
for i in range(len(L)):
	minValue = min([ C[i][j] for j in range(len(C[i])) if i != j ])
	minIndex = C[i].index(minValue)
	for j in range(len(Z[i])):
		variables.append(Z[i][j])
		if(j == minIndex):
			values.append(1)
		else:
			values.append(0)
solver.SetHint(variables, values) # Definir solução inicial.

# Definir tempo máximo de execução.
solver.SetTimeLimit(1000 * MAX_EXECUTION_TIME)



# Resolver utilizando o OR-Tools enquanto calcula o tempo de execução.
start_time = time.time()
status = solver.Solve()
end_time = time.time()

# Ajustar tempo de execução para texto legível.
execution_time = end_time - start_time
if execution_time < 60:
	execution_time = '%.3fs' % (execution_time)
else:
	execution_time = '%dmin %.3fs' % (int(execution_time/60), execution_time%60.0)



# Verificar se há solução e imprimir resultados.
if status != pywraplp.Solver.INFEASIBLE and status != pywraplp.Solver.NOT_SOLVED:
	for i in range(len(L)):
		for j in range(len(L)):
			Z[i][j] = Z[i][j].solution_value()
	print()
	if status == pywraplp.Solver.OPTIMAL:
		print('== Solução Ótima ==')
	else:
		print('== Solução Viável ==')
	print('\tTempo de execução: %s' % (execution_time))
	print('\tValor da função objetivo: %.3f' % (solver.Objective().Value()))
	for i in range(len(L)):
		print('\tDa galáxia %d, vamos para a galáxia %d.' % (i+1, Z[i].index(1) + 1))
	print('\tCaminho: 1', end = '')
	current = Z[0].index(1)
	while current != 0:
		print(' -> %d' % (current + 1), end='')
		current = Z[current].index(1)
	print(' -> 1')
	print('\t(uma imagem desse caminho foi salva em "%s")' % ('em construção...'))
	print()
else:
	print()
	print('== Sem Solução ==')
	print('\tTempo de execução: %s' % (execution_time))
	print('\tO problema não possui uma solução viável encontrada em tempo prático. :(')
	print()


