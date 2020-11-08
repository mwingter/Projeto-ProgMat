#!python3

from ortools.linear_solver import pywraplp
import itertools
import numpy as np
import math
import re



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
REGEX_3NUMBERS = re.compile(r'^\s*[0-9]+\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.
REGEX_2NUMBERS = re.compile(r'^\s*([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.



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
		Z[i].append( solver.IntVar(0.0, INFINITY, 'Z_' + str(i) + '_' + str(j)) ) # Zij.
print('Número de variáveis =', solver.NumVariables())



# Definir restrições usando nossa modelagem.
for i in range(len(L)):
	solver.Add(Z[i][i] == 0) # Zii = 0.
	solver.Add(solver.Sum([Z[i][j] for j in range(len(L))]) == 1) # Somatório de j = 1 até N de (Zij) = 1, para todo i.
	solver.Add(solver.Sum([Z[j][i] for j in range(len(L))]) == 1) # Somatório de j = 1 até N de (Zji) = 1, para todo i.
print('Gerando restrições anti-ciclos... %.2f%%' % (0), end='')
for n in range(2, len(L)):
	print('\rGerando restrições anti-ciclos... %.2f%%' % (100 * (n - 2)/(len(L) - 2)), end='')
	for s in itertools.combinations(range(len(L)), n):
		sum_terms = [Z[x_row][y_row] for (x_row, y_row) in itertools.combinations_with_replacement(s, 2) if x_row != y_row]
		sum_terms += [Z[y_row][x_row] for (x_row, y_row) in itertools.combinations_with_replacement(s, 2) if x_row != y_row]
		solver.Add(solver.Sum(sum_terms) <= n - 1) # Somatório de i,j em Q de (Zji) <= |Q|, para todo Q subconjunto das galáxias.
print('\rGerando restrições anti-ciclos... %.2f%%' % (100))
print('Número de restrições =', solver.NumConstraints())



# Definir a função objetivo usando nossa modelagem.
solver.Minimize(solver.Sum([np.dot(x_row, cost_row) for (x_row, cost_row) in zip(Z, C)])) # Min somatório de i,j em L de (Cij * Zij)



# Resolver utilizando o OR-Tools.
status = solver.Solve()



# Verificar se há solução e imprimir resultados.
if status == pywraplp.Solver.OPTIMAL:
	for i in range(len(L)):
		for j in range(len(L)):
			Z[i][j] = Z[i][j].solution_value()
	print()
	print('== Solução ==')
	print('\tValor da função objetivo: %.3f' % (solver.Objective().Value()))
	for i in range(len(L)):
		print('\tDa galáxia %d, vamos para a galáxia %d.' % (i+1, Z[i].index(1) + 1))
	print('Caminho: 1', end = '')
	current = Z[0].index(1)
	while current != 0:
		print(' -> %d' % (current + 1), end='')
		current = Z[current].index(1)
	print(' -> 1')
	print()
else:
	print()
	print('== Solução ==')
	print('O problema não possui uma solução ótima. :(')
	print()


