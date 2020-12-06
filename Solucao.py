#!python3

from matplotlib import path as pltpath, patches as pltpatches, pyplot as plt
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
RESULT_IMG_PATH = "resultado.png"



# Definição de funções usadas no código.
def minIndex(arr):
	### Encontra o índice do valor mínimo em um vetor. ###
	if len(arr) < 0:
		return -1
	minI = 0
	for i in range(1, len(arr)):
		if arr[i] < arr[minI]:
			minI = i
	return minI

def gerarImagem():
	### Gera uma imagem da solução encontrada, salva em disco, e exibe na tela (se suportado). ###
	global Z, L, path
	path_pos = [ L[path[i]] for i in range(len(path)) ]
	codes = [ pltpath.Path.MOVETO ] + [ pltpath.Path.LINETO for _ in range(1, len(path))]
	p = pltpath.Path(path_pos, codes)
	limx = [L[0][0], L[0][0], 0]
	limy = [L[0][1], L[0][1], 0]
	for pos in L:
		if pos[0] < limx[0]:
			limx[0] = pos[0]
		if pos[0] > limx[1]:
			limx[1] = pos[0]
		if pos[1] < limy[0]:
			limy[0] = pos[1]
		if pos[1] > limy[1]:
			limy[1] = pos[1]
	limx[2] = 0.02 * abs(limx[1] - limx[0])
	limy[2] = 0.02 * abs(limy[1] - limy[0])
	fig, ax = plt.subplots()
	patch = pltpatches.PathPatch(p, facecolor='black', fill=False, lw=2, alpha=0.75)
	ax.add_patch(patch)
	ax.plot([ p[0] for p in L ], [ p[1] for p in L ], 'bo')
	ax.set_xlim(limx[0] - limx[2], limx[1] + limx[2])
	ax.set_ylim(limy[0] - limy[2], limy[1] + limy[2])
	plt.savefig(RESULT_IMG_PATH)
	plt.show()



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
values = [ ]
for i in range(len(L)):
	values.append([ ])
	for j in range(len(L)):
		values[i].append(0)
salesman1 = 0
salesman2 = 0
available = [ i for i in range(len(L)) if i != salesman1 ]
while len(available) > 0:
	next_salesman1 = minIndex([ C[salesman1][i] for i in available ])
	values[salesman1][available[next_salesman1]] = 1
	salesman1 = available[next_salesman1]
	available.pop(next_salesman1)
	if len(available) < 1:
		break
	next_salesman2 = minIndex([ C[i][salesman2] for i in available ])
	values[available[next_salesman2]][salesman2] = 1
	salesman2 = available[next_salesman2]
	available.pop(next_salesman2)
values[salesman1][salesman2] = 1
variables = [ Z[i][j] for j in range(len(L)) for i in range(len(L)) ]
values = [ values[i][j] for j in range(len(L)) for i in range(len(L)) ]
solver.SetHint(variables, values) # Definir a solução inicial encontrada.

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
	execution_time = int(execution_time)
	execution_time = '%dmin %ds' % (int(execution_time/60), execution_time%60)



# Verificar se há solução e imprimir resultados.
if status != pywraplp.Solver.INFEASIBLE and status != pywraplp.Solver.NOT_SOLVED:

	# Gerar resultados.
	for i in range(len(L)):
		for j in range(len(L)):
			Z[i][j] = Z[i][j].solution_value()
	path = [ 0 ]
	current = Z[0].index(1)
	while current != 0:
		path.append(current)
		current = Z[current].index(1)
	path.append(0)

	# Exibir no terminal.
	print()
	if status == pywraplp.Solver.OPTIMAL:
		print('== Solução Ótima ==')
	else:
		print('== Solução Viável ==')
	print('\tTempo de execução: %s' % (execution_time))
	print('\tValor da função objetivo: %.3f' % (solver.Objective().Value()))
	print('\tInterações do Simplex: %d' % (solver.iterations()))
	print('\tNós explorados: %d' % (solver.nodes()))
	print('\tCaminho resultado: L1', end = '')
	for i in range(1, len(path)):
		print(' -> L%d' % (path[i] + 1), end='')
	print()
	gerarImagem()
	print('\t(uma imagem desse caminho foi salva em "%s")' % (RESULT_IMG_PATH))
	print()

else:
	print()
	print('== Sem Solução ==')
	print('\tTempo de execução: %s' % (execution_time))
	print('\tO problema não possui uma solução viável encontrada em tempo prático. :(')
	print()


