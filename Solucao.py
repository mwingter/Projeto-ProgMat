#!python3

from matplotlib import path as pltpath, patches as pltpatches, pyplot as plt
from ortools.linear_solver import pywraplp
import numpy as np
import time
import math
import sys
import os
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



# Definições do programa.
MAX_EXECUTION_TIME = 60*10 # Tempo máximo de execução do solver (em segundos).
MAX_2OPT_EXECUTION_TIME = 60*1 # Tempo máximo de execução da heurística (em segundos). Recomendado não alterar.
RESULT_IMG_PATH = "solucao-encontrada.png" # Nome do arquivo em que a solução será salva como imagem.
INVERT_AXIS = False # Inverter eixos x/y da entrada (só muda efetivamente a imagem final gerada).
ROUND_POINTS = False # Realizar arredondamentos nos pontos ou não (True ou False).
ROUND_EUCLIDIAN = True # Realizar arredondamentos na distância euclidiana ou não (True ou False).
SOLVER_NAME = 'SCIP' # Use 'SCIP' ou 'CP-SAT'.



# Primeiro, vamos instanciar a classe do nosso solver (OR-Tools) e declarar constantes usadas ao longo do código.
solver = pywraplp.Solver.CreateSolver(SOLVER_NAME)
INFINITY = solver.infinity()
REGEX_3NUMBERS = re.compile(r'^\s*[0-9]+\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.
REGEX_2NUMBERS = re.compile(r'^\s*([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s+([+-]?[0-9]+(?:\.[0-9]+)?|[+-]?\.[0-9]+)\s*$') # Usado no código.

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
	_, ax = plt.subplots()
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	patch = pltpatches.PathPatch(p, facecolor='black', fill=False, lw=2, alpha=0.75)
	ax.add_patch(patch)
	ax.plot([ L[p][0] for p in path[1:len(path)-2] ], [ L[p][1] for p in path[1:len(path)-2] ], 'b.')
	ax.plot([ L[path[0]][0] ], [ L[path[0]][1] ], 'b*')
	ax.plot([ L[path[len(p) - 2]][0] ], [ L[path[len(p) - 2]][1] ], 'rx')
	ax.set_xlim(limx[0] - limx[2], limx[1] + limx[2])
	ax.set_ylim(limy[0] - limy[2], limy[1] + limy[2])
	plt.savefig(RESULT_IMG_PATH)
	plt.show()
def swap2opt(i, j):
	### Função suporte para o algoritmo da heurística 2-OPT. ###
	global path, C
	new_path = path.copy()
	if i < j:
		new_path[i:j + 1] = reversed(new_path[i:j + 1])
		if i == 0: # Caso em que troca o nó inicial.
			new_path[len(new_path) - 1] = new_path[0]
	cost_sum = 0
	for i in range(1, len(new_path)):
		cost_sum += C[new_path[i - 1]][new_path[i]]
	return new_path, cost_sum



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

# Inverter eixos x/y da entrada se necessário.
if INVERT_AXIS:
	for i in range(len(L)):
		L[i] = (L[i][1], L[i][0])

# Arredondar valores dos pontos se necessário.
if ROUND_POINTS:
	for i in range(len(L)):
		L[i] = (round(L[i][0]), round(L[i][1]))

# Calcular custos, ou seja, distâncias percorridas pelo telescópio (variável C).
C = [ ] # Armazenará os custos de locomoção do telescópio (da galáxia i para a galáxia j).
if ROUND_EUCLIDIAN:
	euclidian_distance = lambda p1, p2: round(math.sqrt(math.pow(p2[0] - p1[0], 2) + math.pow(p2[1] - p1[1], 2))) # Função de distância euclidiana arredondada.
else:
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
		values[i].append(0) # Preencher solução inicial com tudo como '0'.
salesman1 = 0 # Anda para frente.
salesman2 = 0 # Anda para trás.
available = [ i for i in range(len(L)) if i != salesman1 ] # Pontos disponíveis para andar.
while len(available) > 0:
	# Andar um ponto para frente (o mais próximo).
	next_salesman1 = minIndex([ C[salesman1][i] for i in available ])
	values[salesman1][available[next_salesman1]] = 1
	salesman1 = available[next_salesman1]
	available.pop(next_salesman1)
	if len(available) < 1:
		break
	# Andar um ponto para trás (o mais próximo).
	next_salesman2 = minIndex([ C[i][salesman2] for i in available ])
	values[available[next_salesman2]][salesman2] = 1
	salesman2 = available[next_salesman2]
	available.pop(next_salesman2)
values[salesman1][salesman2] = 1
variables = [ Z[i][j] for j in range(len(L)) for i in range(len(L)) ] # Matriz 2D -> Vetor 1D.
values = [ values[i][j] for j in range(len(L)) for i in range(len(L)) ] # Matriz 2D -> Vetor 1D.
solver.SetHint(variables, values) # Definir a solução inicial encontrada.

# Definir tempo máximo de execução.
solver.SetTimeLimit(1000 * MAX_EXECUTION_TIME) # Multiplica por 1000 porque é em milisegundos.

# Definir processamento paralelo, se o solver for o 'CP-SAT'.
if SOLVER_NAME == 'CP-SAT':
	cpu_count = os.cpu_count()
	if cpu_count is None:
		solver.SetNumThreads(2)
	else:
		solver.SetNumThreads(max(2, min(64, cpu_count)))



# Resolver utilizando o OR-Tools enquanto calcula o tempo de execução.
print("Buscando solução em tempo máximo de %dmin %ds..." % (MAX_EXECUTION_TIME/60, MAX_EXECUTION_TIME%60))
start_time = time.time()
status = solver.Solve()
end_time = time.time()
execution_time = int(end_time - start_time)



# Verificar se não encontrou nenhuma solução, e se for o caso encerrar.
if status == pywraplp.Solver.INFEASIBLE or status == pywraplp.Solver.NOT_SOLVED:
	print()
	print('== Sem Solução ==')
	print('\tO problema não possui uma solução viável encontrada em tempo prático. :(')
	print('\tTempo de execução: %dmin %ds' % (int(execution_time/60), execution_time%60))
	print()
	sys.exit(0)

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

# Melhorar com algoritmo 2-OPT.
print("Otimizando com heurística 2-OPT...")
current_objective_value = solver.Objective().Value()
old_objective_value = current_objective_value + 1
max_2opt_time = max(0, MAX_EXECUTION_TIME - execution_time) + MAX_2OPT_EXECUTION_TIME
start_2opt_time = time.time()
try:
	while old_objective_value > current_objective_value and time.time() - start_2opt_time < max_2opt_time:
		old_objective_value = current_objective_value
		for i in range(len(L) - 1):
			for j in range(i + 1, len(L) - 2):
				new_path, new_objective_value = swap2opt(i, j)
				if new_objective_value < current_objective_value:
					path = new_path
					current_objective_value = new_objective_value
					break
			else:
				continue
			break
	if old_objective_value > current_objective_value:
		print("Otimização 2-OPT esgotada.")
except KeyboardInterrupt:
	print("Otimização 2-OPT interrompida.")
	_, current_objective_value = swap2opt(0, 0)
end_2opt_time = time.time()
execution_2opt_time = int(end_2opt_time - start_2opt_time)

# Exibir no terminal.
print()
if status == pywraplp.Solver.OPTIMAL:
	print('== Solução Ótima ==')
else:
	print('== Solução Viável ==')
print('\tInterações do Simplex: %d' % (solver.iterations()))
print('\tNós explorados: %d' % (solver.nodes()))
print('\tValor da função objetivo: %.3f' % (solver.Objective().Value()))
print('\tValor da função objetivo após heurísticas: %.3f' % (current_objective_value))
print('\tTempo de execução do solver: %dmin %ds' % (int(execution_time/60), execution_time%60))
print('\tTempo de execução da heurística: %dmin %ds' % (int(execution_2opt_time/60), execution_2opt_time%60))
print('\tCaminho resultado: L%d' % (path[0] + 1), end = '')
for i in range(1, len(path)):
	print(' -> L%d' % (path[i] + 1), end='')
print()
gerarImagem()
print('\t(uma imagem desse caminho foi salva em "%s")' % (RESULT_IMG_PATH))
print()


