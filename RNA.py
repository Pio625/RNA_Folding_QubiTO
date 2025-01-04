import dwave
import dimod
import qubovert as qv
import qubovert.utils
import qubovert.utils 
from qubovert import boolean_var
from neal import SimulatedAnnealingSampler
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#here you should put the RNA SEQUENCE
RNA="UAGCAGCGGGAACAGUUCUGCAG"


#INIZIAMO A FARE IL PRE-PROCESSING PER LA CREAZIONE DELLE VARIABILI DA OTTIMIZZARE
#USIAMO COME VARIABILE DI BASE LE QUARTETTE E TOGLIAMO SUBITO LE COMBINAZIONI CHE NON POSSONO ESISTERE POICHÉ NON SI POSSONO LEGARE

#Faccio il  controllo mettendo dentro al x solo le coppie che hanno complementari

x={}

def check(RNA,i,j):
    valid_coupling={"AU","GC","UG"}
    if (RNA[i]+RNA[j] in valid_coupling) or (RNA[j]+RNA[i] in valid_coupling) :
        return True

Possible_x=[]
for i in range(len(RNA)):
    for j in range(len(RNA)):
        if i!=j and j-1>=0 and i+1<len(RNA):
            if check(RNA,i,j) and check(RNA,i+1,j-1):
                if  (i==j-1 and j==i+1 ) :
                    pass 
                #or ([j-1,i+1,j,i] in Possible_x)
                    
                else:
                    Possible_x.append([i,j,i+1,j-1]) 
                         
# da qui abbiamo i nostri bei quartetti filtrati

for i in range(len(Possible_x)):
    a=str(Possible_x[i][0])+"_"+str(Possible_x[i][1])+"_"+str(Possible_x[i][2]) +"_"+ str(Possible_x[i][3])
    x.update( {i:boolean_var('x(%s)' %  a)  } )

model=0
print(len(x))




#Dobbiamo iniziare a formare la cost function


#Trovo i valori di e_


e={"AUAU":0.9, "AUCG":2.1, "AUGC":2.4, "AUUA":1.1, "AUGU":0.6,"AUUG":1.4,
   "CGAU":2.1, "CGCG":3.3, "CGGC":2.4, "CGUA":2.1, "CGGU":1.4,"CGUG":2.1,
    "GCAU":2.4, "GCCG":3.4, "GCGC":3.3, "GCUA":2.2, "GCGU":1.5,"GCUG":2.5,
    "UAAU":1.3, "UACG":2.4, "UAGC":2.1, "UAUA":0.9, "UAGU":1,"UAUG":1.3,
    "GUAU":1.3, "GUCG":2.5, "GUGC":2.1, "GUUA":1.4, "GUGU":0.5,"GUUG":-1.3,
    "UGAU":1, "UGCG":1.5, "UGGC":1.4, "UGUA":0.6, "UGGU":-0.3,"UGUG":0.5 }
 

#primo termine ovvero quello sulle quartette

for i in range(len(Possible_x)):
    e_i=e[RNA[Possible_x[i][0]] + RNA[Possible_x[i][1]]+RNA[Possible_x[i][2]]+RNA[Possible_x[i][3]]]
    model = model +e_i*x[i]

#Ora passiamo invece alla creazione del secondo termine della obj function


r=-1.1
p=-3
t=2

#il quartetto che viene dopo è nella forma i+n,j-n,i+n+1,j-n-1
#vado anche a disincentivare la formazioni di quartette con finale AU e GU


for i in range(len(Possible_x)):
    cont=1
    for j in range(len(Possible_x)):
        if Possible_x[i][0]+cont<len(RNA)-1 and i!=j and Possible_x[i][1]-cont>0 and Possible_x[j]==[Possible_x[i][0]+cont,Possible_x[i][1]-cont,Possible_x[i][0]+cont+1,Possible_x[i][1]-cont-1]:
            model= model + r*x[i]*x[j]
            if (RNA[Possible_x[j][2]]+RNA[Possible_x[j][3]])=="GU" or (RNA[Possible_x[j][2]]+RNA[Possible_x[j][3]])=="UA" :
                model = model +p*x[i]*(1-x[j])
            cont=cont + 1
        
#for i in range(len(Possible_x)):
#    if RNA[Possible_x[i][2]]+RNA[Possible_x[i][3]]=="AU":
#        print(RNA[Possible_x[i][0]],RNA[Possible_x[i][1]],RNA[Possible_x[i][2]],RNA[Possible_x[i][3]],i)

#Cerco di andare a formare dei crossing nodes

def is_crossing_quartet(q1, q2):
    # Unpack the quartets
    i, j, i_plus1, j_minus1 = q1
    i_prime, j_prime, i_prime_plus1, j_prime_minus1 = q2

    # Check if any pair crosses
    # i < i' < j < j'
    condition_1 = i < i_prime and i_prime < j and j < j_prime
    # i < i'+1 < j < j'-1
    condition_2 = i < i_prime_plus1 and i_prime_plus1 < j and j < j_prime_minus1
    # i+1 < i' < j-1 < j'
    condition_3 = i_plus1 < i_prime and i_prime < j_minus1 and j_minus1 < j_prime
    # i+1 < i'+1 < j-1 < j'-1
    condition_4 = i_plus1 < i_prime_plus1 and i_prime_plus1 < j_minus1 and j_minus1 < j_prime_minus1

    # If any of the conditions hold, return True (crossing exists)
    return condition_1 or condition_2 or condition_3 or condition_4



for i in range(len(Possible_x)):
    for j in range(len(Possible_x)):
        if i!=j and is_crossing_quartet(Possible_x[i],Possible_x[j]):
            model =model + t*x[i]*x[j]
    


#modello QUBO finito, ora si deve passare alla parte cicciotta, ovvero runnare su d-wave (Qui lo faccio col simulatore per vedere se dafunzia)

qubo=model.to_qubo()
dwave_qubo = qubo.Q

# solve with D-Wave
res = SimulatedAnnealingSampler().sample_qubo(dwave_qubo, num_reads=10)
qubo_solution = res.first.sample

# convert the qubo solution back to the solution to the model
model_solution = model.convert_solution(qubo_solution)
print(model.value(model_solution))

Strin=""
solution=[]
for i in model_solution:
    if model_solution[i]==1:
        Strin=i[2:len(i)-1].split("_")
        solution.append(list(map(int,Strin)))
            

#Questa parte finale serve per disegnare il grafico del RNA



#Graficare il problema 
from collections import Counter
# Funzione per calcolare tutte le tuple (a, b) e (c, d)
def find_repeated_tuples(array):
    

    # Estrazione delle tuple (ordinate per evitare duplicati come (a, b) e (b, a))
    tuples = []
    for row in array:
        tuples.append(tuple(sorted(row[:2])))  # (a, b)
        tuples.append(tuple(sorted(row[2:])))  # (c, d)

    # Contare la frequenza delle tuple
    tuple_counts = Counter(tuples)

    # Raccogliere solo le tuple che si ripetono almeno una volta
    repeated_tuples = [t for t, count in tuple_counts.items() if count > 1]

    return repeated_tuples

# Calcolo delle tuple ripetute
result = find_repeated_tuples(solution)

#uso il modulo networkx

G=nx.Graph()

color_map=[]
for i in RNA:
    if i=="A":
        color_map.append('#EBD687')
    elif i=="C":
        color_map.append('violet')
    elif i=="G":
        color_map.append('#87CEEB')
    elif i=="U":
        color_map.append('#EBA487')
      
dictRNA={}        
for i in range(len(RNA)):
    #aggiungo i nodi di partenza 
    G.add_node(i, name=RNA[i])   
    #aggiungo i colori
    dictRNA.update({i:RNA[i]}) 
    #aggiungo i collegamenti banali
    if i+1<len(RNA):
        G.add_edge(i,i+1)

for i in result:
    G.add_edge(i[0],i[1])

nx.draw_kamada_kawai(G,node_color=color_map,labels=dictRNA,with_labels=True)
plt.show()
