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

#I have used for the project Qubovert and networkx for graph some problem
#for this project i have used a simulated annealer but is very easy to make it with the real one (the solving part use the D-wave "formalism")


#here you should put the RNA SEQUENCE
RNA="UCUCCGAUCUUCGGUGUCGAGU"


#Pre-processing
#We start creating the quartet that can be physically connected.
# Experimentaly we see that other than the complementary quartete there are also GU\UG couple that are very common 


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
                         
#Here we have our filter quartet

#We put the filtret quartet in the format that Qubovert want
for i in range(len(Possible_x)):
    a=str(Possible_x[i][0])+"_"+str(Possible_x[i][1])+"_"+str(Possible_x[i][2]) +"_"+ str(Possible_x[i][3])
    x.update( {i:boolean_var('x(%s)' %  a)  } )

model=0
print(len(x))




#Our quartet now is formed and set in the qubovert version


#The first part of our quartet is the quartet enerty part (infact every quartet contribute different energetically)
#So the first part of our const function is the quartet multiplied by e_i
#How you can see the direction of the quartet does mattet, infact 
#for example logically we could think thath AUUG is the same AUGU but this is not the case
#Infact the order in where the pairs come can change the energy of the structure

e={"AUAU":0.9, "AUCG":2.1, "AUGC":2.4, "AUUA":1.1, "AUGU":0.6,"AUUG":1.4,
   "CGAU":2.1, "CGCG":3.3, "CGGC":2.4, "CGUA":2.1, "CGGU":1.4,"CGUG":2.1,
    "GCAU":2.4, "GCCG":3.4, "GCGC":3.3, "GCUA":2.2, "GCGU":1.5,"GCUG":2.5,
    "UAAU":1.3, "UACG":2.4, "UAGC":2.1, "UAUA":0.9, "UAGU":1,"UAUG":1.3,
    "GUAU":1.3, "GUCG":2.5, "GUGC":2.1, "GUUA":1.4, "GUGU":0.5,"GUUG":-1.3,
    "UGAU":1, "UGCG":1.5, "UGGC":1.4, "UGUA":0.6, "UGGU":-0.3,"UGUG":0.5 }
 

#We start to create the obj function, infact we usa a for loop for the summatin of all the quartet

for i in range(len(Possible_x)):
    e_i=e[RNA[Possible_x[i][0]] + RNA[Possible_x[i][1]]+RNA[Possible_x[i][2]]+RNA[Possible_x[i][3]]]
    model = model +e_i*x[i]

#Now we start with the creation of the second and third part of our obj function

r=-1.1
p=-3
t=2

#Stacked quartet are in the form (i+n,j-n,i+n+1,j-n-1) and are very more stable than the single 
#in the obj function i try to discourage the formation of GU and UA end pairs


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

#here we see if there are crossed quartet

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
    


#Here we have build all the qubo model, now we have the part where we start to
#see the Hardware part of the system

#Vi convert in a Qubo mode

qubo=model.to_qubo()

#we convert from qubo model to the D-wave format
dwave_qubo = qubo.Q

# solve with D-Wave (we use a simulator but it is easy also to use the real hardware)
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
            

#The last part is used to draw the final structure of the system, 
#i use the Networkx model for be fast, but in the future i would like to make a rapresentation 
#that take in consideration also the angle of the bond, in the end it does no matter but for
#visualization purpose can be usefull



#Here we make this part just for the visualization part
#is my first time using this library and so in the end probabily exist a better way
#for the creation of the graph.


from collections import Counter
# function for extract the  tuple (a, b) e (c, d)
def find_repeated_tuples(array):
    tuples = []
    for row in array:
        tuples.append(tuple(sorted(row[:2])))  # (a, b)
        tuples.append(tuple(sorted(row[2:])))  # (c, d)

    # Tuple frequency
    tuple_counts = Counter(tuples)

    # take all the couple that are taken more that 2 time
    repeated_tuples = [t for t, count in tuple_counts.items() if count > 1]

    return repeated_tuples

result = find_repeated_tuples(solution)

#There are some optimization problem (of the program):
#The result of the solution are 2 times bigger the the one needed
#because as said before there are ordination problem {(a,b,c,d) is different from (c,d,a,b}
#there should be a way to rappresent the system that can make use of the symmetry
#but all my attempt where useless and had make the prediction power weaker

#Networkx module

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
    #add the node
    G.add_node(i, name=RNA[i])   
    #add the color
    dictRNA.update({i:RNA[i]}) 
    #i add the linear connection
    if i+1<len(RNA):
        G.add_edge(i,i+1)

for i in result:
    G.add_edge(i[0],i[1])

nx.draw_kamada_kawai(G,node_color=color_map,labels=dictRNA,with_labels=True)
plt.show()
