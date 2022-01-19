# 2D FRAME ANALYSIS BY USING STIFFNESS METHOD

import numpy as np
import numpy.linalg as lin

n = 4 # Node  Numbers
m = 3 # Member Numbers
XY_values = [0.0, 0.0, 0.0, 3.0, 7.0, 3.0, 7.0, 0.0]
NC_values = [1, 2, 2, 3, 3, 4]
BC_values = [1,1,1,1, 4,1,1,1]
FEXT_values = [2, 10.0, 0, 0]

#calculating member end forces
def calculate_memberEndForce_q(k_prime,T,U):

    matrix = np.dot(k_prime,np.dot(T,U)) 
    return matrix

def calculate_R(Kpf,Df,Kpp,Dp):

    matrix = np.dot(Kpf,Df) + np.dot(Kpp,Dp)

    return matrix

def calculate_Df (Kff,F,Kfp,Dp):
    division_matrix = np.subtract(F,(Kfp.dot(Dp)))
    matrix = lin.solve(Kff,division_matrix)
    return matrix  
    

def assemble_matrix_1_to_2_using_Ce(Ce,k,K):

    size_array_k=k.shape

    row_k=size_array_k[0]

    col_k=size_array_k[1]


    # loop over k matrix's elements 

    for i in range(0,row_k,1):

        for j in range(0,col_k,1):

            # defining indices of K value by using DOF_numbers

            K_index_1 = int(Ce[j]-1.0) 

            K_index_2 = int(Ce[i]-1.0) 

        

            K_value = K[K_index_1,K_index_2] 

            

           # adding previous K value to new one recursively 

            K[K_index_1][K_index_2] = K_value + k[i][j]  

    return K

# defining el_stiff function to create element stiffness matrix
def el_stiff(T,k_prime):
    T_transpose = np.transpose(T) 
    matrix = np.dot(T_transpose,np.dot(k_prime,T))
    return matrix

def T_matrix(x,y):

    delta = find_delta_Values(x,y)

    delta_x = delta[0]

    delta_y = delta[1]


    L = findLength(x,y)

    cos_theta = delta_x / L

    sin_theta = delta_y / L 

    negative_sin_theta =(-1*delta_y/L) 

    matrix = [[cos_theta,sin_theta,0,0,0,0],
              [negative_sin_theta,cos_theta,0,0,0,0],
              [0,0,1,0,0,0],
              [0,0,0,cos_theta,sin_theta,0],
              [0,0,0,negative_sin_theta,cos_theta,0],
              [0,0,0,0,0,1]]

    return matrix


def el_stiff_local(A,E,I,x,y):
    L = findLength(x,y)
    
     
    matrix=[[(E*A/L),0,0,(-(E*A)/L),0,0],
            
            [0,(12*I*E/L**3),(6*I*E)/(L**2),0,(-12*I*E/L**3),(6*I*E/L**2)],
            [0,(6*I*E/L**2),(4*I*E/L),0,(-6*I*E/L**2),(2*I*E/L)],
            [(-E*A/L),0,0,(E*A/L),0,0],
            [0,(-12*E*I/L**3),(-6*E*I/L**2),0,(12*E*I/L**3),(-6*E*I/L**2)],
            [0,(6*E*I/L**2),(2*E*I/L),0,(-6*E*I/L**2),(4*E*I/L)]]

    return matrix

def findLength(x,y):
    delta = find_delta_Values(x,y)
    delta_x=delta[0]
    delta_y=delta[1]
    L_value=np.sqrt((delta_x**2)+(delta_y**2))
    return L_value

def find_delta_Values(x,y):
   start_x= x[0]
   end_x= x[1]
   
   start_y= y[0]
   end_y= y[1]
   
   delta_x=(end_x[0]-start_x[0])
   delta_y=(end_y[0]-start_y[0])
   delta=[delta_x,delta_y]
   return delta


def createMatrix(row,column,values):

    matrix=np.zeros([row, column], dtype =float)

    counter = 0 # counter for element number in values

    for i in range(0,row):

         for j in range(0,column):

             matrix[i][j] = values[counter]

             counter = counter + 1

    return matrix


XY = createMatrix(n,2, XY_values) 
NC = createMatrix(m,2, NC_values)
BC = createMatrix(2,4, BC_values)
F = createMatrix(1,4,FEXT_values)


#DOF----------------------------------------------------------------------------------------------------------------

DOF = np.zeros((n,3))

disp_counter= 0

row_BC = len(BC)
k = 1  # defining k as counter
for i in range(0,n,1):              # loop over nodes
    BC_state="not_exist" 
    for j in range(0,row_BC,1):     # loop over rows of BC array
        if BC[j,0]==i+1:            # Check whether there exists BC on the node
            BC_state = "exist" 
            x_value = j               # defining x_value as in which row BC on the node exists          
           
        
    
    if BC_state=="not_exist":       # if there is no BC on the node
        DOF[i,0]=k
        DOF[i,1]=k+1
        DOF[i,2]=k+2
        k=k+3     
        disp_counter=disp_counter+3
        
    elif BC_state == "exist" and BC[x_value,1]==0:     # elif there is BC on the node AND X dir. is NOT constrained
        DOF[i,0]=k 
        k=k+1 
        disp_counter= disp_counter+1    
    elif BC_state == "exist" and BC[x_value,2]==0:   # elif there is BC on the node AND Y dir. is NOT constrained
        DOF[i,1]=k 
        k=k+1 
        disp_counter= disp_counter+1
    if  BC_state == "exist" and BC[x_value,3]==0:     # elif there is BC on the node AND Z dir. is NOT constrained
        DOF[i,2]=k 
        k=k+1 
        disp_counter= disp_counter+1
        
print(disp_counter)
# filling the remaining 0 entries of the DOF matrix

for i in range(0,n,1):         # loop over nodes
    for j in range(0,3,1):     # Loop over X,Y and Z directions
        if DOF[i][j]==0:
            DOF[i][j]=k
            k=k+1
print(DOF)   



MG=np.zeros((m,3), dtype=float)

# C30
E = 32e5 
#Section Properties
hc = 0.7 #Column Depth
bc = 0.7 #Column Width
hb = 0.5 #Beam Depth
bb = 0.35 #Beam Width

#Moment Of Inertia
Ic = (bc*(hc**3))/12
Ib = (bb*(hb**3))/12

for i in range(0,len(NC),1):
    
    bn_id = int(NC[i,0]-1.0)
    
    en_id = int(NC[i,1]-1.0)
    
    x = np.array([XY[bn_id,0],XY[en_id,0]])
    
    y = np.array([XY[bn_id,1],XY[en_id,1]])
    
    if x[0] == x[1] and y[0]!=y[1]:
    
        MG[i,0] = E
        
        A_col = bc*hc #m^2
        
        MG[i,1] = A_col
        
        I_col = (1/12)*bc*(hc**3) #m^4
        
        MG[i,2] = I_col
        
    else:
        MG[i,0] = E
        
        A_beam = bb*hb #m^2
        
        MG[i,1] = A_beam
        
        I_beam = (1/12)*bb*(hb**3) #m^4
        
        MG[i,2] = I_beam
        

size_array_DOF = DOF.shape

row_DOF=size_array_DOF[0]

col_DOF=size_array_DOF[1]


q=row_DOF*col_DOF


size_array_NC=NC.shape
row_NC=size_array_NC[0]
col_NC=size_array_NC[1]

k_primee = []
T_valuess = []
K=np.zeros((q,q))
#assembling K matrix , get global Stifness matrix.
for e in range(len(NC)): #bn_id başlangıç id
    bn_id=int(NC[e,0]-1)
    
    en_id=int(NC[e,1]-1)
    x=np.array([[XY[bn_id,0]], [XY[en_id,0]]])
    y=np.array([[XY[bn_id,1]], [XY[en_id,1]]])
    E=MG[e,0]
    A=MG[e,1]
    I=MG[e,2]
    
    k_prime = el_stiff_local(A,E,I,x,y)
    T = T_matrix(x,y)                              
    T_valuess.append(T)
    k_primee.append(k_prime)
    k=el_stiff(T,k_prime)
    
    Ce=[DOF[bn_id,0],DOF[bn_id,1],DOF[bn_id,2] ,DOF[en_id,0] ,DOF[en_id,1],DOF[en_id,2]] 
    
  
    K=assemble_matrix_1_to_2_using_Ce(Ce,k,K)    
#print(K)  
np.set_printoptions(linewidth=np.inf)

# creating all zero Q_values array with q elements

Q_values=np.empty([q], dtype=float) 



for i in range(0,q,1):

    Q_values[i] = 0

Q=createMatrix(q,1,Q_values)


size_array_F=F.shape

row_F=size_array_F[0]

col_F=size_array_F[1]

for i in range(0,row_F,1):

    joint_id = int(F[i,0]-1)  # Joint number on which an external force is applied

    Cn = DOF[joint_id,:] # dof numbers of the node
    
    first_index = int(Cn[0]-1) 

    second_index = int(Cn[1]-1)

    third_index = int(Cn[2]-1)
    
    Q[first_index] = Q[first_index] + F[i,1]    # assembly of Fx to Q

    Q[second_index] = Q[second_index] + F[i,2]  # assembly of Fy to Q

    Q[third_index] = Q[third_index] + F[i,3]    # assembly of M to Q
    
    

rf = disp_counter  # rf: number of rows of Kff,displacement counter 

Kff = K[0:rf,0:rf]    

Kfp = K[0:rf,rf:]

Kpf = K[rf:,0:rf]

Kpp = K[rf:,rf:]

size_array_Cn=Cn.shape
row_Cn=1
col_Cn=size_array_Cn[0]   

F = Q[0:rf]

print(F)
# PRESCRIBED DISPLACEMENTS : ALL OF THEM ARE ZERO
Dp_values=np.empty([q-rf], dtype=float)
for i in range(0,q-rf,1): 
    Dp_values[i] = 0
Dp=createMatrix(q-rf,1,Dp_values) 



Df = calculate_Df(Kff,F,Kfp,Dp)


len_Df=len(Df)
len_Dp = len(Dp)
D=np.empty([(len_Df+len_Dp),1], dtype=float)
counter_D_matrix=0

for i in range(0,len_Df,1):
    D[counter_D_matrix]=Df[i]
    counter_D_matrix = counter_D_matrix+1
for j in range(0,len_Dp,1):
    D[counter_D_matrix]=Dp[j]
    counter_D_matrix=counter_D_matrix+1
print("D=")

print (D)


R=calculate_R(Kpf,Df,Kpp,Dp)
print ("R",R)

row_NC= len(NC)
for i in range(0,row_NC,1):
    D_values=np.zeros(row_NC*col_NC, dtype=float)
    counter_D=0
    
    for j in range(0,col_NC,1):
        row_DOF=int(NC[i,j])
        
        D1=int(DOF[row_DOF-1,0])
        D_values[counter_D]=D[D1-1]
        counter_D=counter_D+1
        
        D2=int(DOF[row_DOF-1,1])
        D_values[counter_D]=D[D2-1]
        counter_D=counter_D+1
        
        D3=int(DOF[row_DOF-1,2])
        D_values[counter_D]=D[D3-1]
        counter_D=counter_D+1
    
    U = createMatrix(6,1,D_values)
    
    f_prime = calculate_memberEndForce_q(k_primee[i],T_valuess[i],U)
    
    print("fprime",(i+1),"=")
    print(f_prime)