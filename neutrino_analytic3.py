###+++++++++ python rutine to get the anlytical neutrino eigenvalues +++++++++++++++++++++++++++++

import numpy as np

def MATRIXCHIDIAG(vevSM, MN, MPsi, Meta, l1, l2, y1, y2):    
    
    dcOut={}
    
    #Diagonalization of Mchi matrix by the bi-unitary transfortion V and U
    
    MX0 = np.matrix( [[MN,0,1./np.sqrt(2.)*vevSM*y2],
                      [0,MPsi,l1],
                      [1./np.sqrt(2.)*vevSM*y1,l2,Meta]])
    
    #squared eigenvalues e eigenvectors for the V MATRIX
    (MVdiag2,V)=np.linalg.eig(MX0*np.transpose(MX0))
    #squared eigenvalues e eigenvectors for the U MATRIX
    (MUdiag2,U)=np.linalg.eig(np.transpose(MX0)*MX0)


    M1=np.sqrt(np.abs(MVdiag2[0]))
    M2=np.sqrt(np.abs(MVdiag2[1]))
    M3=np.sqrt(np.abs(MVdiag2[2]))
    V11=V[0,0]
    V12=V[0,1]
    V13=V[0,2]
    V21=V[1,0]
    V22=V[1,1]
    V23=V[1,2]
    V31=V[2,0]
    V32=V[2,1]
    V33=V[2,2] 
    U11=U[0,0]
    U12=U[0,1]
    U13=U[0,2]
    U21=U[1,0]
    U22=U[1,1]
    U23=U[1,2]
    U31=U[2,0]
    U32=U[2,1]
    U33=U[2,2]   

    dcOut['M1']= M1
    dcOut['M2']= M2  
    dcOut['M3']= M3  
    dcOut['V11']= V11 
    dcOut['V12']= V12 
    dcOut['V13']= V13 
    dcOut['V21']= V21 
    dcOut['V22']= V22 
    dcOut['V23']= V23 
    dcOut['V31']= V31 
    dcOut['V32']= V32 
    dcOut['V33']= V33 
    dcOut['U11']= U11 
    dcOut['U12']= U12 
    dcOut['U13']= U13 
    dcOut['U21']= U21 
    dcOut['U22']= U22 
    dcOut['U23']= U23 
    dcOut['U31']= U31 
    dcOut['U32']= U32 
    dcOut['U33']= U33 
            
    return dcOut

#Loop factor
def LAMBDA(mNk,mSk,Vk2,Uk1):
    
    mk = 1./(16.*np.pi**2)*Vk2*Uk1*(mNk**3/(mNk**2-mSk**2))*np.log(mNk**2/mSk**2)
    
    return mk   

#Mab matrix. sum over i and k is expanded
def Mab(YB1b,YB2b,YA1a,YA2a,m1,m2,ms1,ms2,V12,V22,U11,U21):
    sumS1= (LAMBDA(m1, ms1, V12, U11)+LAMBDA(m2, ms1, V22, U21))*(YB1b*YA1a) 
    
    sumS2= (LAMBDA(m1, ms2, V12, U11)+LAMBDA(m2, ms2, V22, U21))*(YB2b*YA2a) 
    
    return sumS1 + sumS2

#Compute the neutrino eigenvalues
def MATRIX_NU_DIAG(YB11,YB12,YB13,YB21,YB22,YB23,YA11,YA12,YA13,YA21,YA22,YA23,mS1,mS2,m1,m2,m3,V12,V22,U11,U21):    
    
    #Matrix elements
    M11 = Mab(YB11,YB21,YA11,YA21,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M12 = Mab(YB11,YB21,YA12,YA22,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M13 = Mab(YB11,YB21,YA13,YA23,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M21 = Mab(YB12,YB22,YA11,YA21,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M22 = Mab(YB12,YB22,YA12,YA22,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M23 = Mab(YB12,YB22,YA13,YA23,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M31 = Mab(YB13,YB23,YA11,YA21,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M32 = Mab(YB13,YB23,YA12,YA22,m1,m2,mS1,mS2,V12,V22,U11,U21)
    M33 = Mab(YB13,YB23,YA13,YA23,m1,m2,mS1,mS2,V12,V22,U11,U21)


    Mvij = np.matrix( [[M11, M12, M13],
                       [M21, M22, M23],
                       [M31, M32, M33]] )

    #eigenvalues e eigenvectors
    (Mdiag2,V)=np.linalg.eig(Mvij*np.transpose(Mvij))
    
    #took eigenvalues
    MX1 = np.sqrt(np.abs(Mdiag2[0]))
    MX2 = np.sqrt(np.abs(Mdiag2[1]))
    MX3 = np.sqrt(np.abs(Mdiag2[2]))
    
    ## reorganize the eigenvalues (neutrino masses)
    mn1 = 0.0
    mn2 = 0.0
    mn3 = 0.0

    if MX1 < MX2 and MX1 < MX3:
        mn1 = MX1
        #print "Hola1"

        if MX2 < MX3:
            mn2 = MX2
            mn3 = MX3
        else:
            mn2 = MX3
            mn3 = MX2  

    if MX2 < MX1 and MX2 < MX3:
        mn1 = MX2
        #print "Hola2" 

        if MX1 < MX3:
            mn2 = MX1
            mn3 = MX3
        else:
            mn2 = MX3
            mn3 = MX1   

    if MX3 < MX1 and MX3 < MX2:
        mn1 = MX3
        #print "Hola3"  

        if MX1 < MX2:
            mn2 = MX1
            mn3 = MX2
        else:
            mn2 = MX2
            mn3 = MX1

    #print("Theoretical values found:")        
    #print(mn1, mn2,mn3)   

    return mn1, mn2, mn3

#run all dataframe
MatrixDiag_new=np.vectorize(MATRIX_NU_DIAG,excluded={'vev':246.2,'LAMBDA':1E16},doc='Input for pyfunc below:\
			    YB11,YB12,YB13,YB21,YB22,YB23,YA11,YA12,YA13,YA21,YA22,YA23,mS1,mS2')
