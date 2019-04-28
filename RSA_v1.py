
# coding: utf-8

# In[ ]:


import random as rd
import numpy as np
from sys import getsizeof
from math import sqrt, floor, log2


# Implémentation de RSA sous forme de classe.
# 

# In[175]:


class RSA:
    
    # Huge prime number
    _PRIME = 4879228431756730165299385010232837291120885737358336711467496639025522618199839237653918283371408014017658368720147197955013865183708637800199520868317779
    
    def __init__(self):
        pass
    @staticmethod
    def getPrime(): return RSA._PRIME
    
    @staticmethod
    def SnM(a , b , n : int ): 
        if type(a) == str : a = int(a, 2)
        if type(b) == str : b = int(b, 2)
            
        return pow(a, b, n)
        
    
    
    @staticmethod
    def bigPrime(l : int) -> int:


        # Génère un nombre aléatoire de l-bits (impair)
        n = rd.getrandbits(l) | 1
        n |= (1 << l - 1)

        # Boucle pour s'assurer que n est premier
        while not RSA.millerRabin(n, 5):
            # Génère un nombre aléatoire de l-bit et fait un OR avec 1
            n = rd.getrandbits(l) | 1   # Obligatoirement impair
            n |= (1 << l - 1)
        return n
    
    @staticmethod
    def EA(a : int, b : int) -> int:


        # Cas de base
        if (b == 0): return a
        if (a == b): return b
        
        # Appel récursif selon le cas 
        if (a > b): return RSA.EA(b, a%b)
        else : return RSA.EA(a, b%a)
    
    @staticmethod
    def hash_function(m : int, l : int):
        """hash_function Retourne un hash de longueur 'l' de 'm'
        
        Args:
            m (int): Message à hasher
            l (int): Longueur du hash 
        """
        if (type(m) == str):
            m = int(m, 2)
        # Convertie le message en nombre
        mhash = bin((m * RSA.getPrime()) % pow(2, l-1))
        return mhash[2:].zfill(l)

    @staticmethod
    def EEA(a : int, b : int) -> int:

        # S'assure que le nombre 'a' est plus grand
        if (a < b): RSA.EEA(b, a)

        # Initialisation des coefficients := [[x0, x1], [y0, y1]]
        coeff = [[1, 0], [0, 1]]

        # Initialise un tableau pour contenir les coefficients 'x' et 'y' ainsi que le PGCD
        result = [None, None, None]

        # Boucle principale
        while(True) :
            # Détermine q et r
            r = a % b   
            if (r == 0): return result
            q = a // b 

            # Obtient les résultats de la ie itération
            result[0] = coeff[0][0] - q * coeff[0][1]
            result[1] = coeff[1][0] - q * coeff[1][1]
            result[2] = r

            # Modifie les coefficients 
            coeff[0] = [coeff[0][1], result[0]]
            coeff[1] = [coeff[1][1], result[1]]

            # Modifie les valerus de 'a' et 'b' qui seront utilisées pour la prochaine itération
            a = b
            b = r

    @staticmethod
    def millerRabin(n : int, s : int) -> bool:

        # Vérifie si le nombre 'p' est pair; retourne vrai si 'p' vaut 2 et faux sinon
        if (n % 2 == 0) : return n == 2
        if n == 3: return True

        # Initialise r, u et _p
        u, _n = 0, n-1
        # Détermine la valeur de 'u'
        while _n & 1 == 0 :
            u += 1      # Augmente la valeur de 'u'
            _n //= 2    # Division entière de _n par 2

        # Boucle principale
        for i in range(s):
            # Détermine un a aléatoire [2, n-1[
            a = rd.randrange(2, n-1)
            # Détermine 'z' initiale
            z = RSA.SnM(a, _n, n)

            if (z != 1 and z != n-1): 
                # Boucle pour déterminer si 'n' est composé
                j = 1
                while j < u and z != n-1:
                    # Nouvelle valeur de z
                    z = RSA.SnM(z, _n, n)
                    if z == 1: return False
                    # Augmente le compteur
                    j += 1
                if z != n-1: return False
        return True
    
    @staticmethod
    def OAEP(msg : str, k0 : int, l : int):

        k1 = l - k0 - len(msg)
        # Genere un padding
        padding = "0" * k1
        # Ajoute le padding
        m = msg + padding
        
        # Génère un nombre aléatoire
        r = bin(rd.randrange(0, pow(2, k0-1)))[2:]
        r = r.zfill(k0)
        
        # Hash 'r' à n-k0 bits
        hash_r = RSA.hash_function(r, l-k0)[2:].zfill(l-k0)
        # Hash 'm' à k0 bits
        hash_m = RSA.hash_function(m, k0)[2:].zfill(k0)
        
        # X == m XOR hash_r
        X = ''.join(str(int(a) ^ int(b)) for a,b in zip(m, hash_r))
        # Y == hash_m XOR r
        Y = ''.join(str(int(a) ^ int(b)) for a,b in zip(hash_m, r))
        return X + Y
    
    
    @staticmethod
    def genKeys(l : int): #l longeur voulu pour N

        e = 257
        q, p = e + 1, e + 1
    
        while (p-1)%e == 0 : p = RSA.bigPrime(l//2) 
        while (q-1)%e == 0 : q = RSA.bigPrime(l//2)

        n = q*p
        phi_n = (q-1)*(p-1)
        d = RSA.SnM(e, phi_n-1, n)
        
        PK = [n, e]
        SK = [p, q, d]
        
        return PK, SK

    
    @staticmethod
    def exp_CRT(C : str, SK : list):
        
        C_num = int(C, 2)
        #print("C_num :",C_num)
        
        p, q, d = SK
        N = p*q
        n = len(bin(N)[2:])-1
        d_p, d_q = d%(p-1), d%(q-1)
        k_p, k_q = RSA.EEA(q, p)[0], RSA.EEA(p, q)[0]
        
        x_p, x_q = C_num%p, C_num%q
        y_p, y_q = RSA.SnM(x_p, d_p, p), RSA.SnM(x_q, d_q, q)
        
        msg = bin(((q*k_p)*y_p + (p*k_q)*y_q)%N)[2:].zfill(n)
        msg_int = ((q*k_p)*y_p + (p*k_q)*y_q)%N
        
        #print("MSG int :",msg_int,"len MSG :",len(msg))
        
        return msg #un bitstring
    
   # @staticmethod
   # def block_encrypt( block : str, PK : list):
   #     return RSA.SnM(block, PK[1], PK[0])
    
    @staticmethod
    def encrypt(msg : str , PK : list): #msg un bitstring
        N, e = PK
        n = len(bin(N)[2:])-1
        cipher = ''
        #msg_bit = ''.join(map(msg,bytearray(msg,'utf8'))) #si on veut mettre un vrai string..

        for i in range(len(msg)//n) :
            
            block = msg[i*n :(i+1)*n]
            #print("block :",block,"  ",len(block))
            cipher+=(bin(RSA.SnM(block, PK[1], PK[0]))[2:].zfill(n))
            
            #print("cipher block :",cipher,"  ",len(bin(RSA.SnM(block, PK[1], PK[0]))[2:].zfill(n)))
        
        reminder = len(msg)%n #je sais pas trop quoi faire pour reminder très près de n... il faut un certain lousse pour padder non? ou on fait juste changer la longeur de output du hash... en réalité ils utilise sha-x avec un output fixe...
        
        return cipher #un bitsrting
    
    @staticmethod
    def decrypt(C : str, SK : list):
        n = len(bin(SK[0]*SK[1])[2:])
        #print("D n :",n)
        nblocks = len(C)//n
        msg = ''
        #print(len(C))
        #print(nblocks)
        for i in range(nblocks) :
            block = C[i*n :(i+1)*n]
            #print("block :",block)
            msg+=(RSA.exp_CRT(block, SK))
            #print(msg)
       
        #reminder OAEP
        return msg #un bitstring...

    def Sign( M : str, SK : list, hash_fct): #pour être tight? avec un time stamp
        pass
        
    def checkSign( S : str, PK : list, hash_fct ):
        
        pass
        
    def MAC( M : str, Bob_PK : list, A_SK : list, hash_fct ):
        pass




"""
msg = bin(rd.randrange(0, pow(2, 520))+ pow(2, 520))[2:] 
print(len(msg))

A_PK, A_SK = RSA.genKeys(512)
print(len(bin(A_PK[0])[2:]))
Cipher = RSA.encrypt(msg, A_PK)
#print(len(Cipher))

MSG = RSA.decrypt(Cipher, A_SK)
#print("MSG :",MSG) #msg n,a pas la bonne taille pour être decrypt sans OAEP.. TBC

#print(msg)
print("msg :",msg[:10])
print("MSG :",MSG[:10])

#print(''.join(str(int(a) ^ int(b)) for a,b in zip(msg[:(len(msg)-len(MSG))], MSG)))

 """       

PK = [91, 5]
SK = [7, 13, 29]

n = len(bin(PK[0])[2:])
print("n :", n)
M = bin(10)[2:].zfill(n-1)
print("M :",M)

C = RSA.encrypt(M, PK)


#print("cipher :", C)
#print("R_cipher :", bin(pow(10,5,91))[2:].zfill(n-1))
print("CRT :",RSA.exp_CRT(C, SK))
print("Decrypt :", RSA.decrypt(C, SK))

