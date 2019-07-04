import random as rd
import numpy as np
from sys import getsizeof
from math import sqrt, floor, log2

class RSA:
	
	# Huge prime number
	_PRIME = 4879228431756730165299385010232837291120885737358336711467496639025522618199839237653918283371408014017658368720147197955013865183708637800199520868317779
	def __init__(self, n, k0, e):
		self.n = n
		self.k0 = k0
		self.e = e
		pass
	
	@staticmethod
	def getPrime(): return RSA._PRIME
	
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
	def SnM(a , b , n : int ): 
		if (type(a) == str):
			a = int(a, 2)
		if (type(b) == str):
			b = int(a, 2)
		return pow(a, b, n)
	
	@staticmethod
	def EA(a : int, b : int) -> int:


		# Cas de base
		if (b == 0): return a
		if (a == b): return b
		
		# Appel récursif selon le cas 
		if (a > b): return RSA.EA(b, a%b)
		else : return RSA.EA(a, b%a)
	
	@staticmethod 
	def Inverse( a : int, n : int):
		inv = RSA.EEA(a, n)[0]
		if inv < 0: inv += n
		return inv
	
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
			z = pow(a, _n, n)

			if (z != 1 and z != n-1): 
				# Boucle pour déterminer si 'n' est composé
				j = 1
				while j < u and z != n-1:
					# Nouvelle valeur de z
					z = pow(z, _n, n)
					if z == 1: return False
					# Augmente le compteur
					j += 1
				if z != n-1: return False
		return True
	
	#staticmethod
	def OAEP( msg : str, n : int , k0 : int) -> str:

		k1 = n - k0 - len(msg)
		# Genere un padding
		padding = "0" * k1
		# Ajoute le padding
		m = msg + padding
		
		# Génère un nombre aléatoire
		r = bin(rd.randrange(0,pow(2, k0-1)))[2:]
		r = r.zfill(k0)
		
		# Hash 'r' à n-k0 bits
		hash_r = RSA.hash_function(r, n-k0)
		
		# X == m XOR hash_r
		X = ''.join(str(int(a) ^ int(b)) for a,b in zip(m, hash_r))
		# Hash 'X'
		hash_X = RSA.hash_function(X, k0)

		# Y == r XOR hash_X
		Y = ''.join(str(int(a) ^ int(b)) for a,b in zip(r, hash_X))
		
		return X + Y
	
	#@staticmethod
	def OAEP_inv( XY : str, n : int, k0 : int) -> str:


		# Sépare le message en deux parties: X || Y
		X = XY[:n-k0]
		Y = XY[n-k0:]
		
		# Hash 'X'
		hash_X = RSA.hash_function(X, k0)#[2:].zfill(k0)
		# Retrouve 'r'
		r = ''.join(str(int(a)^int(b)) for a, b in zip(Y, hash_X))
		# Hash 'r'
		hash_r = RSA.hash_function(r, n-k0)#[2:].zfill(l-k0)
		# Retrouve 'msg + padding'
		m = ''.join(str(int(a)^int(b)) for a, b in zip(X, hash_r))
		
		return m
	
	#@staticmethod
	def genKeys(l : int):
		
		e = 217
		
		# Obtient les bonnes valeurs de q et p
		while True :
			p, q = e+1, e+1
			while (p-1)%e == 0: p = RSA.bigPrime(l//2) 
			while (q-1)%e == 0: q = RSA.bigPrime(l//2+1)
			
			phi_n = (p-1) * (q-1)
			
			if RSA.EA(e, phi_n) == 1 : break
		
		# Détermine la valeur de 'n' et 'phi_n'
		n = p * q
		phi_n = (p-1) * (q-1)
		
		#print(RSA.EEA(e, phi_n))
		
		# Détermine la valeur de 'd'
		d = RSA.Inverse(e, phi_n)
		#print("d*e:", (d*e)%phi_n)
		#d = RSA.SnM(e, phi_n-1, n)
		
		
		# Forme les clés
		PK = [n, e]
		SK = [p, q, d]
		
		return PK, SK

	@staticmethod
	def exp_CRT(C : str, SK : list):
		
		C_num = int(C, 2)
		
		p, q, d = SK
		N = p*q
		n = len(bin(N)[2:])-1
		d_p, d_q = d%(p-1), d%(q-1)
		k_p, k_q = RSA.Inverse(q, p), RSA.Inverse(p, q)
		
		x_p, x_q = C_num%p, C_num%q
		y_p, y_q = RSA.SnM(x_p, d_p, p), RSA.SnM(x_q, d_q, q)
		
		msg = bin(((q*k_p)*y_p + (p*k_q)*y_q)%N)[2:].zfill(n-1)
		msg_int = ((q*k_p)*y_p + (p*k_q)*y_q)%N
		
		
		return msg
		 
	
	
	#self, @staticmethod
	def encrypt(msg : str, PK : list):
		"""encrypt Méthode pour encrypter un message avec une clé publique
		
		Args:
			msg (str): Message à encrypter
			PK (list): Clé publique
			
		Returns:
			cipher (str): Message encrypté
		"""
		# Obtient 'n' et 'e' de la clé publique
		N, e = PK
		# Obtient la taille
		n = len(bin(N)[2:])
		_n = n - 1
		
		nBlock = len(msg)//_n
		print("reste :",len(msg)%_n)

		# Initialise le cipher
		msg_info = RSA.msg_info_block(msg, _n)
		msg_info_block = RSA.OAEP(msg_info, _n, k0)
		cipher = bin(RSA.SnM(msg_info_block, e, N))[2:].zfill(n)

		# Encrypte le message par block de n-bits
		for i in range(nBlock):

			block = msg[i*_n : (i+1)*_n]
			cipher += bin(RSA.SnM(block, e, N))[2:].zfill(n)
					
		f_block = ''

		if ((msg_info[0] == '0') and (len(msg)%_n != 0)) : 
			#un pad de zero comme m est très grand
			pad = "0"*(_n-(len(msg)%_n))
			f_block = msg[nBlock*_n :]+pad
			cipher += bin(RSA.SnM(f_block, e, N))[2:].zfill(n)
		
		if (len(msg)%_n != 0 and msg_info[0]=='1'):
			
			f_block = RSA.OAEP(msg[nBlock*_n :], _n, k0)
			cipher += bin(RSA.SnM(f_block, e, N))[2:].zfill(n)

		return cipher
	
	#@staticmethod
	def decrypt(cipher : str, SK : list):
		"""decrypt Méthode pour décrypter un cipher selon une clé privée
		
		Args:
			cipher (str): Meesage à décrypter
			SK (list): Clé privée
		
		Returns:
			msg (str): Message décrypté
		"""
		# Taille des blocks encryptés
		n = len(bin((SK[0] * SK[1]))[2:])
		# Détermine le nombre de block à décrypter
		nblock = len(cipher) // n
		# Initialise le message
		msg = ""

		# Décrypte le block d'information
		msg_info_D = RSA.exp_CRT(cipher[0:n], SK).zfill(n-1)
		#print(len(msg_info_D))
		msg_info = RSA.OAEP_inv(msg_info_D, n-1, k0)#.zfill(n-1)
		#print(len(msg_info),"msg_info:    ", msg_info)

		# Détermine si un padding a été fait
		padded = (msg_info[0] == "1")
		print("padded :",padded)
		
		"""
		# Obtient la taille du message
		index = 1;
		while msg_info[index] == '0': 
			index += 1
			if index == len(msg_info): break
				
		if index < len(msg_info):
			length_bin = msg_info[index:]
		else:
			length_bin = "0"
			
		print("length_bin:",length_bin)
		length = int(length_bin, 2)
		print("length:",length)
		"""	
		
		# Boucle pour décrypter les blocks
		for i in range(1, nblock-1):
			
			block = cipher[i*n: (i+1)*n]
			#print(i," cipher len :", len(cipher[i*n: (i+1)*n]), "block len : ",len(RSA.exp_CRT(block, SK)))
			msg += RSA.exp_CRT(block, SK).zfill(n-1)
			
		
		block = RSA.exp_CRT(cipher[-n:], SK).zfill(n-1)
				
		if padded:
			block = RSA.OAEP_inv(block, n-1, k0)
		
		msg+=block
		
		return msg


	def msg_info_block( msg : str, _n : int):
		#BLOCK débutant le message, toujours OAEP-padded. Le premier bit indique si
		#le dernier block du message est OAEP-padded (BOOL) ou non, suivi de la 
		#longeur du message.
		
		#taille du message 
		length = len(msg)
		
		# Détermine la taille du restant
		reminder_size = len(msg)%_n
		
		#ajoute le OAEP-padded BOOL
		if reminder_size == 0: msg_info = '0'
		elif (reminder_size <= _n-k0): msg_info = '1'	
		else : msg_info = '0'
		
		msg_info+= bin(length)[2:]		
		return msg_info
"""		
	@staticmethod
	def sign(msr : str, SK : list):
		hash_m = RSA.hash_function(int(msr, 2), 160)
		return RSA.exp_CRT(hash_m, SK ).zfill(SK[0]*SK[1])
	
	@staticmethod
	def verify_sign(msr : str , sign_m : str, PK : list):
		hash_m1 = hash_m = RSA.hash_function(int(msr, 2), 160)
		hash_m2 = bin(RSA.SnM(sign_m, PK[1], PK[0]))[2:].zfill(PK[0])
		if hash_m1 == hash_m2 : return True
		else : return False
		
"""		

"""
SIGNATURE : encrypt(sign(hash(C),timestamp) ; C) avec sign qui est juste exp_CRT()
et (on pourait la mettre dans le bloc initial) et alors il faut juste aller la chercher, 
et la vérifier à la fin de decrypt()
"""
"""

k0 = 128 #l'ai défini à l'exterieur

PK, SK = RSA.genKeys(512)
n = len(bin(PK[0])[2:])
print("n :", n)

m = rd.randrange(0, pow(2, 120))
M = bin(m)[2:]

C = RSA.encrypt(M, PK)
D = RSA.decrypt(C, SK)
print(len(M)," M :",M[:50])
print(len(D)," D :",D[:50])

XOR = ''.join(str(int(a) ^ int(b)) for a,b in zip(M, D))

print("ZERO IS SUCESS : ",int(XOR))

"""
"""
M= '100001111011010100111'
S = RSA.sign(M,SK)
V = RSA.verify_sign(M , S , PK )

print(V)
"""
PK, SK = RSA.genKeys(1024)
N = PK[0]
print(N)
d=5
m = floor(pow(N, 1/d))
print("N :",N, "; n :", len(bin(N)[2:]), "; d =", d, "; m :",m , "m^d :", pow(m,d), "m^d+1 :")
print(pow(m, d+1))

