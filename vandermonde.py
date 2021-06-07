from sage.all import matrix,SR

def Vandermonde_first_second_third(x_num,x_den,y_num,y_den):
	M=matrix(SR,8,8)
	x=SR(x_num)/x_den
	y=SR(y_num)/y_den
		
	for i in range(7):
		M[0,i]=x**(7-i)
		M[1,i]=y**(7-i)
	M[0,7]=SR(1)
	M[1,7]=SR(1)
	for i in range(6):
		M[2,i]=(7-i)*x**(6-i)
		M[3,i]=(7-i)*y**(6-i)
	M[2,6]=1
	M[3,6]=1
	for i in range(5):
		M[4,i]=(6-i)*(7-i)*x**(5-i)
		M[5,i]=(6-i)*(7-i)*y**(5-i)
	M[4,5]=2
	M[5,5]=2
	for i in range(4):
		M[6,i]=(5-i)*(6-i)*(7-i)*x**(4-i)
		M[7,i]=(5-i)*(6-i)*(7-i)*y**(4-i)	
	M[6,4]=6
	M[7,4]=6
	return M

def Vandermonde_first_second(x_num,x_den,y_num,y_den):
	M=matrix(SR,6,6)
	x=SR(x_num)/x_den
	y=SR(y_num)/y_den
		
	for i in range(5):
		M[0,i]=x**(5-i)
		M[1,i]=y**(5-i)
	M[0,5]=SR(1)
	M[1,5]=SR(1)
	for i in range(4):
		M[2,i]=(5-i)*x**(4-i)
		M[3,i]=(5-i)*y**(4-i)
	M[2,4]=1
	M[3,4]=1
	for i in range(3):
		M[4,i]=(5-i)*(4-i)*x**(3-i)
		M[5,i]=(5-i)*(4-i)*y**(3-i)
	M[4,3]=2
	M[5,3]=2
	return M
