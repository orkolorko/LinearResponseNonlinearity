from scipy.sparse import *
from scipy import io
from os import mkdir
import time
import numpy as np
#import joblib as jl
import gc
from tempfile import mkdtemp
import os.path as path

class huge_matrix:
	r""" Class representing huge matrices by storing on the hard drive
	horizontal blocks which are sparse matrices, works best for matrices
	of row_size> 65536
	Please note that the matrix are stored on disc in the COO format,
	the access is granted by converting them to the lil format.
	Everything is a big bloat but since the number of conversion operations
	should be relatively small, we hope it is not going to be too slow.
	"""
	
	
	def block_filename(self,i):
		"""Given a block number, returns the associated filename"""
		return str.format("{}/slice_{}_{}.npy",self._directory_name,i,self._timestamp)
	
	def __init__(self,size_rows,size_columns,split_rows=131072, directory_name=None, storage_method='joblib'):
		self._number_of_blocks=size_rows//split_rows
		print "number of blocks", self._number_of_blocks
		self._remaining_lines=size_rows%split_rows
		print "remaining lines", self._remaining_lines
		self._current_block=0
		self._directory_name=None
		self._split_rows=split_rows
		self._size_rows=size_rows
		self._size_columns=size_columns
		self._timestamp = int(time.time())
		self._storage_method=storage_method
		self._P=None
		
		if (directory_name is None):
			self._directory_name=str.format("Huge_matrix_{}x{}_{}_{}",size_rows,size_columns,split_rows,self._timestamp)
		else:
			self._directory_name=directory_name
		
		mkdir(self._directory_name)
		
		for i in range(self._number_of_blocks):
			self._P=coo_matrix((self._split_rows,self._size_columns))
			file_name=self.block_filename(i)
			if (self._storage_method=='joblib'):
				jl.dump(self._P,file_name)
			else:
				io.mmwrite(file_name,self._P)
	#		print "block ", i , "row, columns ", (self._P).shape
		if (self._remaining_lines != 0):
			self._P=coo_matrix((self._remaining_lines,size_columns))
			file_name=self.block_filename(self._number_of_blocks)
			if (self._storage_method=='joblib'):
				jl.dump(self._P,file_name)
			else:
				io.mmwrite(file_name,self._P)
			
		self._current_block=0
		file_name=str.format("{}/slice_{}_{}.npy",self._directory_name,0,self._timestamp)
		
		if (self._storage_method=='joblib'):
			self._P=(jl.load(file_name)).tolil()
		else:
			self._P=(io.mmread(file_name)).tolil()	
	
	def access_current_block(self):
		r""" This function returns the current block and its index"""
	
		return self._P,self._current_block
	
	def access_block_as_csr(self,i):
		r"""Given a block index this functions returns the block with this index as a csr matrix 
		good for matrix column multiplication"""
	
		self._current_block=i
		file_name=self.block_filename(i)
		
		if (self._storage_method=='joblib'):
			Q=(jl.load(file_name)).tocsr()
		else:
			Q=(io.mmread(file_name)).tocsr()
		return Q
	
	
	def access_block_as_lil(self,i):
		r"""Given a block index this functions returns the block with this index as a lil matrix 
		good for accessing and modyfying the coefficients"""
		self._current_block=i
		file_name=self.block_filename(i)
		if (self._storage_method=='joblib'):
			Q=(jl.load(file_name)).tolil()
		else:
			Q=(io.mmread(file_name)).tolil()
		return Q
	
	
	def write_current_block(self):
		r"""Writes the current block to disk"""
		file_name=self.block_filename(self._current_block)
		if (self._storage_method=='joblib'):
			jl.dump(self._P,file_name)
		else:
			io.mmwrite(file_name,self._P)
				
	
	def line_to_block_and_row_in_block(self,i):
		r"""
		Converts a line to the block index and its index inside the block
		"""
		block=i//self._split_rows
		row_in_block=i%self._split_rows
		return block,row_in_block
	
	def read_coefficients(self,i,j):
		r"""
		This function reads a coefficient from the matrix
		"""
		block,row_in_block=self.line_to_block_and_row_in_block(i)
		if (self._current_block==block):
			return self._P[row_in_block,j]
		else:
			# Writing the former current block on disk
			self.write_current_block()
			# Loading the new current block in memory, as a lil matrix
			self._P=self.access_block_as_lil(block)
			return (self._P[row_in_block,j]).todense()
	
	def write_coefficients(self,i,j,value):
		r""" This function writes a coefficient into the matrix. 
		If the block where we are writing is different from the current block
		it changes the block.
		Therefore it is better for performance to write as many coefficients
		in the same block as possible.
		This is going to be substituted by a buffer in memory that writes
		the coefficients in bigger chunks"""
		block,row_in_block=self.line_to_block_and_row_in_block(i)
		if (self._current_block==block):
			self._P[row_in_block,j]=value
		else:
			# Writing the former current block on disk
			self.write_current_block()
			# Loading the new current block in memory
			self._P=self.access_block_as_lil(block)
			self._P[row_in_block,j]=value
	
	def sum_coefficients(self,i,j,value):
		r""" This function writes a coefficient into the matrix,
		summing it to the current value. 
		If the block where we are writing is different from the current block
		it changes the block.
		Therefore it is better for performance to write as many coefficients
		in the same block as possible.
		This is going to be substituted by a buffer in memory that writes
		the coefficients in bigger chunks"""
		block,row_in_block=self.line_to_block_and_row_in_block(i)
		if (self._current_block==block):
			self._P[row_in_block,j]+=value
		else:
			# Writing the former current block on disk
			self.write_current_block()
			# Loading the new current block in memory
			self._P=self.access_block_as_lil(block)
			self._P[row_in_block,j]+=value
	
	def finished_assembling(self):
		r""" Call this function when the assembly is finished to insure
		that all the blocks are written on disc"""
		self.write_current_block()
	
	
	def __mul__(self,w):
		r"""
		Overloading the multiplication operator
		"""
		v=np.zeros(self._size_rows)
		for i in range(self._number_of_blocks):
			#print "product with block", i
			Q=self.access_block_as_csr(i)
			v[i*self._split_rows:(i+1)*self._split_rows]=Q*w
			del Q
			
		if (self._remaining_lines != 0):
			Q=self.access_block_as_csr(self._number_of_blocks)
			v[self._number_of_blocks*self._split_rows:self._number_of_blocks*self._split_rows+self._remaining_lines]=Q*w
			del Q
		
		#gc.collect()
		
		return v

	def convert_to_memmapped_csr(self):
		
		len_data=np.zeros(self._number_of_blocks+1)
		len_indices=np.zeros(self._number_of_blocks+1)
		len_indptr=np.zeros(self._number_of_blocks+1)
		
		for i in range(self._number_of_blocks):
			#print "product with block", i
			Q=self.access_block_as_csr(i)
			len_data[i]=len(Q.data)
			len_indices[i]=len(Q.indices)
			len_indptr[i]=len(Q.indptr)
			del Q
			
		if (self._remaining_lines != 0):
			Q=self.access_block_as_csr(self._number_of_blocks)
			len_data[self._number_of_blocks]=len(Q.data)
			len_indices[self._number_of_blocks]=len(Q.indices)
			len_indptr[self._number_of_blocks]=len(Q.indptr)
			
			del Q
		
		total_len_data=int(np.linalg.norm(len_data,1))
		total_len_indices=int(np.linalg.norm(len_indices,1))
		
		filename_data = path.join(mkdtemp(), 'data.dat')
		filename_indices = path.join(mkdtemp(), 'indices.dat')
		filename_indptr = path.join(mkdtemp(), 'indptr.dat')

		fp_data = np.memmap(filename_data, dtype='float64', mode='w+', shape=total_len_data)
		fp_indices = np.memmap(filename_indices, dtype='uint64', mode='w+', shape=total_len_indices)
		fp_indptr = np.memmap(filename_indptr, dtype='uint64', mode='w+', shape=self._size_rows+1)
		
		nnz_preceeding_bloc=0
		accumulated_nnz=0
		
		for i in range(self._number_of_blocks):
			#print "product with block", i
			Q=self.access_block_as_csr(i)
			
			nnz_current_block=Q.indptr[self._split_rows]
			print "nnz_current_bloc", nnz_current_block
			print "len data", len(Q.data)
			
			fp_data[accumulated_nnz:(accumulated_nnz+nnz_current_block)]=Q.data[:]
			fp_indices[accumulated_nnz:(accumulated_nnz+nnz_current_block)]=Q.indices[:]
			fp_indptr[i*self._split_rows:(i+1)*self._split_rows]=Q.indptr[0:self._split_rows]+accumulated_nnz
			
			accumulated_nnz+=nnz_current_block
			print "accumulated_nnz", accumulated_nnz
			
			del Q
			
		if (self._remaining_lines != 0):
			Q=self.access_block_as_csr(self._number_of_blocks)
			
			nnz_current_block=Q.indptr[self._remaining_rows]
			print "nnz_current_bloc (remaining rows)", nnz_current_block
			print "len data", len(Q.data)
			
			fp_data[accumulated_nnz:(accumulated_nnz+nnz_current_block)]=Q.data[:]
			fp_indices[accumulated_nnz:(accumulated_nnz+nnz_current_block)]=Q.indices[:]
			fp_indptr[(self._number_of_blocks-1)*self._split_rows:self_size_rows]=Q.indptr[0:self._remaining_rows]+accumulated_nnz
			
			del Q
		
		fp_indptr[self._size_rows]=accumulated_nnz
		
		P=csr_matrix((fp_data, fp_indices, fp_indptr),shape=(self._size_rows,self._size_columns))
		
		return P
		
		
def test_conversion_to_csr(P,Q):
	size=P._size_rows
	for i in range(10):
		v=np.random.uniform(-1,1,size)
		w=P*v-Q*v
		print np.linalg.norm(w,1)
	
	return 0	
	
		

def check_matrices(P_coarse,P_fine):
	a=[]
	n=P_coarse.ncols()
	print n
	for i in range(0,n):
		for j in range(0,n):
			if (P_fine.read_coefficients(i,j)-P_coarse[i,j])>0.0000001: 
				print 'Caraio'
				print P_fine.read_coefficients(i,j)-P_coarse[i,j]
				a.append((i,j))
	return a

