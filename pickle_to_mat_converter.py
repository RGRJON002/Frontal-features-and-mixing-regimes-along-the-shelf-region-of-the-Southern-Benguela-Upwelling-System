import numpy, scipy.io
import pickle, sys

source_name = 'SBUS_300m.pckl'
dest_name = 'SBUS_300m_MAT.mat'
	
a=pickle.load( open( source_name, "rb" ) )
scipy.io.savemat(dest_name, mdict={'pickle_data': a})
print("Data successfully converted to .mat file with variable name \"pickle_data\"")

	

