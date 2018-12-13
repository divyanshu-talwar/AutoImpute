from __future__ import division, print_function, absolute_import
import os
import argparse
import datetime
import matplotlib
import numpy as np
import scipy.io
import tensorflow as tf
from sklearn.metrics import mean_absolute_error, mean_squared_error
	
matplotlib.use('Agg')

import matplotlib.pyplot as plt

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

parser = argparse.ArgumentParser(description = 'AutoImpute')

# Print debug statements
parser.add_argument('--debug', type = bool, default=True, nargs = '+', help = "Want debug statements.")
parser.add_argument('--debug_display_step', type=int, default=1, help="Display loss after.")

# Hyper-parameters
parser.add_argument('--hidden_units', type=int, default=2000, help="Size of hidden layer or latent space dimensions.")
parser.add_argument('--lambda_val', type=int, default=1, help="Regularization coefficient, to control the contribution of regularization term in the cost function.")
parser.add_argument('--initial_learning_rate', type=float, default=0.0001, help="Initial value of learning rate.")
parser.add_argument('--iterations', type=int, default=7000, help="Number of iterations to train the model for.")
parser.add_argument('--threshold', type=int, default=0.0001, help="To stop gradient descent after the change in loss function value in consecutive iterations is less than the threshold, implying convergence.")

# Data
parser.add_argument('--data', type = str, default='blakeley.csv', help = "Dataset to run the script on. In the paper we choose from : ['blakeley.csv', 'jurkat-293T.mat', 'kolodziejczyk.csv', 'PBMC.csv', 'preimplantation.mat', 'quake.csv', 'usoskin.csv', 'zeisel.csv']")

# Run the masked matrix recovery test
parser.add_argument('--masked_matrix_test', type = bool, default=False, nargs = '+', help = "Run the masked matrix recovery test?")
parser.add_argument('--masking_percentage', type = float, default=10, nargs = '+', help = "Percentage of masking required. Like 10, 20, 12.5 etc")

# Model save and restore options
parser.add_argument('--save_model_location', type=str, default='checkpoints/model1.ckpt', help="Location to save the learnt model")
parser.add_argument('--load_model_location', type=str, default='checkpoints/model0.ckpt', help="Load the saved model from.")
parser.add_argument('--log_file', type=str, default='log.txt', help="text file to save training logs")
parser.add_argument('--load_saved', type=bool, default=False, help="flag to indicate if a saved model will be loaded")

# masked and imputed matrix save location / name
parser.add_argument('--imputed_save', type=str, default='imputed_matrix', help="save the imputed matrix as")
parser.add_argument('--masked_save', type=str, default='masked_matrix', help="save the masked matrix as")

FLAGS = parser.parse_args()

if __name__ == '__main__':
	# started clock
	start_time = datetime.datetime.now()

	if not os.path.exists('checkpoints'):
		os.makedirs('checkpoints')

	if FLAGS.debug:
	    if not FLAGS.load_saved:
	        with open(FLAGS.log_file, 'w') as log:
	            log.write('Step\tLoss\tLoss per Cell\t Change \n')

	# reading dataset
	try:
		extn = FLAGS.data.split('.')[1]
		if(extn == 'mat'):
			print("[!data read] Reading from data/" + FLAGS.data)
			processed_count_matrix = scipy.io.mmread("data/" + FLAGS.data)
			processed_count_matrix = processed_count_matrix.toarray()
			processed_count_matrix = np.array(processed_count_matrix)
		else:
			print("[!data read] Reading from data/" + FLAGS.data)
			with open("data/" + FLAGS.data) as f:
			    ncols = len(f.readline().split(','))
			processed_count_matrix = np.loadtxt(open("data/" + FLAGS.data, "rb"), delimiter=",", skiprows=1, usecols=range(1,ncols+1))
	except :
		print("[!data read] Please make sure that your processed dataset is in data/ in .mat or .csv format and you have entered the filename as data parameter. e.g. blakeley.csv")
		exit()

	dataset = FLAGS.data.split('.')[0]

	if(FLAGS.masked_matrix_test):
		masking_percentage = FLAGS.masking_percentage
		masking_percentage = masking_percentage/100.0

		idxi, idxj = np.nonzero(processed_count_matrix)

		ix = np.random.choice(len(idxi), int(np.floor(masking_percentage * len(idxi))), replace = False)
		store_for_future = processed_count_matrix[idxi[ix], idxj[ix]]
		indices = idxi[ix], idxj[ix]

		processed_count_matrix[idxi[ix], idxj[ix]] = 0  # making masks 0
		matrix_mask = processed_count_matrix.copy()
		matrix_mask[matrix_mask.nonzero()] = 1 

		if(FLAGS.masked_save):
			scipy.io.savemat(FLAGS.masked_save + str(masking_percentage*100) + ".mat", mdict = {"arr" : processed_count_matrix})

		mae = []
		rmse = []
		nmse = []

	# finding number of genes and cells.
	genes = processed_count_matrix.shape[1]
	cells = processed_count_matrix.shape[0]
	print("[info] Genes : {0}, Cells : {1}".format(genes, cells))

	# placeholder definitions
	X = tf.placeholder("float32", [None, genes])
	mask = tf.placeholder("float32", [None, genes])

	matrix_mask = processed_count_matrix.copy()
	matrix_mask[matrix_mask.nonzero()] = 1

	print("[info] Hyper-parameters")
	print("\t Hidden Units : " + str(FLAGS.hidden_units))
	print("\t Lambda : {0}".format(FLAGS.lambda_val))
	print("\t Threshold : " + str(FLAGS.threshold))
	print("\t Iterations : " + str(FLAGS.iterations))
	print("\t Initial learning rate : " + str(FLAGS.initial_learning_rate))

	# model definition
	weights = {
		'encoder_h': tf.Variable(tf.random_normal([genes, FLAGS.hidden_units])),
		'decoder_h': tf.Variable(tf.random_normal([FLAGS.hidden_units, genes])),
		}
	biases = {
		'encoder_b': tf.Variable(tf.random_normal([FLAGS.hidden_units])),
		'decoder_b': tf.Variable(tf.random_normal([genes])),
	}

	def encoder(x):
		layer_1 = tf.nn.sigmoid(tf.add(tf.matmul(x, weights['encoder_h']), biases['encoder_b']))
		return layer_1

	def decoder(x):
		layer_1 = tf.add(tf.matmul(x, weights['decoder_h']), biases['decoder_b'])
		return layer_1

	encoder_op = encoder(X)
	decoder_op = decoder(encoder_op)

	# loss definition
	y_pred = decoder_op
	y_true = X 
	rmse_loss = tf.pow(tf.norm(y_true - y_pred * mask), 2)
	regularization = tf.multiply(tf.constant(FLAGS.lambda_val/2.0, dtype="float32"), tf.add(tf.pow(tf.norm(weights['decoder_h']), 2), tf.pow(tf.norm(weights['encoder_h']), 2)))
	loss = tf.add(tf.reduce_mean(rmse_loss), regularization)
	optimizer = tf.train.RMSPropOptimizer(FLAGS.initial_learning_rate).minimize(loss)

	init = tf.global_variables_initializer()

	saver = tf.train.Saver()

	with tf.Session() as sess:
		if(FLAGS.load_saved):
			saver.restore(sess, FLAGS.load_model_location)
			print("[info] model restored.")
		else:
			sess.run(init)
		prev_loss = 0
		for k in range(1, FLAGS.iterations+1):
			_, loss = sess.run([optimizer, rmse_loss], feed_dict={X: processed_count_matrix, mask: matrix_mask})
			lpentry = loss/cells
			change = abs(prev_loss - lpentry)
			if ( change <= FLAGS.threshold ):
				print("Reached the threshold value.")
				break
			prev_loss = lpentry
			if(FLAGS.debug):
				if (k - 1) % FLAGS.debug_display_step == 0:
					print('Step %i : Total loss: %f, Loss per Cell : %f, Change : %f' % (k, loss, lpentry, change))
					with open(FLAGS.log_file, 'a') as log:
						log.write('{0}\t{1}\t{2}\t{3}\n'.format(
						    k,
						    loss,
						    lpentry,
						    change
						))
			if((k-1) % 5 == 0):
				save_path = saver.save(sess, FLAGS.save_model_location)
		imputed_count_matrix = sess.run([y_pred], feed_dict={X: processed_count_matrix, mask: matrix_mask})
		scipy.io.savemat(FLAGS.imputed_save + ".mat", mdict = {"arr" : imputed_count_matrix})
		
		if(FLAGS.masked_matrix_test):	
			predictions = []

			for idx, value in enumerate(store_for_future):
				prediction = imputed_count_matrix[0][indices[0][idx], indices[1][idx]]
				predictions.append(prediction)

	if(FLAGS.masked_matrix_test):
		store_for_future = np.array(store_for_future)
		predictions = np.array(predictions)

		print("<------------------------------------Statistics for " + FLAGS.data.split('.')[0] + " dataset------------------------------------>")
		print("MAE = {0}".format( mean_absolute_error(store_for_future, predictions) ))
		print( "RMSE = {0}".format( mean_squared_error(store_for_future, predictions) ** 0.5 ))
		print("NMSE = {0}".format( (np.linalg.norm(store_for_future - predictions) / np.linalg.norm(store_for_future)) ))

	finish_time = datetime.datetime.now()
	print("[info] Total time taken = {0}".format(finish_time - start_time))