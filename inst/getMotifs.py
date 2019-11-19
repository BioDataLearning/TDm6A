#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on 09/06/2019

@author: Jun Wang
"""
def TDm6A_getMotifs(x_train,y_train, x_test, y_test,cellType):
    import numpy as np

    x_train = np.array(x_train)
    y_train = np.array(y_train)
    x_test = np.array(x_test)
    y_test = np.array(y_test)

    import keras
    from keras import backend as K
    from keras.models import Sequential
    from keras.layers import Conv1D, MaxPooling1D, LSTM, Bidirectional
    from keras.layers import Dense, Dropout, Flatten
    from keras import initializers
    from keras import optimizers
    from keras import losses

    np.random.seed(420)

    model = Sequential()

    conv_layer = Conv1D(filters=75,
                    kernel_size=6,
                    strides=1,
                    padding='valid',
                    activation='relu',
                    input_shape=(1001,4))

    model.add(conv_layer)
    model.add(Dropout(0))

    model.add(MaxPooling1D(pool_size=4, strides=4))
    model.add(Dropout(0.5))
    
    model.add(Bidirectional(LSTM(35, return_sequences = True)))
    model.add(Dropout(0.6))

    model.add(Flatten())
    model.add(Dense(2,activation='sigmoid'))

    model.summary()

    model.compile(loss=losses.binary_crossentropy,
              optimizer=optimizers.Adamax(),
              metrics=['accuracy'])

    model.fit(x_train, y_train,
          batch_size=128,
          epochs=60,
          validation_split=0.2)

    # motif visualization:
    ### motif visualization for BS model human:
    conv_output = conv_layer.get_output_at(0)
    f = K.function([model.input], [K.argmax(conv_output, axis=1), K.max(conv_output, axis=1)])

    motifs = np.zeros((75, 6, 4))
    nsites = np.zeros(75)

    # select the positive samples from test set:
    row = y_test[:,1]==1.0
    y_test_positive = y_test[np.ix_(row)]
    x_test_positive = x_test[np.ix_(row)]

    for i in range(0, len(x_test_positive), 100):
        x = x_test_positive[i:i+100]
        z = f([x])
        max_inds = z[0] # N x M matrix, where M is the number of motifs
        max_acts = z[1]
        for m in range(75):
            for n in range(len(x)):
                # Forward strand
                if max_acts[n, m] > 0:
                    nsites[m] += 1
                    motifs[m] += x[n, max_inds[n, m]:max_inds[n, m] + 6, :]

    print('Making motifs')
    file_name = cellType+'_Position_Probability_Matrix.txt'
    motifs_file = open(file_name, 'w')
    motifs_file.write('MEME version 4.9.0\n\n'
                  'ALPHABET= ACGU\n\n'
                  'strands: + -\n\n'
                  'Background letter frequencies (from uniform background):\n'
                  'A 0.25000 C 0.25000 G 0.25000 U 0.25000\n\n')

    for m in range(75):
        if nsites[m] == 0:
            continue
        motifs_file.write('MOTIF M_n%i O%i\n' % (m, m))
        motifs_file.write("letter-probability matrix: alength= 4 w= %i nsites= %i E= 1337.0e-6\n" % (6, nsites[m]))
        for j in range(6):
            motifs_file.write("%f %f %f %f\n" % tuple(1.0 * motifs[m, j, 0:4] / np.sum(motifs[m, j, 0:4])))
        motifs_file.write('\n')

    motifs_file.close()
    print("Done")

